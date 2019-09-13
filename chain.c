/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2019 Jeremy Wang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "klib/ksort.h"
#include "klib/khash.h"
#include "klib/kvec.h"
#include "chain.h"

#define score_pos_gt(a,b) ((a).score > (b).score)
KSORT_INIT(score_pos_cmp, score_pos, score_pos_gt)

#define pos_pair_lt(a,b) ((a).tpos < (b).tpos)
KSORT_INIT(pos_pair_cmp, posPair, pos_pair_lt)

chain* do_chain(khash_t(matchHash) *hits, khash_t(matchHash) *rv_hits, int max_chains, int match_score, int max_gap, int min_chain_length) {
  khash_t(matchHash) **ht = malloc(sizeof(khash_t(matchHash)*) * 2);
  ht[0] = hits;
  ht[1] = rv_hits;

  // assess hits for each target
  pairVec anchors;
  khint_t bin;
  uint32_t target;

  scoreVec scores;
  kv_init(scores);
  score_pos s;

  int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
  int h = 50; // number of previous anchors to check
  int ref_offset = 0; // offset into the scores vector of the current target - sets the backward limit for finding chained anchors

  // iterate through hits for each target, and append them to the same scores vector
  int rv;
  for(rv = 0; rv <= 1; rv++) {
    for (bin = kh_begin(ht[rv]); bin != kh_end(ht[rv]); ++bin) {
      if (!kh_exist(ht[rv], bin)) continue;
      target = kh_key(ht[rv], bin);
      anchors = kh_val(ht[rv], bin);
      //fprintf(stderr, "Ref %u has %u anchors\n", target, kv_size(anchors));

      //sort anchor pairs by target pos increasing
      ks_mergesort(pos_pair_cmp, kv_size(anchors), anchors.a, 0);

      s.score = match_score;
      s.anchor_idx = 0;
      s.score_idx = ref_offset;
      s.prev = -1;
      s.used = 0;
      s.ref = target;
      s.rv = rv;
      kv_push(score_pos, scores, s); // this actually copies the score_pos struct values, so we can reuse 's'

      int st = 0;

      for(i = 1; i < kv_size(anchors); i++) {

        //s.score = match_score;
        s.score = 1;
        s.anchor_idx = i;
        s.score_idx = i + ref_offset;
        s.prev = -1; // no predecessor

        /*
        while(st < i && kv_A(anchors, i).tpos - kv_A(anchors, st).tpos > 100000) st++;
        s.score = i - st;
        */

        for(j = i < h ? 0 : i-h; j < i; j++) {
          qdiff = (int)kv_A(anchors, i).qpos - (int)kv_A(anchors, j).qpos;
          tdiff = (int)kv_A(anchors, i).tpos - (int)kv_A(anchors, j).tpos;
          if(tdiff <= 0 || qdiff <= 0 || qdiff > max_gap || tdiff > max_gap) { // we often have multiple hits to the same target pos, so we can't let them chain with either other
            continue;
          }
          /*
          diffdiff = (qdiff > tdiff ? qdiff - tdiff : tdiff - qdiff);
          gap_cost = (diffdiff == 0 ? 0 : 0.01 * match_score * diffdiff + 0.5 * log2(diffdiff)); // affine gap penalty a la minimap2
          qdiff = (qdiff > tdiff ? tdiff : qdiff); // now not qdiff but the minimum difference
          score = kv_A(scores, j + ref_offset).score + (qdiff > match_score ? match_score : qdiff) - gap_cost;
          */
          score = kv_A(scores, j + ref_offset).score + 1; // scores are now just the # of properly ordered anchors
          if(score > s.score) {
            s.score = score;
            s.prev = j + ref_offset;
          }
        }

        kv_push(score_pos, scores, s);
      }
      ref_offset = kv_size(scores);
    }
  }

  //fprintf(stderr, "%u scores\n", kv_size(scores));

  // copy unsorted scores:
  score_pos* anchor_scores = malloc(sizeof(score_pos) * kv_size(scores));
  memcpy(anchor_scores, scores.a, sizeof(score_pos) * kv_size(scores));
  //sort scores decreasing
  ks_mergesort(score_pos_cmp, kv_size(scores), scores.a, 0);

  // build non-overlapping chains from highest to lowest score
  chain* chains = calloc(max_chains, sizeof(chain));
  int c; // chain index
  i = 0; // scores index
  for(c = 0; c < max_chains && i < kv_size(scores); c++) {
    //fprintf(stderr, "building chain %d\n", c);
    int chain_len = 0;
    int chain_pos = kv_A(scores, i).score_idx; // get the pos which should be index into the anchor_scores array
    //fprintf(stderr, "chain pos: %d (of %u)\n", chain_pos, kv_size(scores));
    // get anchors associated with this score/ref
    bin = kh_get(matchHash, ht[kv_A(scores, i).rv], kv_A(scores, i).ref);
    if(bin == kh_end(ht[kv_A(scores, i).rv])) { // key not found, *shouldn't* happen
      fprintf(stderr, "something went very wrong with the chain computation!");
      return (chain*)NULL;
    }
    anchors = kh_val(ht[kv_A(scores, i).rv], bin);
    //fprintf(stderr, "pos %d used: %u\n", chain_pos, anchor_scores[chain_pos].used);
    while(anchor_scores[chain_pos].used == 0) {
      //fprintf(stderr, "backtracking from %d (%d:%d) to %d\n", chain_pos, kv_A(anchors, anchor_scores[chain_pos].anchor_idx).qpos, kv_A(anchors, anchor_scores[chain_pos].anchor_idx).tpos, anchor_scores[chain_pos].prev);
      chain_len++;
      if(anchor_scores[chain_pos].prev == -1)
        break;
      chain_pos = anchor_scores[chain_pos].prev;
    }
    if(chain_len < min_chain_length) {
      c--;
      i++;
      continue;
    }
    kv_init(chains[c].anchors);
    chains[c].ref = kv_A(scores, i).ref;
    chains[c].rv = kv_A(scores, i).rv;
    chains[c].score = kv_A(scores, i).score;
    fprintf(stderr, "chain %d (%c) (score %d, %d anchors) is from ref %u anchor %d\n", c, "+-"[kv_A(scores, i).rv], kv_A(scores, i).score, chain_len, chains[c].ref, kv_A(scores, i).anchor_idx);
    kv_resize(posPair, chains[c].anchors, chain_len);
    chains[c].anchors.n = chain_len; // set the size explicitly, then we'll set the values explicitly
    fprintf(stderr, "creating chain %d of length %d\n", c, kv_size(chains[c].anchors));
    chain_pos = kv_A(scores, i).score_idx;
    for(j = 0; j < chain_len; j++) {
      anchor_scores[chain_pos].used = 1;
      kv_A(chains[c].anchors, chain_len-1-j) = kv_A(anchors, anchor_scores[chain_pos].anchor_idx);
      chain_pos = anchor_scores[chain_pos].prev;
    }
    i++;
  }
  fprintf(stderr, "made %d chains\n", c);
  if(c < max_chains) kv_init(chains[c].anchors); // an empty vector will indicate the end of the chains array if there are fewer than max_chains results

  free(anchor_scores);
  kv_destroy(scores);
  return chains;
}

void free_chains(chain* chs) {
  int i;
  for(i = 0; ; i++) {
    if(kv_size(chs[i].anchors) > 0) {
      break;
    }
    kv_destroy(chs[i].anchors);
  }
  free(chs);
}

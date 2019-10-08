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

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "klib/kvec.h" // C dynamic vector
#include "klib/khash.h"
#include <zlib.h>
#include <math.h>
#include "chain.h"

#ifndef _kseq_
#define _kseq_

#include "klib/kseq.h"

// init kseq struct
KSEQ_INIT(gzFile, gzread)

#endif

#define score_gap(a) (a == 0 ? 0 : ((GAP_OPEN_1 + (a-1) * GAP_EXTEND_1 > GAP_OPEN_2 + (a-1) * GAP_EXTEND_2) ? GAP_OPEN_1 + (a-1) * GAP_EXTEND_1 : GAP_OPEN_2 + (a-1) * GAP_EXTEND_2))

typedef struct {
  char* s;
  uint32_t l;
} fa_seq;

// map seq name to index into *seq vec
KHASH_MAP_INIT_STR(faHash, uint32_t)

void version() {
  printf("Simple dynamic programming aligner\n");
}

void usage() {
  printf("Usage: aln [options]\n");
  printf("Commands:\n");
  printf("  global: end-to-end query and ref - useful for full-length amplicons, etc.\n");
  printf("  read-local: end-to-end ref, but query may be partial\n");
  printf("  ref-local: end-to-end reads, but ref may be partial - this is most useful for read -> reference mapping\n");
  printf("  local: both may be partial - often reasonable for long-read -> reference mapping\n");
  printf("Options:\n");
  printf("  -q: FASTA/Q[.gz] file with reads\n");
  printf("  -r: Reference FASTA/Q[.gz]\n");
  printf("  -v, --verbose: verbose\n");
  printf("  -h, --help: show this\n");
  printf("  -l, --limit N: stop processing after N query seqs\n");
  printf("Advanced options:\n");
  printf("  -b, --bandwidth: DP band width\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "verbose",                no_argument,       0, 'v' },
  { "help",                   no_argument,       0, 'h' },
  { "limit",                  required_argument, 0, 'l' },
  { "bandwidth",              required_argument, 0, 'b' },
  { 0, 0, 0, 0}
};

typedef kvec_t(unsigned char) charvec;

typedef struct {
  uint32_t ref_id;
  uint32_t pos;
} ref_pos;

typedef kvec_t(ref_pos) rpvec;

// string (k-mer) to ref/pos
KHASH_MAP_INIT_STR(kmerHash, rpvec)

typedef struct aln_result {
  int score;
  int qstart;
  int qend;
  int tstart;
  int tend;
  int failed; // boolean flag
} result;

unsigned char MATCH = 0, INS = 1, DEL = 2, MISMATCH = 3;
int MATCH_SCORE = 2;
int MISMATCH_SCORE = -4;
int GAP_OPEN_1 = -4;
int GAP_EXTEND_1 = -2;
int GAP_OPEN_2 = -24;
int GAP_EXTEND_2 = -1;
int LOW = -2000000000; // almost the lowest 32-bit integer
char compl[256];

float cigar_accuracy(unsigned char *cigar, int cigar_len) {
  int matches = 0, total = 0;
  int i;
  for(i = 0; i < cigar_len; i++) {
    if(cigar[i] == MATCH) {
      matches++;
    }
    total++;
  }
  return ((float)matches) / total;
}

char* rc(char* s, int l, char compl[256]) {
  char* r = malloc((l+1) * sizeof(char));
  int i;
  for(i = 0; i < l; i++)
    r[i] = compl[s[l-1-i]];
  r[i] = '\0';
  return r;
}

typedef struct {
  int score;
  unsigned char direction;
  uint16_t gap_size; // supports gaps up to 65Kbp
} dp_cell;

// absolute x,y offset of each band
typedef struct {
  uint32_t x;
  uint32_t y;
} os;

// returns a status code, 0 if OK, 1 if the alignment didn't make sense
int output_human_aln(char* qseq, int ql, char* tseq, int rv, charvec path, FILE* o) {
  int q = 0;
  int t = 0;
  uint32_t qst = 0;
  uint32_t tst = 0;
  char* qs = malloc(151 * sizeof(char));
  char* as = malloc(151 * sizeof(char));
  char* ts = malloc(151 * sizeof(char));
  qs[150] = '\0';
  as[150] = '\0';
  ts[150] = '\0';
  int i;
  //fprintf(o, "path of length %u\n", kv_size(path));
  for(i = 0; i < kv_size(path); i++) {
    if(kv_A(path, i) == MATCH) {
      qs[i%150] = rv ? compl[qseq[ql-1-q++]] : qseq[q++];
      as[i%150] = '|';
      ts[i%150] = tseq[t++];
      if(qs[i%150] != ts[i%150]) {
        fprintf(stderr, "ERROR: bad alignment at position %d in the path\n", i);
        qs[i%150+1] = '\0';
        as[i%150+1] = '\0';
        ts[i%150+1] = '\0';
        if(o) {
          fprintf(o, "%10u %s\n", qst, qs);
          fprintf(o, "%10s %s\n", " ", as);
          fprintf(o, "%10u %s\n", tst, ts);
          fprintf(o, "\n");
        }
        free(qs);
        free(as);
        free(ts);
        return 1;
      }
    } else if(kv_A(path, i) == MISMATCH) {
      qs[i%150] = rv ? compl[qseq[ql-1-q++]] : qseq[q++];
      as[i%150] = ' ';
      ts[i%150] = tseq[t++];
    } else if(kv_A(path, i) == INS) {
      qs[i%150] = rv ? compl[qseq[ql-1-q++]] : qseq[q++];
      as[i%150] = ' ';
      ts[i%150] = '-';
    } else if(kv_A(path, i) == DEL) {
      qs[i%150] = '-';
      as[i%150] = ' ';
      ts[i%150] = tseq[t++];
    }
    //fprintf(stderr, "i: %u\n", i);
    if(i%150 == 149 || i == kv_size(path)-1) {
      if(i == kv_size(path)-1) {
        qs[i%150+1] = '\0';
        as[i%150+1] = '\0';
        ts[i%150+1] = '\0';
      }
      if(o) {
        fprintf(o, "%10u %s\n", qst, qs);
        fprintf(o, "%10s %s\n", " ", as);
        fprintf(o, "%10u %s\n", tst, ts);
        fprintf(o, "\n");
      }
      qst = q;
      tst = t;
    }
  }
  free(qs);
  free(as);
  free(ts);
  return 0;
}

void output_paf(FILE* o, char* qn, int ql, int qs, int qe, chain ch, char* tn, int tl, int ts, int te, charvec fullpath, int score, int f1, int f2) {
  // minimap2 flags: NM:i:10828    ms:i:21848  AS:i:21066      nn:i:0 tp:A:P                                     cm:i:535   s1:i:5023   s2:i:0 de:f:0.1848 rl:i:0 cg:Z:100M1D100M...
  //                 edit distance             alignment score        type (P/primary, S/secondary, I/inversion) minimizers chain score                           CIGAR
  fprintf(o, "%s\t", qn);
  fprintf(o, "%d\t", ql);
  fprintf(o, "%d\t", qs);
  fprintf(o, "%d\t", qe);
  fprintf(o, "%c\t", "+-"[ch.rv]);
  fprintf(o, "%s\t", tn);
  fprintf(o, "%d\t", tl);
  fprintf(o, "%d\t", ts);
  fprintf(o, "%d\t", te);
  int matches = 0, total = 0;
  int i;
  for(i = 0; i < kv_size(fullpath); i++) {
    if(kv_A(fullpath, i) == MATCH) {
      matches++;
    }
    total++;
  }
  fprintf(o, "%d\t", matches);
  fprintf(o, "%d\t", total);
  // to match minimap2: mapQ = 40 * (1-f2/f1) * min(1,m/10) * ln(f1) where m is the number of anchors, f1 is the chaining score, f2 is the second best chain score
  int map_qual = 255;
  //map_qual = (int)(40.0 * (1-(f2 > 0 ? f2/f1 : 0)) * (1 < kv_size(ch.anchors)/10.0 ? 1 : kv_size(ch.anchors)/10.0) * log(f1));
  fprintf(o, "%d\t", map_qual);
  fprintf(o, "NM:i:%d\t", total-matches);
  fprintf(o, "AS:i:%d\t", score);
  fprintf(o, "tp:A:P\t");
  fprintf(o, "s1:i:%d\t", f1);
  fprintf(o, "s2:i:%d\t", f2);
  fprintf(o, "cg:Z:");
  char op;
  int ct = 0;
  for(i = 0; i < kv_size(fullpath); i++) {
    if(i == 0 || kv_A(fullpath, i) != op) {
      if(i > 0) {
        fprintf(o, "%d", ct);
        fprintf(o, "%c", "MIDX"[op]);
      }
      op = kv_A(fullpath, i);
      ct = 0;
    }
    ct++;
  }
  fprintf(o, "%d", ct);
  fprintf(o, "%c", "MIDX"[op]);
  fprintf(o, "\n");
}

void print_matrix(dp_cell **m, int q0, int q1, int t0, int t1, char* query, char* target) {
  int x, y;
  fprintf(stderr, "        ");
  for(x = t0+1; x <= t1; x++) {
    fprintf(stderr, "    %c", target[x-1]);
  }
  fprintf(stderr, "\n");
  for(y = q0; y <= q1; y++) {
    if(y > q0)
      fprintf(stderr, "%c |", query[y-1]);
    else
      fprintf(stderr, "  |");
    for(x = 0; x <= q1; x++) {
      fprintf(stderr, m[y][x].score < 0 ? " %c%3d" : "  %c%2d", "MIDX"[m[y][x].direction], m[y][x].score % 100);
    }
    fprintf(stderr, "\n");
  }
}

void print_bands(dp_cell **b, os* offset, int band_width, int q0, int q1, int t0, int t1, char* query, char* target) {
  int x, y;
  fprintf(stderr, "       ");
  for(x = t0; x < t1; x++) {
    fprintf(stderr, " %4d", x);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "       ");
  for(x = t0; x < t1; x++) {
    fprintf(stderr, "    %c", target[x]);
  }
  fprintf(stderr, "\n");
  for(y = q0; y < q1; y++) {
    fprintf(stderr, "%4d %c |", y, query[y]);
    for(x = t0; x < t1; x++) {
      dp_cell* c = x-offset[x+y].x >= 0 && x-offset[x+y].x < band_width ? &b[x+y][x-offset[x+y].x] : NULL;
      if(c)
        fprintf(stderr, c->score < 0 ? " %c%3d" : "  %c%2d", "MIDX"[c->direction], c->score % 100);
      else
        fprintf(stderr, "     ");
    }
    fprintf(stderr, "\n");
  }
}

// dynamic banded alignment, maybe w/SSE in the future
// query along y-axis, target along x
result align_banded(char* query, char* target, int qlen, int tlen, charvec *path, int rv, int verbose, int band_width) {
  if(verbose) {
    fprintf(stderr, "aligning q%c(%d) to t(%d)\n", "+-"[rv], qlen, tlen);
  }

  if(tlen == 0 || qlen == 0) {
    result res;
    res.failed = 1;
    // other fields are unset and unreliable
    return res;
  }

  dp_cell** bands = malloc((qlen+tlen-1) * sizeof(dp_cell*));
  int i, j;
  for(i = 0; i < qlen+tlen-1; i++) {
    bands[i] = malloc(band_width * sizeof(dp_cell));
  }

  os* offset = malloc((qlen+tlen-1) * sizeof(os));
  offset[0].x = 0;
  offset[0].y = 0;

  int x, y; // position of the LEFT END of the band
  int max_in_band;
  int match_score, ins_score, del_score;

  bands[0][0].score = 0;

  for(i = 0; i < qlen+tlen-1; i++) {
    x = offset[i].x;
    y = offset[i].y;
    //fprintf(stderr, "band starts at %d, %d\n", x, y);
    max_in_band = 0;

    for(j = 0; j < band_width; j++) {
      // find the appropriate cell for each 3 possible sources, if they exist
      dp_cell* prev = (i >= 2 && j + (offset[i].x - offset[i-2].x)/2 - (offset[i].y - offset[i-2].y)/2 < band_width && j + (offset[i].x - offset[i-2].x)/2 - (offset[i].y - offset[i-2].y)/2 >= 0 && y > 0 && x > 0 ? &bands[i-2][j + (offset[i].x - offset[i-2].x)/2 - (offset[i].y - offset[i-2].y)/2] : NULL);
      dp_cell* up = (i >= 1 && j + (offset[i].x - offset[i-1].x) < band_width && j + (offset[i].x - offset[i-1].x) >= 0 && y > 0 ? &bands[i-1][j + (offset[i].x - offset[i-1].x)] : NULL);
      dp_cell* left = (i >= 1 && j - (offset[i].y - offset[i-1].y) < band_width && j - (offset[i].y - offset[i-1].y) >= 0 && x > 0 ? &bands[i-1][j - (offset[i].y - offset[i-1].y)] : NULL);

      // do x,y
      // match
      if((rv ? compl[query[qlen-y-1]] : query[y]) == target[x]) {
        match_score = MATCH_SCORE;
      } else {
        match_score = MISMATCH_SCORE;
      }
      // construct the full ins and del gaps that would be needed to get to this spot if there is no prev
      match_score = (prev ? prev->score : score_gap(x) + score_gap(y)) + match_score; // gap score should be appropriate even for 0,0 and other edge cells

      // ins
      if(up && up->direction == INS) {
        ins_score = up->score - score_gap(up->gap_size) + score_gap(up->gap_size+1);
      } else {
        ins_score = up ? up->score + score_gap(1) : score_gap(x) + score_gap(y+1);
      }

      // del
      if(left && left->direction == DEL) {
        del_score = left->score - score_gap(left->gap_size) + score_gap(left->gap_size+1);
      } else {
        del_score = left ? left->score + score_gap(1) : score_gap(x+1) + score_gap(y);
      }

      // compare
      if(match_score >= ins_score && match_score >= del_score) {
        bands[i][j].score = match_score;
        if((rv ? compl[query[qlen-y-1]] : query[y]) == target[x]) {
          bands[i][j].direction = MATCH;
        } else {
          bands[i][j].direction = MISMATCH;
        }
        bands[i][j].gap_size = 0;
      } else if(ins_score >= del_score) {
        bands[i][j].score = ins_score;
        bands[i][j].direction = INS;
        bands[i][j].gap_size = up && up->direction == INS ? up->gap_size + 1 : 1;
      } else {
        bands[i][j].score = del_score;
        bands[i][j].direction = DEL;
        bands[i][j].gap_size = left && left->direction == DEL ? left->gap_size + 1 : 1;
      }

      // keep track of the best scoring position
      if(j > 0 && bands[i][j].score > bands[i][max_in_band].score) {
        max_in_band = j;
      }

      if(++x >= tlen)
        break;
      if(--y < 0)
        break;
    }

    //fprintf(stderr, "maximum score pos is %d, (x: %d, y: %d)\n", max_in_band, offset[i].x + max_in_band, offset[i].y - max_in_band);
    // recompute offsets, shifting 1 either right or down - right now we can't shift by more than 1 every round
    if(i < qlen+tlen-2) {
      if((max_in_band > band_width/2 && offset[i].y >= band_width) || offset[i].y == qlen-1) { // don't move horizontally until the full band is in bounds or we hit the bottom of the query
        offset[i+1].x = offset[i].x + 1;
        offset[i+1].y = offset[i].y;
      } else {
        offset[i+1].x = offset[i].x;
        offset[i+1].y = offset[i].y + 1;
      }
    }
  }

  if(verbose) {
    print_bands(bands, offset, band_width, 0, (30 < qlen ? 30 : qlen), 0, (30 < tlen ? 30 : tlen), query, target);
    if(qlen > 30)
      print_bands(bands, offset, band_width, 30, (60 < qlen ? 60 : qlen), 0, (30 < tlen ? 30 : tlen), query, target);
    if(tlen > 30)
      print_bands(bands, offset, band_width, 0, (30 < qlen ? 30 : qlen), 30, (60 < tlen ? 60 : tlen), query, target);
    if(qlen > 30 && tlen > 30)
      print_bands(bands, offset, band_width, 30, (60 < qlen ? 60 : qlen), 30, (60 < tlen ? 60 : tlen), query, target);
  }

  // backtrack - we must have ended at bands[qlen+tlen-1][0]
  i = qlen+tlen-2; // band index
  j = 0; // index into band of best pos
  x = offset[i].x;
  y = offset[i].y;
  while(i >= 0) {
    /*
    fprintf(stderr, "i = %d, j = %d, direction = %u, x = %d, y = %d, offset=%d,%d\n", i, j, bands[i][j].direction, x, y, offset[i].x, offset[i].y);
    char tmp = query[qlen];
    query[qlen] = '\0';
    fprintf(stderr, "query:  %s\n", query);
    query[qlen] = tmp;
    tmp = target[tlen];
    target[tlen] = '\0';
    fprintf(stderr, "target: %s\n", target);
    target[tlen] = tmp;
    fprintf(stderr, "q: %c, t: %c\n", (rv ? compl[query[qlen-y]] : query[y]), target[x]);
    */
    kv_push(unsigned char, *path, bands[i][j].direction);
    if(bands[i][j].direction == MATCH || bands[i][j].direction == MISMATCH) {
      // check that this is true...
      if(bands[i][j].direction == MATCH && (rv ? compl[query[qlen-y-1]] : query[y]) != target[x]) {
        fprintf(stderr, "these do NOT match!\n");
      }
      x--;
      y--;
      if(i >= 2)
        j = j + (offset[i].x - offset[i-2].x)/2 - (offset[i].y - offset[i-2].y)/2;
      i -= 2;
    } else if(bands[i][j].direction == INS) {
      y--;
      if(i >= 1)
        j = j + (offset[i].x - offset[i-1].x);
      i--;
    } else if(bands[i][j].direction == DEL) {
      x--;
      if(i >= 1)
        j = j - (offset[i].y - offset[i-1].y);
      i--;
    }
    if(j < 0 || j >= band_width || x < 0 || y < 0) { // fell off the side, everything else must be INDEL, then quit
      for(i = 0; i <= x; i++)
        kv_push(unsigned char, *path, DEL);
      for(i = 0; i <= y; i++)
        kv_push(unsigned char, *path, INS);
      break;
    }
  }

  result res;
  res.score = bands[qlen+tlen-2][0].score;
  res.qstart = 0;
  res.qend = offset[qlen+tlen-2].y;
  res.tstart = 0;
  res.tend = offset[qlen+tlen-2].x;
  // end positions are INCLUSIVE

  for(i = 0; i < qlen+tlen-1; i++) {
    free(bands[i]);
  }
  free(bands);
  free(offset);

  if(verbose) {
    fprintf(stderr, "  done aligning q%c(%d) to t(%d) with score %d\n", "+-"[rv], qlen, tlen, res.score);
  }

  return res;
} // align_banded

// query along y-axis, target along x
result align_full_matrix(char* query, char* target, int qlen, int tlen, charvec *path, int rv, int verbose) {
  if(verbose) {
    fprintf(stderr, "aligning q%c(%d) to t(%d)\n", "+-"[rv], qlen, tlen);
  }

  if(tlen == 0 || qlen == 0) {
    result res;
    res.failed = 1;
    // other fields are unset and unreliable
    return res;
  }

  if(verbose) {
    fprintf(stderr, "  allocating matrix with %d cells (%d bytes)\n", ((uint64_t)qlen+1)*(tlen+1), ((uint64_t)qlen+1)*(tlen+1)*7);
  }
  dp_cell** dp_matrix = malloc((qlen+1) * sizeof(dp_cell*));
  int i;
  for(i = 0; i <= qlen; i++) {
    dp_matrix[i] = malloc((tlen+1) * sizeof(dp_cell));
  }

  int x, y;
  int max_x = 0, max_y = 0;
  int match_score, ins_score, del_score;

  dp_matrix[0][0].score = 0;

  // MUST MATCH ENTIRETY OF BOTH query AND ref
  for(y = 1; y <= qlen; y++) {
    dp_matrix[y][0].score = score_gap(y);
    dp_matrix[y][0].direction = INS;
    dp_matrix[y][0].gap_size = y;
  }
  for(x = 1; x <= tlen; x++) {
    dp_matrix[0][x].score = score_gap(x);
    dp_matrix[0][x].direction = DEL;
    dp_matrix[0][x].gap_size = x;
  }

  for(y = 1; y <= qlen; y++) {
    for(x = 1; x <= tlen; x++) {
      // match
      if((rv ? compl[query[qlen-y]] : query[y-1]) == target[x-1]) {
        match_score = MATCH_SCORE;
      } else {
        match_score = MISMATCH_SCORE;
      }
      match_score = dp_matrix[y-1][x-1].score + match_score;

      // ins
      if(dp_matrix[y-1][x].direction == INS) {
        // check whether gap scores 1 or 2 are better for this size gap
        ins_score = dp_matrix[y-dp_matrix[y-1][x].gap_size][x].score + score_gap(dp_matrix[y-1][x].gap_size+1);
      } else {
        ins_score = dp_matrix[y-1][x].score + score_gap(1);
      }

      // del
      if(dp_matrix[y][x-1].direction == DEL) {
        // check whether gap scores 1 or 2 are better for this size gap
        del_score = dp_matrix[y][x-dp_matrix[y][x-1].gap_size].score + score_gap(dp_matrix[y][x-1].gap_size+1);
      } else {
        del_score = dp_matrix[y][x-1].score + score_gap(1);
      }

      // compare
      if(match_score >= ins_score && match_score >= del_score) {
        dp_matrix[y][x].score = match_score;
        if((rv ? compl[query[qlen-y]] : query[y-1]) == target[x-1]) {
          dp_matrix[y][x].direction = MATCH;
        } else {
          dp_matrix[y][x].direction = MISMATCH;
        }
        dp_matrix[y][x].gap_size = 0;
      } else if(ins_score >= del_score) {
        dp_matrix[y][x].score = ins_score;
        dp_matrix[y][x].direction = INS;
        dp_matrix[y][x].gap_size = dp_matrix[y-1][x].direction == INS ? dp_matrix[y-1][x].gap_size + 1 : 1;
      } else {
        dp_matrix[y][x].score = del_score;
        dp_matrix[y][x].direction = DEL;
        dp_matrix[y][x].gap_size = dp_matrix[y][x-1].direction == DEL ? dp_matrix[y][x-1].gap_size + 1 : 1;
      }
    }
  }

  /*
  int qst = qlen >= 71 ? 71-30 : (qlen >= 30 ? qlen - 30 : 0);
  int tmax = tlen > 30 ? 30 : tlen;
  int qmax = qst + 30 > qlen ? qlen : qst + 30;
  print_matrix(dp_matrix, qst, qmax, 0, tmax);
  */

  max_x = tlen;
  max_y = qlen;
  // compute maximum score position
  /*
  max_x = 0;
  max_y = qlen;
  for(x = 1; x < tlen; x++) { // check last row
    if(dp_matrix[qlen - 1][x].score > dp_matrix[max_y][max_x].score) {
      max_y = qlen - 1;
      max_x = x;
    }
  }
  for(y = 0; y < qlen; y++) { // check last column
    if(dp_matrix[y][tlen - 1].score > dp_matrix[max_y][max_x].score) {
      max_y = y;
      max_x = tlen - 1;
    }
  }
  */

  x = max_x;
  y = max_y;
  while(x > 0 || y > 0) {
    kv_push(unsigned char, *path, dp_matrix[y][x].direction);
    if(dp_matrix[y][x].direction == MATCH || dp_matrix[y][x].direction == MISMATCH) {
      x--;
      y--;
    } else if(dp_matrix[y][x].direction == INS) {
      y--;
    } else if(dp_matrix[y][x].direction == DEL) {
      x--;
    }
    if(verbose) {
      //fprintf(stderr, "  x: %d, y: %d\n", x, y);
    }
  }

  result res;
  res.score = dp_matrix[max_y][max_x].score;
  res.qstart = y;
  res.qend = max_y;
  res.tstart = x;
  res.tend = max_x;
  // end positions are INCLUSIVE

  // free these up in case we'll be doing this repeatedly
  for(i = 0; i <= qlen; i++) {
    free(dp_matrix[i]);
  }
  free(dp_matrix);

  if(verbose) {
    //fprintf(stderr, "  done aligning q%c(%d) to t(%d) with score %d\n", "+-"[rv], qlen, tlen, res.score);
  }

  return res;
}


int main(int argc, char *argv[]) {
  char* read_fasta = NULL;
  char* ref_fasta = NULL;
  char* method = "kmer";
  int verbose = 0;
  int limit = 0;
  int band_width = 500;

  // ---------- options ----------
  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "q:r:l:b:vh", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'q':
        read_fasta = optarg;
        break;
      case 'r':
        ref_fasta = optarg;
        break;
      case 'l':
        limit = atoi(optarg);
        break;
      case 'b':
        band_width = atoi(optarg);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'q' || optopt == 'r')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        // as long as all the long arguments have characters too, I don't think this section will be used
        if (long_idx == 0) verbose = 1; // --verbose
        else if (long_idx == 1) {usage(); return 0;} // --help
        else if (long_idx == 2) limit = atoi(optarg); // --limit
        else if (long_idx == 3) band_width = atoi(optarg); // --bandwidth
        break;
      default:
        usage();
        return 1;
    }
  }

  int index;
  char* command = NULL;
  for (index = optind; index < argc; index++) {
    if(index == optind) {
      command = argv[index];
    }
  }
  if(command == NULL) {
    usage();
    return 1;
  }

  if(read_fasta == NULL) {
    fprintf(stderr, "-q read FASTA/Q[.gz] is required\n");
    usage();
    return 1;
  }

  if(ref_fasta == NULL) {
    fprintf(stderr, "-r ref FASTA/Q[.gz] is required\n");
    usage();
    return 1;
  }

  int l;
  uint32_t i, j;
  int k = 21; // this should be a parameter

  // we'll reuse this vector for each alignment
  charvec path; // will just point to one of fw/rv
  charvec fw_path;
  kv_init(fw_path);
  charvec rv_path;
  kv_init(rv_path);

  // set up complement array
  for(i = 0; i < 256; i++) {
    compl[i] = 'N';
  }
  compl[65] = 'T';
  compl[67] = 'G';
  compl[71] = 'C';
  compl[84] = 'A';
  compl[97] = 't';
  compl[99] = 'g';
  compl[103] = 'c';
  compl[116] = 'a';

  // ---------- load ref FASTA file ----------
  gzFile f = gzopen(ref_fasta, "r");
  kseq_t* seq = kseq_init(f);
  fprintf(stderr, "Loading ref FASTA file: %s\n", ref_fasta);

  khash_t(faHash) *refmap = kh_init(faHash);
  kvec_t(fa_seq) refs;
  kv_init(refs);
  kvec_t(char*) refnames;
  kv_init(refnames);

  khint_t bin, bin2; // hash bin (result of kh_put)
  int absent;
  char* a;
  fa_seq fs;
  fa_seq fn;

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //printf("Reading %s (%i bp).\n", seq->name.s, l);

    // make <seq>
    fs.s = malloc(sizeof(char)*l);
    memcpy(fs.s, seq->seq.s, sizeof(char)*l);
    fs.l = l;

    // add <seq> to refs vector
    kv_push(fa_seq, refs, fs);

    // add <name>:<refs idx> to refmap
    a = malloc(strlen(seq->name.s)+1);
    strcpy(a, seq->name.s);
    bin = kh_put(faHash, refmap, a, &absent);
    kh_val(refmap, bin) = kv_size(refs)-1;
    kv_push(char*, refnames, a);
  }

  fprintf(stderr, "Loaded %d ref sequences.\n", kv_size(refs));

  kseq_destroy(seq);
  gzclose(f);

  // ---------- build k-mer hash (kmer : ref/pos) ------------
  khash_t(kmerHash) *h = kh_init(kmerHash);
  char *kmer;
  ref_pos rp;
  for(i = 0; i < kv_size(refs); i++) {
    rp.ref_id = i;
    for(j = 0; j < kv_A(refs, i).l - k + 1; j++) {
      kmer = malloc((k+1) * sizeof(char));
      kmer[k] = '\0';
      memcpy(kmer, kv_A(refs, i).s+j, k);
      //fprintf(stderr, "%s at %u:%u\n", kmer, i, j);
      bin = kh_put(kmerHash, h, kmer, &absent);
      if(absent) {
        kv_init(kh_val(h, bin));
      } else {
        free(kmer);
      }
      rp.pos = j;
      kv_push(ref_pos, kh_val(h, bin), rp);
    }
  }

  // ---------- load reads ----------
  f = gzopen(read_fasta, "r");
  seq = kseq_init(f);
  fprintf(stderr, "Reading fasta file: %s\n", read_fasta);

  char* rv_kmer;
  result aln, rv_aln;
  char strand;
  char tmp;
  int n = 0;

  posPair pp;
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    if(l < k) {
      if(verbose)
        fprintf(stderr, "%s is too short (%d bp)\n", seq->name.s, l);
      continue;
    }
    if(verbose) {
      fprintf(stderr, "Aligning %s (%i bp).\n", seq->name.s, l);
    }

    if(strcmp(method, "kmer") == 0) {
      khash_t(matchHash) *hits = kh_init(matchHash);

      for(j = 0; j < l - k + 1; j++) {
        if(j < l - k) { // except for the very last k-mer, make a fake substring
          tmp = seq->seq.s[j+k];
          seq->seq.s[j+k] = '\0';
        }

        //fprintf(stderr, "%s at %u:%u\n", kmer, i, j);
        bin = kh_get(kmerHash, h, seq->seq.s+j);
        if(j < l - k)
          seq->seq.s[j+k] = tmp;
        if(bin != kh_end(h)) { // hit something
          pp.qpos = j;
          for(i = 0; i < kv_size(kh_val(h, bin)); i++) {
            bin2 = kh_put(matchHash, hits, kv_A(kh_val(h, bin), i).ref_id, &absent);
            if(absent) {
              kv_init(kh_val(hits, bin2));
            }
            pp.tpos = kv_A(kh_val(h, bin), i).pos;
            kv_push(posPair, kh_val(hits, bin2), pp);
          }
        }

        // same for rev/compl kmer
        rv_kmer = rc(seq->seq.s+j, k, compl);
        bin = kh_get(kmerHash, h, rv_kmer);
        if(bin != kh_end(h)) { // hit something
          pp.qpos = l - k - j; // reverse the k-mer position so it represents the start position if the entire query were r/c
          for(i = 0; i < kv_size(kh_val(h, bin)); i++) {
            bin2 = kh_put(matchHash, hits, kv_A(kh_val(h, bin), i).ref_id + (1<<31), &absent);
            if(absent) {
              kv_init(kh_val(hits, bin2));
            }
            pp.tpos = kv_A(kh_val(h, bin), i).pos;
            kv_push(posPair, kh_val(hits, bin2), pp);
          }
        }
        free(rv_kmer);
      }
      if(verbose) {
        for(i = kh_begin(hits); i < kh_end(hits); i++) {
          if(kh_exist(hits, i)) {
            fprintf(stderr, "ref %u (%c) has %u hits\n", kh_key(hits, i) << 1 >> 1, "+-"[kh_key(hits, i)>>31], kv_size(kh_val(hits, i)));
          }
        }
      }

      chain *ch = do_chain(hits, 10000, 4, 10000, 3);
      for(i = kh_begin(hits); i < kh_end(hits); i++) {
        if(kh_exist(hits, i)) {
          kv_destroy(kh_value(hits, i));
        }
      }
      kh_destroy(matchHash, hits);

      // print all anchor sets
      if(verbose) {
        for(i = 0; ; i++) {
          if(kv_size(ch[i].anchors) > 0) {
            fprintf(stderr, "ref %u (%c), score %d, anchors:", ch[i].ref, "+-"[ch[i].rv], ch[i].score);
            for(j = 0; j < kv_size(ch[i].anchors); j++) {
              fprintf(stderr, " %u:%u", kv_A(ch[i].anchors, j).qpos, kv_A(ch[i].anchors, j).tpos);
            }
            fprintf(stderr, "\n");
          } else {
            break;
          }
        }
      }

      if(kv_size(ch[0].anchors) == 0) {
        if(verbose)
          fprintf(stderr, "No chain (no k-mer hits) for read '%s'\n", seq->name.s);
        free_chains(ch);
        continue;
      }

      uint32_t q = 0, t = 0;

      int score = 0;
      charvec fullpath;
      kv_init(fullpath);
      uint32_t qe, te;
      for(j = 0; j <= kv_size(ch[0].anchors); j++) {
        qe = j < kv_size(ch[0].anchors) ? kv_A(ch[0].anchors, j).qpos : l;
        te = j < kv_size(ch[0].anchors) ? kv_A(ch[0].anchors, j).tpos : kv_A(refs, ch[0].ref).l;
        if(q > 0 && qe <= q && te <= t) {
          // both overlap, but they are NOT necessarily identical, but they must overlap a homopolymer stretch
          // add a match for every base they *don't* overlap AND are the same length, up to the end of the read
          uint32_t min_len = qe + k - q < te + k - t ? qe + k - q : te + k - t;
          uint32_t max_len = qe + k - q > te + k - t ? qe + k - q : te + k - t;
          if(verbose && (qe+k-q > 1 || te+k-t > 1)) {
            fprintf(stderr, "inter-kmer match of length  of %u (q) and %u (t)\n", qe + k - q, te + k - t);
          }
          for(i = 0; i < min_len; i++)
            kv_push(unsigned char, fullpath, MATCH);
          for(i = 0; i < max_len - min_len; i++) {
            if(qe + k - q > te + k - t)
              kv_push(unsigned char, fullpath, INS);
            if(qe + k - q < te + k - t)
              kv_push(unsigned char, fullpath, DEL);
          }
        } else {
          if(qe <= q && te > t) {
            if(j > 0)
              for(i = 0; i < (k+qe-q); i++)
                kv_push(unsigned char, fullpath, MATCH);
            if(verbose) {
              fprintf(stderr, "DEL of length %u\n", ((te-t) - (qe-q)));
            }
            for(i = 0; i < (te-t) - (qe-q); i++)
              kv_push(unsigned char, fullpath, DEL);
          } else if(te <= t && qe > q) {
            if(j > 0)
              for(i = 0; i < (k+te-t); i++)
                kv_push(unsigned char, fullpath, MATCH);
            if(verbose) {
              fprintf(stderr, "INS of length %u\n", ((qe-q) - (te-t)));
            }
            for(i = 0; i < (qe-q) - (te-t); i++)
              kv_push(unsigned char, fullpath, INS);
          } else {
            if(j > 0)
              for(i = 0; i < k; i++)
                kv_push(unsigned char, fullpath, MATCH);
            fw_path.n = 0; // reset vector to reuse
            rv_path.n = 0;
            if(verbose) {
              fprintf(stderr, "aligning q %u-%u to t %u-%u\n", q, qe, t, te);
              char tmp = *(seq->seq.s + (ch[0].rv ? l-q : qe));
              *(seq->seq.s + (ch[0].rv ? l-q : qe)) = '\0';
              fprintf(stderr, "q: %s\n", seq->seq.s + (ch[0].rv ? l-qe : q));
              *(seq->seq.s + (ch[0].rv ? l-q : qe)) = tmp;
              tmp = *(kv_A(refs, ch[0].ref).s+te);
              *(kv_A(refs, ch[0].ref).s+te) = '\0';
              fprintf(stderr, "t: %s\n", kv_A(refs, ch[0].ref).s+t);
              *(kv_A(refs, ch[0].ref).s+te) = tmp;
            }
            //aln = align_full_matrix(seq->seq.s + (ch[0].rv ? l-qe : q), kv_A(refs, ch[0].ref).s+t, qe - q, te - t, &fw_path, ch[0].rv, verbose);
            aln = align_banded(seq->seq.s + (ch[0].rv ? l-qe : q), kv_A(refs, ch[0].ref).s+t, qe - q, te - t, &fw_path, ch[0].rv, verbose, band_width);
            score += aln.score;
            for(i = 0; i < kv_size(fw_path); i++) { // path is reversed from alignment - we have to do it this way instead of i-- because i is unsigned!!
              kv_push(unsigned char, fullpath, kv_A(fw_path, kv_size(fw_path)-1-i));
              if(verbose)
                fprintf(stderr, "%u", kv_A(fw_path, kv_size(fw_path)-1-i));
            }
            if(verbose)
              fprintf(stderr, "\n");
          }
        }
        if(j < kv_size(ch[0].anchors)) {
          q = kv_A(ch[0].anchors, j).qpos + k;
          t = kv_A(ch[0].anchors, j).tpos + k;
          //fprintf(stderr, "q %u :: t %u\n", q, t);
        }
      }

      int f1 = ch[0].score;
      int f2 = ch[1].score;
      output_paf(stdout, seq->name.s, l, 0, l, ch[0], kv_A(refnames, ch[0].ref), kv_A(refs, ch[0].ref).l, 0, kv_A(refs, ch[0].ref).l, fullpath, score, f1, f2);
      int res = output_human_aln(seq->seq.s, l, kv_A(refs, ch[0].ref).s, ch[0].rv, fullpath, (verbose ? stdout : NULL));

      kv_destroy(fullpath);
      free_chains(ch);

      if(res != 0) {
        break;
      }
      if(limit > 0 && ++n >= limit) {
        break;
      }
    }

    else if(strcmp(method, "full_dp") == 0) {
      fw_path.n = 0; // reset vector to reuse
      rv_path.n = 0;

      aln = align_full_matrix(seq->seq.s, kv_A(refs, 0).s, l, kv_A(refs, 0).l, &fw_path, 0, verbose);

      //char* rv = rc(seq->seq.s, l, compl);
      //rv_aln = align_full_matrix(rv, kv_A(refs, 0).s, l, kv_A(refs, 0).l, &rv_path, 1, verbose);
      //free(rv);
      rv_aln = align_full_matrix(seq->seq.s, kv_A(refs, 0).s, l, kv_A(refs, 0).l, &rv_path, 1, verbose);

      strand = '+';
      path = fw_path;
      if(rv_aln.score > aln.score) {
        aln = rv_aln;
        path = rv_path;
        strand = '-';
      }

      printf("Score: %d\n", aln.score);
      printf("Strand: %c\n", strand);
      printf("Query: %d - %d\n", aln.qstart, aln.qend);
      printf("Target: %d - %d\n", aln.tstart, aln.tend);
      printf("Path length: %d\n", kv_size(path));
      float acc = cigar_accuracy(path.a, kv_size(path));
      printf("Accuracy: %f\n", acc);
    }

    else {
      fprintf(stderr, "ERROR: unknown method '%s'\n", method);
      usage();
      return 1;
    }
  }

  kv_destroy(fw_path);
  kv_destroy(rv_path);
  for(i = 0; i < kv_size(refs); i++) {
    free(kv_A(refs, i).s);
  }
  kv_destroy(refs);
  kv_destroy(refnames); // name strings themselves are reused in refmap and are freed below
  for(i = kh_begin(refmap); i < kh_end(refmap); i++) {
    if(kh_exist(refmap, i)) {
      free((void*)kh_key(refmap, i));
    }
  }
  kh_destroy(faHash, refmap);
  for(i = kh_begin(h); i < kh_end(h); i++) {
    if(kh_exist(h, i)) {
      kv_destroy(kh_value(h, i));
      free((void*)kh_key(h, i));
    }
  }
  kh_destroy(kmerHash, h);

  kseq_destroy(seq);
  gzclose(f);
}

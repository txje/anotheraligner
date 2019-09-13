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

#ifndef __CHAIN_H__
#define __CHAIN_H__

typedef struct {
  uint32_t qpos;
  uint32_t tpos;
} posPair;

typedef kvec_t(posPair) pairVec;

// creates uint32(ref id):kvec<qpos,tpos> hash
KHASH_MAP_INIT_INT(matchHash, pairVec);

typedef struct score_pos {
  int score;
  int anchor_idx;
  int score_idx;
  uint32_t ref;
  int prev;
  uint8_t used;
  uint8_t rv;
} score_pos;

typedef kvec_t(score_pos) scoreVec;

typedef struct chain {
  int score;
  pairVec anchors;
  uint32_t ref;
  uint8_t rv;
} chain;

chain* do_chain(khash_t(matchHash) *hits, khash_t(matchHash) *rv_hits, int max_chains, int match_score, int max_gap, int min_chain_length);

void free_chains(chain* chs);

#endif /* __CHAIN_H__ */

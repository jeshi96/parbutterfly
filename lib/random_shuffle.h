// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once
#include <math.h>
#include <stdio.h>
#include <cstdint>
#include "utilities.h"
#include "random.h"
#include "counting_sort.h"

namespace pbbs {

  template <typename Seq>
  void seq_random_shuffle(Seq A, random r = default_random) {
    size_t n = A.size();
    // the Knuth shuffle
    if (n < 2) return;
    for (size_t i=n-1; i > 0; i--)
      std::swap(A[i],A[r.ith_rand(i)%(i+1)]);
  }
  
  template <typename Seq>
  void random_shuffle(Seq A, random r = default_random) {
    size_t n = A.size();
    if (n < SEQ_THRESHOLD) {
      seq_random_shuffle(A);
      return;
    }

    size_t bits = 10;
    if (n < (1 << 27))
      bits = (log2_up(n) - 7)/2;
    else bits = (log2_up(n) - 17);
    size_t num_buckets = (1<<bits);
    size_t mask = num_buckets - 1;
    auto rand_pos = [&] (size_t i) -> size_t {
      return r.ith_rand(i) & mask;};

    auto get_pos = make_sequence<size_t>(n, rand_pos);

    // first randomly sorts based on random values [0,num_buckets)
    sequence<size_t> bucket_offsets = count_sort(A, A, get_pos, num_buckets);
	
    // now sequentially randomly shuffle within each bucket
    parallel_for_1(size_t i = 0; i < num_buckets; i++) {
      size_t start = bucket_offsets[i];
      size_t end = bucket_offsets[i+1];
      seq_random_shuffle(A.slice(start,end), r.fork(i));
    }
  }
}

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
#include <stdio.h>
#include <math.h>
#include "utilities.h"
#include "sequence_ops.h"
#include "transpose.h"

namespace pbbsa {

  // the following parameters can be tuned
  constexpr const size_t SEQ_THRESHOLD = 4196;
  constexpr const size_t BUCKET_FACTOR = 32;
  constexpr const size_t LOW_BUCKET_FACTOR = 16;

  // Sequential internal version
  template <typename b_size_t, typename s_size_t, typename InSeq, typename KeySeq>
  void _seq_count_sort(InSeq In, typename InSeq::T* Out, KeySeq Keys,
		       s_size_t* counts, size_t num_buckets) {
    size_t n = In.size();
    sequence<b_size_t> tmp(n);

    for (size_t i = 0; i < num_buckets; i++) {
      counts[i] = 0;
    }
    for (size_t j = 0; j < n; j++) {
      size_t k = tmp[j] = Keys[j];
      if (k >= num_buckets) abort();
      counts[k]++;
    }
    // generate offsets
    size_t s = 0;
    for (size_t i = 0; i < num_buckets; i++) {
      size_t t = counts[i];
      counts[i] = s;
      s += t;
    }
    for (size_t j = 0; j < n; j++) {
      size_t k = counts[tmp[j]]++;
      // needed for types with self defined assignment or initialization
      // otherwise equivalent to: Out[k+start] = In[j+start];
      move_uninitialized(Out[k], In[j]);
    }
    // convert offsets back to counts
    s = 0;
    for (size_t i = 0; i < num_buckets; i++) {
      size_t t = counts[i];
      counts[i] = t - s;
      s = t;
    }
  }

  // Sequential version
  template <typename b_size_t, typename InS, typename OutS, typename KeyS>
  sequence<size_t> seq_count_sort(InS In, OutS Out, KeyS Keys, size_t num_buckets) {
    using T = typename InS::T;
    size_t n = In.size();

    size_t* counts = new_array_no_init<size_t>(num_buckets+1);
    T* B = new_array_no_init<T>(n);
    _seq_count_sort<b_size_t,size_t>(In, B, Keys, counts, num_buckets);
    for (size_t i = 0; i < n ; i++)
      Out.update(i,B[i]);
    free(B);
    size_t c = 0;
    for (size_t i=0; i < num_buckets; i++) {
      size_t x = counts[i];
      counts[i] = c;
      c = c + x;
    }
    counts[num_buckets] = n;
    return sequence<size_t>(counts,num_buckets+1);
  }

  // Parallel internal version
  template <typename b_size_t, typename s_size_t,
    typename InS, typename OutS, typename KeyS>
  sequence<size_t> _count_sort(InS In, OutS Out, KeyS Keys, size_t num_buckets) {
    using T = typename InS::T;
    size_t n = In.size();

    // pad to 16 buckets to avoid false sharing (does not affect results)
    num_buckets = std::max(num_buckets, (size_t) 16);

    //s_size_t num_blocks = (n/(num_buckets*BUCKET_FACTOR) + 1);
    size_t sqrt = (size_t) ceil(pow(n,0.5));
    size_t num_blocks =
      (size_t) (n < (1 << 24)) ? (sqrt/32) : ((n < (1 << 28)) ? sqrt/2 : sqrt);
    s_size_t num_threads = __cilkrts_get_nworkers();
    num_blocks = 1 << log2_up(num_blocks);

    // if insufficient parallelism, sort sequentially
    if (n < SEQ_THRESHOLD || num_blocks == 1 || num_threads == 1)
      return seq_count_sort<b_size_t>(In,Out,Keys,num_buckets);

    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;

    T *B = new_array_no_init<T>(n);
    s_size_t *counts = new_array_no_init<s_size_t>(m,1);

    // sort each block
    parallel_for_1 (size_t i = 0; i < num_blocks; ++i) {
      s_size_t start = std::min(i * block_size, n);
      s_size_t end =  std::min(start + block_size, n);
      _seq_count_sort<b_size_t,s_size_t>(In.slice(start,end), B+start,
					 Keys.slice(start,end),
					 counts + i*num_buckets, num_buckets);
    }

    OutS C = Out;
    size_t* bucket_offsets = transpose_buckets(B, Out.as_array(),
					       counts, n, block_size,
					       num_blocks, num_buckets);
    free(counts); free(B);
    return sequence<size_t>(bucket_offsets,num_buckets+1);
  }

  // Sequential internal version
  template <typename b_size_t, typename s_size_t, typename InSeq, typename KeySeq>
  void _seq_count_sort_no_transpose(InSeq In, typename InSeq::T* Out, KeySeq Keys,
      size_t start, size_t end, s_size_t* counts, size_t num_buckets) {

    size_t n = end-start;
    s_size_t* offsets = new_array_no_init<s_size_t>(num_buckets);
    b_size_t* tmp = new_array_no_init<b_size_t>(n);

    for (size_t i = 0; i < num_buckets; i++) {
      offsets[i] = 0;
    }
    for (size_t j = 0; j < n; j++) {
      size_t k = tmp[j] = Keys[j + start];
      offsets[k]++;
    }
    // generate offsets
    size_t s = 0;
    for (size_t i = 0; i < num_buckets; i++) {
      counts[i] = s;
      s += offsets[i];
      offsets[i] = s;
    }
    for (long j = ((long) n)-1; j >= 0; j--) {
      s_size_t k =  --offsets[tmp[j]];
      // needed for types with self defined assignment or initialization
      // otherwise equivalent to: Out[k+start] = In[j+start];
      move_uninitialized(Out[k+start], In[j+ start]);
    }
    free(offsets); free(tmp);
  }

  // Returns the untransposed matrix of size num_blocks*num_buckets, and the
  // number of blocks used.
  template <typename b_size_t, typename s_size_t,
    typename InS, typename KeyS>
  tuple<typename InS::T*, s_size_t*, size_t> _count_sort_no_transpose(InS In, KeyS Keys, size_t num_buckets) {
    using T = typename InS::T;
    size_t n = In.size();

    // pad to 16 buckets to avoid false sharing (does not affect results)
    num_buckets = std::max(num_buckets, (size_t) 16);

    size_t sqrt = (size_t) ceil(pow(n,0.5));
    size_t num_blocks = sqrt/32; // magic.
    s_size_t num_threads = __cilkrts_get_nworkers();
    num_blocks = 1 << log2_up(num_blocks);

    // if insufficient parallelism, sort sequentially
    if (n < SEQ_THRESHOLD || num_blocks == 1 || num_threads == 1) {
      T* Out = new_array_no_init<T>(n);
      s_size_t* counts = new_array_no_init<s_size_t>(num_buckets+1);
      _seq_count_sort_no_transpose<b_size_t,s_size_t>(In,Out,Keys,0, n, counts,num_buckets);
      return make_tuple(Out, counts, 1);
    }

    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;

//    cout << "num_blocks = " << num_blocks << " num_buckets = " << num_buckets << endl;
    T *B = new_array_no_init<T>(n);
    s_size_t *counts = new_array_no_init<s_size_t>(m,1);

    // sort each block
    parallel_for_1 (size_t i = 0; i < num_blocks; ++i) {
      size_t start = std::min(i * block_size, (size_t)n);
      size_t end =  std::min(start + block_size, (size_t)n);
      _seq_count_sort_no_transpose<b_size_t,s_size_t>(In, B,
					 Keys, start, end,
					 counts + i*num_buckets, num_buckets);
    }

    return make_tuple(B, counts, num_blocks);
  }

  template <typename s_size_t, typename InS, typename KeyS>
  tuple<typename InS::T*, s_size_t*, size_t>
  _count_sort_no_transpose_size(InS In, KeyS Keys, size_t num_buckets) {
    using T = typename InS::T;
    if (num_buckets <= 256) {
      return _count_sort_no_transpose<uint8_t, T>(In, Keys, num_buckets);
    } else if (num_buckets <= (1 << 16)) {
      return _count_sort_no_transpose<uint16_t, T>(In, Keys, num_buckets);
    } else {
      return _count_sort_no_transpose<T, T>(In, Keys, num_buckets);
    }
  }

//  template <typename InS, typename KeyS>
//  tuple<typename InS::T*, sequence<size_t>, size_t>
//  count_sort_no_transpose(InS In, KeyS Keys, size_t num_buckets) {
//    using T = typename InS::T;
//    size_t n = In.size();
//    size_t max32 = ((size_t) 1) << 32;
//    if (n < max32 && num_buckets < max32)
//      // use 4-byte counters when larger ones not needed
//      return _count_sort_no_transpose<uint32_t>(In, Out, Keys, num_buckets);
//    return _count_sort_no_transpose<size_t>(In, Out, Keys, num_buckets);
//  }

  template <typename s_size_t, typename InS, typename OutS, typename KeyS>
  sequence<size_t> _count_sort_size(InS In, OutS Out, KeyS Keys, size_t num_buckets) {
    if (num_buckets <= 256)
      return _count_sort<uint8_t,s_size_t>(In, Out, Keys, num_buckets);
    else if  (num_buckets <= (1 << 16))
      return _count_sort<uint16_t,s_size_t>(In, Out, Keys, num_buckets);
    else
      return _count_sort<s_size_t,s_size_t>(In, Out, Keys, num_buckets);
  }

  // Parallel version
  template <typename InS, typename OutS, typename KeyS>
  sequence<size_t> count_sort(InS In, OutS Out, KeyS Keys, size_t num_buckets) {
    size_t n = In.size();
    size_t max32 = ((size_t) 1) << 32;
    if (n < max32 && num_buckets < max32)
      // use 4-byte counters when larger ones not needed
      return _count_sort_size<uint32_t>(In, Out, Keys, num_buckets);
    return _count_sort_size<size_t>(In, Out, Keys, num_buckets);
  }
}

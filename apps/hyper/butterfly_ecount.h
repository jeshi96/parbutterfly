#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "math.h"
#include "graph.h"
#include "vertex.h"
#include "sequence.h"
#include "sparseSet.h"
#include "sampleSort.h"
#ifndef OPENMP
#include "../../lib/histogram.h"
#include "../../lib/sample_sort.h"
#endif
#include "../../radixsort/RadixSort/radixSort.h"
#include "../../radixsort/common/blockRadixSort.h"

#include "butterfly_utils.h"

using namespace std;

/*
 *  Computes butterfly counts per edge, given a ranked graph and using
 *  histogram aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHistCE(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
                  intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<long,uintE>;
  using T = tuple<uintE, long>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair  = getWedges<long>(cs.wedges_seq_int, GA, UWedgeIntRankCons(GA.n), max_wedges, curr_idx,
                                                  num_wedges, wedge_idxs);
  intT next_idx = wedges_pair.second;
  long num_wedges_list = wedges_pair.first;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(cs.wedges_seq_int.A, wedges_pair.first);

  // Construct a histogram to aggregate wedges
  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  // Set up a hash table to hold wedge counts
  cs.wedges_hash.resize(wedge_freqs_n);
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

  // Set up an array to hold butterfly counts
  if (cs.butterflies_seq_intt.n < 2*num_wedges_list) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(T, 2*num_wedges_list);
    cs.butterflies_seq_intt.n = 2*num_wedges_list;
  }

  // Retrieve aggregated wedge counts and store butterfly count per edge in another array
  // Iterate through all wedges with endpoint i
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    long idx = 0;
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    for(intT j=0;j<u_deg;j++) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
#ifdef INVERSE
      {
#else
      if (v > i) {
#endif
	// Iterate through all wedges with endpoints i and u2, and center v
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2 = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
    if (u2 < i && u2 < v) {
#else
	  if (u2 > i) {
#endif
	    // Find the number of wedges with endpoints i and u2
	    long to_find = ((((long) i) *GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	    long num_butterflies = cs.wedges_hash.find(to_find).second;
	    if (num_butterflies > 1) {
	      // Store butterfly counts per edge, (i, v) and (u2, v)
	      cs.butterflies_seq_intt.A[2*(wedge_idx+idx)] = make_tuple((v_offset+k), num_butterflies-1);
	      cs.butterflies_seq_intt.A[2*(wedge_idx+idx)+1] = make_tuple((u_offset+j), num_butterflies-1);
	    }
	    else {
	      cs.butterflies_seq_intt.A[2*(wedge_idx+idx)] = make_tuple(UINT_E_MAX, 0);
	      cs.butterflies_seq_intt.A[2*(wedge_idx+idx)+1] = make_tuple(UINT_E_MAX, 0);
	    }
	    ++idx;
	  }
	  else break;
	}
      }
    }
  }

  // Aggregate butterfly counts by constructing another histogram
  pbbsa::sequence<T> wedge_freqs_i_seq = pbbsa::sequence<T>(cs.butterflies_seq_intt.A,2*num_wedges_list);
  tuple<size_t, T*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(wedge_freqs_i_seq,GA.numEdges, getAdd<uintE,long>, getAddReduce<uintE,long>,
                                          cs.tmp_uint, cs.out_uint);
  T* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  // Store butterfly counts per edge
  parallel_for (long i=0; i < butterflies_n; ++i) {
    if (get<0>(butterflies_l[i]) != UINT_E_MAX)
      butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  return next_idx;
#endif
}

/*
 *  Computes butterfly counts per edge, using by side ranking and histogram
 *  aggregation (for wedge and butterfly aggregation).
 * 
 *  cs           : Holds all array space needed, to be reused between batches
 *  GA           : Bipartite graph in CSR format
 *  use_v        : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                 produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges   : Number of wedges produced by the specified bipartition
 *  butterflies  : Array to store butterfly counts per edge, considering edges derived from one orientation (as given
 *                 by use_v)
 *  butterflies_u: Array to store butterfly counts per edge, considering edges derived from the other orientation (as 
 *                 given by use_v)
 *  max_wedges   : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                 memory
 *  wedges_idxs  : Wedge indices to allow for wedge retrieval in parallel
 *  eti          : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  curr_idx     : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHistCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
                  long* wedge_idxs, uintE* eti, intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  using X = tuple<long,uintE>;
  using T = tuple<uintE, long>;
  auto cons = UVertexPairIntCons(nu);

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, cons, max_wedges, curr_idx, num_wedges,
                                                 wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_list = wedges_pair.first;
  intT next_idx = wedges_pair.second;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_list);

  // Construct a histogram to aggregate wedges
  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  // Set up a hash table to hold wedge counts
  cs.wedges_hash.resize(wedge_freqs_n);
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

  // Set up an array to hold butterfly counts
  if (cs.butterflies_seq_intt.n < 2*num_wedges_list) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(T, 2*num_wedges_list);
    cs.butterflies_seq_intt.n = 2*num_wedges_list;
  }

  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve aggregated wedge counts and store butterfly count per edge in another array
  // Iterate through all wedges with endpoint i
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    long idx = 0;
    for(intT j=0;j<u_deg;j++) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Iterate through all wedges with endpoints i and u2, and center v
      for (intT k=0; k < v_deg; ++k) {
	uintE u2 = edgesV[v_offset+k];
	if (u2 < i) {
	  // Find the number of wedges with endpoints i and u2
	  long num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
	  // Store butterfly counts per edge, (i, v) and (u2, v)
	  cs.butterflies_seq_intt.A[2*(wedge_idx+idx)] = make_tuple(eti[u_offset + j], num_butterflies);
	  cs.butterflies_seq_intt.A[2*(wedge_idx+idx)+1] = make_tuple((v_offset + k), num_butterflies);
	  ++idx;
	}
	else break;
      }
    }
  }

  // Aggregate butterfly counts by constructing another histogram
  pbbsa::sequence<T> wedge_freqs_i_seq = pbbsa::sequence<T>(cs.butterflies_seq_intt.A,2*num_wedges_list);
  tuple<size_t, T*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(wedge_freqs_i_seq,GA.numEdges,getAdd<uintE,long>, getAddReduce<uintE,long>,
                                          cs.tmp_uint, cs.out_uint);
  T* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  // Store butterfly counts per edge
  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  return next_idx;
#endif
}

/*
 *  Computes butterfly counts per edge, using by side ranking and histogram
 *  aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs           : Holds all array space needed, to be reused between batches
 *  GA           : Bipartite graph in CSR format
 *  use_v        : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                 produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges   : Number of wedges produced by the specified bipartition
 *  butterflies  : Array to store butterfly counts per edge, considering edges derived from one orientation (as given
 *                 by use_v)
 *  butterflies_u: Array to store butterfly counts per edge, considering edges derived from the other orientation (as 
 *                 given by use_v)
 *  max_wedges   : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                 memory
 *  wedges_idxs  : Wedge indices to allow for wedge retrieval in parallel
 *  eti          : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  curr_idx     : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHist(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long* butterflies_u,
                long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  using X = tuple<long,uintE>;
  auto cons = UVertexPairIntCons(nu);

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, cons, max_wedges, curr_idx, num_wedges,
                                                 wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_list = wedges_pair.first;
  intT next_idx = wedges_pair.second;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_list);

  // Construct a histogram to aggregate wedges
  tuple<size_t, X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  // Set up a hash table to hold wedge counts
  cs.wedges_hash.resize(wedge_freqs_n);
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve aggregated wedge counts and store butterfly count per edge
  // Iterate through all wedges with endpoint i
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Iterate through all wedges with endpoints i and u2, and center v
      for (intT k=0; k < v_deg; ++k) {
	uintE u2 = edgesV[v_offset+k];
	if (u2 < i) {
	  // Find the number of wedges with endpoints i and u2
	  long num_butterflies = cs.wedges_hash.find(i * nu + u2).second - 1;
	  // Store butterfly counts per edge, (i, v) and (u2, v)
	  writeAdd(&butterflies_u[eltsPerCacheLine*(u_offset + j)],num_butterflies); 
	  writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)],num_butterflies); 
	}
	else break;
      }
    }
  }

  return next_idx;
#endif
}

/*
 *  Computes butterfly counts per edge, given a ranked graph and using
 *  histogram aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHist(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
                intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<long,uintE>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair  = getWedges<long>(cs.wedges_seq_int, GA, UWedgeIntRankCons(GA.n), max_wedges, curr_idx,
                                                  num_wedges, wedge_idxs);
  intT next_idx = wedges_pair.second;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(cs.wedges_seq_int.A, wedges_pair.first);

  // Construct a histogram to aggregate wedges
  tuple<size_t, X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp,
                                                                        cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  // Set up a hash table to hold wedge counts
  cs.wedges_hash.resize(wedge_freqs_n);
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

  // Retrieve aggregated wedge counts and store butterfly count per edge
  // Iterate through all wedges with endpoint i
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
#ifdef INVERSE
      {
#else
      if (v > i) {
#endif
	// Iterate through all wedges with endpoints i and u2, and center v
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2 = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
    if (u2 < i && u2 < v) {
#else
	  if (u2 > i) {
#endif
	    // Find the number of wedges with endpoints i and u2
	    long to_find = ((((long) i) *GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	    long num_butterflies = cs.wedges_hash.find(to_find).second;
	    if (num_butterflies > 1) {
	      // Store butterfly counts per edge, (i, v) and (u2, v)
	      writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], num_butterflies - 1);
	      writeAdd(&butterflies[eltsPerCacheLine*(u_offset+j)], num_butterflies - 1);
	    }
	  }
	  else break;
	}
      }
    }
  }

  return next_idx;
#endif
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per edge, using by side ranking and sort
 *  aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *               produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  eti        : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountESortCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
                  long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, use_v, UWedgeCons(), max_wedges, curr_idx,
                                                    num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;

  using X = tuple<uintE,long>;
  if (cs.butterflies_seq_intt.n < 2*num_wedges_f) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*num_wedges_f);
    cs.butterflies_seq_intt.n = 2*num_wedges_f;
  }

  // Retrieve frequency counts for all wedges with the same key
  // First, sort wedges to aggregate
#ifdef OPENMP
  sampleSort(wedges, num_wedges_f, UWedgeCmp());
#else
  pbbs::sample_sort(wedges, num_wedges_f, UWedgeCmp());
#endif
  // Retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  // Store butterfly counts per edge in another array
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num = freq_arr[i+1] - freq_arr[i] - 1;
    // Iterate through all wedges
    granular_for (j, freq_arr[i], freq_arr[i+1], (freq_arr[i+1]-freq_arr[i] > 1000), {
	// Store butterfly counts per edge (v_offset + k and u_offset + j respectively)
	cs.butterflies_seq_intt.A[(long)2*j] = make_tuple(eti[offsetsU[wedges[j].v1] + wedges[j].j], num);
	cs.butterflies_seq_intt.A[(long)2*j+1] = make_tuple(offsetsV[wedges[j].u] + wedges[j].k, num);
      });
  }

  free(freq_arr);

  // Aggregate butterfly counts
  // Sort butterfly counts to aggregate
#ifdef OPENMP
  sampleSort(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>());
#else
  pbbs::sample_sort(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>());
#endif
  // Retrieve a list of indices where consecutive wedges have different keys
  auto b_freq_pair = getFreqs(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>(), tupleEq<uintE,long>(),
                              LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;

  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Store butterfly counts per edge
  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    long num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]), num_freq, tupleAdd<uintE,long>());
    // These are our butterfly counts
    butterflies[eltsPerCacheLine*get<0>(reduce)] += get<1>(reduce);
  }

  free(b_freq_arr);

  return wedges_pair.second;
}

/*
 *  Computes butterfly counts per edge, using by side ranking and sort
 *  aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs           : Holds all array space needed, to be reused between batches
 *  GA           : Bipartite graph in CSR format
 *  use_v        : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                 produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges   : Number of wedges produced by the specified bipartition
 *  butterflies  : Array to store butterfly counts per edge, considering edges derived from one orientation (as given
 *                 by use_v)
 *  butterflies_u: Array to store butterfly counts per edge, considering edges derived from the other orientation (as 
 *                 given by use_v)
 *  max_wedges   : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                 memory
 *  wedges_idxs  : Wedge indices to allow for wedge retrieval in parallel
 *  eti          : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  curr_idx     : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountESort(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long* butterflies_u,
                long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, use_v, UWedgeCons(), max_wedges, curr_idx,
                                                    num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, sort wedges to aggregate
#ifdef OPENMP
  sampleSort(wedges, num_wedges_f, UWedgeCmp());
#else
  pbbs::sample_sort(wedges, num_wedges_f, UWedgeCmp());
#endif
  // Retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Store butterfly counts per edge
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num = freq_arr[i+1] - freq_arr[i] - 1;
    if (num > 0) {
      // Iterate through all wedges
      granular_for (j, freq_arr[i], freq_arr[i+1], (freq_arr[i+1]-freq_arr[i] > 1000), {
	  // Store butterfly counts per edge (v_offset + k and u_offset + j respectively)
	  writeAdd(&butterflies_u[eltsPerCacheLine*(offsetsU[wedges[j].v1] + wedges[j].j)], num);
	  writeAdd(&butterflies[eltsPerCacheLine*(offsetsV[wedges[j].u] + wedges[j].k)], num);
	});
    }
  }

  free(freq_arr);
  return wedges_pair.second;
}

/*
 *  Computes butterfly counts per edge, given a ranked graph and using
 *  sort aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountESortCE(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges,
                  long* wedge_idxs, intT curr_idx=0) {
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx,
                                                    num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, sort wedges to aggregate
#ifdef OPENMP
  sampleSort(wedges, num_wedges_f, UWedgeCmp());
#else
  pbbs::sample_sort(wedges, num_wedges_f, UWedgeCmp());
#endif
  // Retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  using X = tuple<uintE, long>;
  if (cs.butterflies_seq_intt.n < 2*num_wedges_f) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*num_wedges_f);
    cs.butterflies_seq_intt.n = 2*num_wedges_f;
  }

  // Store butterfly counts per edge in another array
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];
    if (num_butterflies > 1){
      // Iterate through all wedges
      parallel_for(long j=freq_arr[i]; j<freq_arr[i+1]; ++j) {
	// Store butterfly counts per edge (v_offset + k and u_offset + j respectively)
	cs.butterflies_seq_intt.A[2*j] = make_tuple((GA.offsets[wedges[j].u]+wedges[j].k), num_butterflies - 1);
	cs.butterflies_seq_intt.A[2*j+1] = make_tuple((GA.offsets[wedges[j].v1]+wedges[j].j), num_butterflies - 1);
      }
    }
    else {
      parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
	cs.butterflies_seq_intt.A[2*j] = make_tuple(UINT_E_MAX, 0);
	cs.butterflies_seq_intt.A[2*j+1] = make_tuple(UINT_E_MAX, 0);
      }
    }
  }

  free(freq_arr);

  // Aggregate butterfly counts
  // Sort butterfly counts to aggregate
#ifdef OPENMP
  sampleSort(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>());
#else
  pbbs::sample_sort(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>());
#endif
  // Retrieve a list of indices where consecutive counts have different keys
  auto b_freq_pair = getFreqs<long>(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>(),
                                    tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;

  // Store butterfly counts per edge
  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    if (get<0>(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]) != UINT_E_MAX) {
      long num_freq = b_freq_arr[i] - b_freq_arr[i-1];
      // Reduce to sum the butterflies over the necessary range
      X reduce = sequence::reduce(&(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]), num_freq, tupleAdd<uintE,long>());
      // These are our butterfly counts
      butterflies[eltsPerCacheLine*get<0>(reduce)] += get<1>(reduce);
    }
  }

  free(b_freq_arr);
  return wedges_pair.second;
}

/*
 *  Computes butterfly counts per edge, given a ranked graph and using
 *  sort aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountESort(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs,
                intT curr_idx=0) {
  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx,
                                                    num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, sort wedges to aggregate
#ifdef OPENMP
  sampleSort(wedges, num_wedges_curr, UWedgeCmp());
#else
  pbbs::sample_sort(wedges, num_wedges_curr, UWedgeCmp());
#endif
  // Retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());

  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Store butterfly counts per edge
  parallel_for (long i = 0; i < freq_pair.second-1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];
    if (num_butterflies > 1) {
      // Iterate through all wedges
      granular_for (j, freq_pair.first[i], freq_pair.first[i+1], (freq_pair.first[i+1]-freq_pair.first[i] > 1000), {
	  // Store butterfly counts per edge (v_offset + k and u_offset + j respectively)
	  writeAdd(&butterflies[eltsPerCacheLine*(GA.offsets[wedges[j].u] + wedges[j].k)], num_butterflies - 1);
	  writeAdd(&butterflies[eltsPerCacheLine*(GA.offsets[wedges[j].v1] + wedges[j].j)], num_butterflies - 1);
	});
    }
  }

  free(freq_pair.first);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per edge, given a ranked graph and using
 *  hash aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHash(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
                intT curr_idx=0) {
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Construct a hash table to aggregate wedges
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges,
                                wedge_idxs);

  // Retrieve aggregated wedge counts and store butterfly count per edge
  // Iterate through all wedges with endpoint i
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
#ifdef INVERSE
      {
#else
      if (v > i) {
#endif
	// Iterate through all wedges with endpoints i and u2, and center v
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2 = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
    if (u2 < i && u2 < v) {
#else
	  if (u2 > i) {
#endif
	    // Find the number of wedges with endpoints i and u2
	    long to_find = ((((long) i) *GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	    long num_butterflies = cs.wedges_hash.find(to_find).second;
	    if (num_butterflies > 1) {
	      // Store butterfly counts per edge, (i, v) and (u2, v)
	      writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], num_butterflies - 1);
	      writeAdd(&butterflies[eltsPerCacheLine*(u_offset+j)], num_butterflies - 1);
	    }
	  }
	  else break;
	}
      }
    }
  }
  return next_idx;
}

/*
 *  Computes butterfly counts per edge, given a ranked graph and using
 *  hash aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHashCE(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
                  intT curr_idx=0) {
  using T = pair<uintE, long>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Construct a hash table to aggregate wedges
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges,
                                wedge_idxs);

  // Retrieve aggregated wedge counts and store butterfly count per edge
  // Iterate through all wedges with endpoint i
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
#ifdef INVERSE
      {
#else
      if (v > i) {
#endif
	// Iterate through all wedges with endpoints i and u2, and center v
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2 = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
    if (u2 < i && u2 < v) {
#else
	  if (u2 > i) {
#endif
	    // Find the number of wedges with endpoints i and u2
	    long to_find = ((((long) i) * GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	    long num_butterflies = cs.wedges_hash.find(to_find).second;
	    if (num_butterflies > 1) {
	      // Insert butterfly counts per edge, (i, v) and (u2, v)
	      cs.butterflies_hash.insert(T((v_offset+k), num_butterflies-1));
	      cs.butterflies_hash.insert(T((u_offset+j), num_butterflies-1));
	    }
	  }
	  else break;
	}
      }
    }
  }

  size_t num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.wedges_seq_intp);

  // Store butterfly counts per edge
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterfly_pair = cs.wedges_seq_intp.A[i];
    butterflies[eltsPerCacheLine*butterfly_pair.first] += butterfly_pair.second;
  }
  return next_idx;
}

/*
 *  Computes butterfly counts per edge, using by side ranking and hash
 *  aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs           : Holds all array space needed, to be reused between batches
 *  GA           : Bipartite graph in CSR format
 *  use_v        : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                 produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges   : Number of wedges produced by the specified bipartition
 *  butterflies  : Array to store butterfly counts per edge, considering edges derived from one orientation (as given
 *                 by use_v)
 *  butterflies_u: Array to store butterfly counts per edge, considering edges derived from the other orientation (as 
 *                 given by use_v)
 *  max_wedges   : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                 memory
 *  wedges_idxs  : Wedge indices to allow for wedge retrieval in parallel
 *  eti          : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  curr_idx     : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHash(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long* butterflies_u,
                long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Construct a hash table to aggregate wedges
  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);

  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve aggregated wedge counts and store butterfly count per edge
  // Iterate through all wedges with endpoint i
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Iterate through all wedges with endpoints i and u2, and center v
      for (intT k=0; k < v_deg; ++k) {
	uintE u2 = edgesV[v_offset + k];
	if (u2 < i) {
	  // Find the number of wedges with endpoints i and u2
	  long num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
	  // Store butterfly counts per edge, (i, v) and (u2, v)
	  writeAdd(&butterflies_u[eltsPerCacheLine*(u_offset + j)],num_butterflies); 
	  writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)],num_butterflies); 
	}
	else break;
      }
    }
  }
  
  return next_idx;
}

/*
 *  Computes butterfly counts per edge, using by side ranking and hash
 *  aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *               produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per edge
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  eti        : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountEHashCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
                  long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Construct a hash table to aggregate wedges
  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);

  // Retrieve aggregated wedge counts and store butterfly count per edge in another hash table
  // Iterate through all wedges with endpoint i
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Iterate through all wedges with endpoints i and u2, and center v
      for (intT k=0; k < v_deg; ++k) {
	uintE u2 = edgesV[v_offset + k];
	if (u2 < i) {
	  // Find the number of wedges with endpoints i and u2
	  long num_butterflies = cs.wedges_hash.find(((long) i) * nu + (long) u2).second - 1;
	  // Insert butterfly counts per edge, (i, v) and (u2, v)
	  if (num_butterflies > 0) {
	    cs.butterflies_hash.insert(make_pair(eti[u_offset + j], num_butterflies));
	    cs.butterflies_hash.insert(make_pair(v_offset + k, num_butterflies));
	  }
	}
	else break;
      }
    }
  }

  size_t num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.wedges_seq_intp);
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Store butterfly counts per edge
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterflies_pair = cs.wedges_seq_intp.A[i];
    butterflies[eltsPerCacheLine*butterflies_pair.first] += butterflies_pair.second;
  }

  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per edge, using by side ranking and wedge-aware
 *  batching.
 * 
 *  GA             : Bipartite graph in CSR format
 *  butterflies    : Array to store butterfly counts per edge, considering edges derived from one orientation (as given
 *                   by use_v)
 *  butterflies_u  : Array to store butterfly counts per edge, considering edges derived from the other orientation (as 
 *                   given by use_v)
 *  use_v          : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                   produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest
 *                   wedges.
 *  max_array_size :
 *  wedgesPrefixSum: Wedge indices to identify which vertices to process per batch
 */
void CountEOrigCompactParallel_WedgeAware(bipartiteCSR& GA, long* butterflies, long* butterflies_u, bool use_v,
                                          long max_array_size, long* wedgesPrefixSum) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu); // tunable parameter

  timer t1,t2,t3,t4;
  t1.start();

  // Stores number of wedges per endpoint
  uintE* wedges = newA(uintE, nu*stepSize);
  // Keeps track of entries used in the wedges array
  uintE* used = newA(uintE, nu*stepSize);

  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();

  // Initialize wedge array
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  const size_t eltsPerCacheLine = 64/sizeof(long);

  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");
  t2.start();
#endif

  // Consider vertices i in batches, as given by wedgesPrefixSum
  for(intT step = 0; step < (nu+stepSize-1)/stepSize; step++) {
    std::function<void(intT,intT)> recursive_lambda =
      [&]
      (intT start, intT end){
      if ((start == end-1) || (wedgesPrefixSum[end]-wedgesPrefixSum[start] < 1000)){ 
	for (intT i = start; i < end; i++){
	  intT used_idx = 0;
	  long shift = nu*(i-step*stepSize);
	  intT u_offset  = offsetsU[i];
	  intT u_deg = offsetsU[i+1]-u_offset;
	  for (long j=0; j < u_deg; ++j ) {
	    uintE v = edgesU[u_offset+j];
	    intT v_offset = offsetsV[v];
	    intT v_deg = offsetsV[v+1]-offsetsV[v];
	    // Iterate through all 2-hop neighbors of i (with appropriate rank)
	    for (long k=0; k < v_deg; ++k) { 
	      uintE u2_idx = edgesV[v_offset+k];
	      if (u2_idx < i) {
		// Count wedges on the second endpoint
		wedges[shift+u2_idx]++;
		// Keep track of used second endpoints
		if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	      }
	      else break;
	    }
	  }

	  // Iterate through wedges to store butterfly counts on the edges
	  for(long j=0; j < u_deg; ++j) {
	    uintE v = edgesU[u_offset+j];
	    intT v_offset = offsetsV[v];
	    intT v_deg = offsetsV[v+1] - v_offset;
	    for(intT k=0; k < v_deg; ++k) {
	      uintE u2_idx = edgesV[v_offset+k];
	      if (u2_idx < i) {
		writeAdd(&butterflies_u[eltsPerCacheLine*(u_offset + j)], (long) wedges[shift+u2_idx]-1);
		writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)], (long) wedges[shift+u2_idx]-1);
	      }
	      else break;
	    }
	  }

	  // Clear wedges array for reuse (only need to clear on used second
	  // endpoints)
	  for(long j=0; j < used_idx; ++j) { wedges[used[shift+j]+shift] = 0; }
	}
      } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,nu));
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
}

/*
 *  Computes butterfly counts per edge, using by side ranking and simple
 *  batching.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  butterflies   : Array to store butterfly counts per edge, considering edges derived from one orientation (as given
 *                  by use_v)
 *  butterflies_u : Array to store butterfly counts per edge, considering edges derived from the other orientation (as 
 *                  given by use_v)
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_array_size:
 */
void CountEOrigCompactParallel(uintE* eti, long* butterflies, long* butterflies_u, bipartiteCSR& GA, bool use_v,
                               long max_array_size) {
  timer t1,t2,t3;
  t1.start();

  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = min<long>(getWorkers() * 15, max_array_size/nu); // tunable parameter

  // Stores number of wedges per endpoints
  uintE* wedges = newA(uintE, nu*stepSize);
  // Keeps track of entries used in the wedges array
  uintE* used = newA(uintE, nu*stepSize);

  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();

  // Initialize wedges array
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  const size_t eltsPerCacheLine = 64/sizeof(long);

  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");
  t2.start();
#endif

  // Consider vertices i in batches
  for(long step = 0; step < (nu+stepSize-1)/stepSize; step++) {
    parallel_for_1(long i=step*stepSize; i < min((step+1)*stepSize,nu); ++i){
      intT used_idx = 0;
      long shift = nu*(i-step*stepSize);
      intT u_offset  = offsetsU[i];
      intT u_deg = offsetsU[i+1]-u_offset;
      for (long j=0; j < u_deg; ++j ) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1]-offsetsV[v];
	// Iterate through all 2-hop neighbors of i (with appropriate rank)
	for (long k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx < i) {
	    // Count wedges on the second endpoint
	    wedges[shift+u2_idx]++;
	    // Keep track of used second endpoints
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	  }
	  else break;
	}
      }

      // Iterate through wedges to store butterfly counts on the edges
      for(long j=0; j < u_deg; ++j) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1] - v_offset;
	for(intT k=0; k < v_deg; ++k) {
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx < i) {
	    writeAdd(&butterflies_u[eltsPerCacheLine*(u_offset + j)], (long) wedges[shift+u2_idx]-1);
	    writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)], (long) wedges[shift+u2_idx]-1);
	  }
	  else break;
	}
      }

      // Clear wedges array for reuse (only need to clear on used second
      // endpoints)
      for(long j=0; j < used_idx; ++j) { wedges[used[shift+j]+shift] = 0; }
    }
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
}

/*
 *  Computes butterfly counts per edge, given a ranked graph and using
 *  wedge-aware batching.
 * 
 *  GA             : Ranked graph in CSR format
 *  butterflies    : Array to store butterfly counts per edge
 *  max_array_size :
 *  wedgesPrefixSum: Wedge indices to identify which vertices to process per batch
 */
void CountEOrigCompactParallel_WedgeAware(graphCSR& GA, long* butterflies, long max_array_size, long* wedgesPrefixSum) {
  long stepSize = min<long>(getWorkers() * 20, max_array_size/GA.n);

  timer t1,t2,t3,t4;
  t1.start();

  // Store number of wedges per endpoints
  uintE* wedges = newA(uintE, GA.n*stepSize);
  // Keeps track of entries used in the wedges array
  uintE* used = newA(uintE, GA.n*stepSize);

  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();

  // Initialize wedges array
  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  const size_t eltsPerCacheLine = 64/sizeof(long);

  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");
  t2.start();
#endif

  // Consider vertices i in batches, as given by wedgesPrefixSum
  for(intT step = 0; step < (GA.n+stepSize-1)/stepSize; step++) {
    std::function<void(intT,intT)> recursive_lambda =
      [&]
      (intT start, intT end){
      if ((start == end-1) || (wedgesPrefixSum[end]-wedgesPrefixSum[start] < 1000)){ 
	for (intT i = start; i < end; i++){
	  intT used_idx = 0;
	  long shift = GA.n*(i-step*stepSize);
	  intT u_offset  = GA.offsets[i];
	  intT u_deg = GA.offsets[i+1]-u_offset;
	  for (intT j=0; j < u_deg; ++j ) {
	    uintE v = GA.edges[u_offset+j] >> 1;
	    intT v_offset = GA.offsets[v];
	    intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
	    if (v <= i) break;
#endif
	    // Iterate through all 2-hop neighbors of i (with appropriate rank)
	    for (intT k=0; k < v_deg; ++k) { 
	      uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
        if (u2_idx < i && u2_idx < v) {
#else
	      if (u2_idx > i) {
#endif
		// Count wedges on the second endpoint
		wedges[shift+u2_idx]++;
		// Keep track of used second endpoints
		if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = GA.edges[v_offset+k];
	      }
	      else break;
	    }
	  }

	  // Iterate through wedges to store butterfly counts on the edges
	  for (long j=0; j < u_deg; ++j ) {
	    uintE v = GA.edges[u_offset+j] >> 1;
	    intT v_offset = GA.offsets[v];
	    intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
	    if (v <= i) break;
#endif
	    for (long k=0; k < v_deg; ++k) { 
	      uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
        if (u2_idx < i && u2_idx < v) {
#else
	      if (u2_idx > i) {
#endif
		if (wedges[shift+u2_idx] > 1) {
		  writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], (long)(wedges[shift+u2_idx]-1));
		  writeAdd(&butterflies[eltsPerCacheLine*(u_offset+j)], (long)(wedges[shift+u2_idx]-1));
		}
	      }
	      else break;
	    }
	  }

	  // Clear wedges array for reuse (only need to clear on used second
	  // endpoints)
	  for(long j=0; j < used_idx; ++j) {
	    uintE u2_idx = used[shift+j] >> 1;
	    wedges[shift+u2_idx] = 0;
	  }
	}
      } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,GA.n));
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
}

/*
 *  Computes butterfly counts per edge, given a ranked graph and using simple
 *  batching.
 * 
 *  GA            : Ranked graph in CSR format
 *  butterflies   : Array to store butterfly counts per edge
 *  max_array_size:
 */
void CountEWorkEfficientParallel(graphCSR& GA, long* butterflies, long max_array_size) {
  timer t1,t2,t3;
  t1.start();

  long stepSize = min<long>(getWorkers() * 7.5, max_array_size/GA.n);
  // Stores number of wedges per endpoints
  uintE* wedges = newA(uintE, GA.n*stepSize);
  // Keeps track of entries used in the wedges array
  uintE* used = newA(uintE, GA.n*stepSize);

  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();

  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  const size_t eltsPerCacheLine = 64/sizeof(long);

  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");
  t2.start();
#endif

  // Consider vertices i in batches
  for(long step = 0; step < (GA.n+stepSize-1)/stepSize; step++) {
    parallel_for_1(long i=step*stepSize; i < min((step+1)*stepSize,GA.n); ++i){
      intT used_idx = 0;
      long shift = GA.n*(i-step*stepSize);
      intT u_offset  = GA.offsets[i];
      intT u_deg = GA.offsets[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
	if (v <= i) break;
#endif
	// Iterate through all 2-hop neighbors of i (with appropriate rank)
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
    if (u2_idx < i && u2_idx < v) {
#else
	  if (u2_idx > i) {
#endif
	    // Count wedges on the second endpoint
	    wedges[shift+u2_idx]++;
	    // Keep track of used second endpoints
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = GA.edges[v_offset+k];
	  }
	  else break;
	}
      }

      // Iterate through wedges to store butterfly counts on the edges
      for (long j=0; j < u_deg; ++j ) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
	if (v <= i) break;
#endif
	for (long k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
    if (u2_idx < i && u2_idx < v) {
#else
	  if (u2_idx > i) {
#endif
	    if (wedges[shift+u2_idx] > 1) {
	      writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], (long)(wedges[shift+u2_idx]-1));
	      writeAdd(&butterflies[eltsPerCacheLine*(u_offset+j)], (long)(wedges[shift+u2_idx]-1));
	    }
	  }
	  else break;
	}
      }

      // Clear wedges array for reuse (only need to clear on used second
      // endpoints)
      for(long j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j] >> 1;
        wedges[shift+u2_idx] = 0;
      }
    }
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
}

/*
 *  Computes butterfly counts per edge, given a ranked graph and using a
 *  serial work-efficient algorithm.
 * 
 *  GA         : Ranked graph in CSR format
 *  butterflies: Array to store butterfly counts per vertex
 */
void CountEWorkEfficientSerial(graphCSR& GA, long* butterflies) {
  timer t1,t2,t3;
  t1.start();

  // Stores number of wedges per endpoints
  uintE* wedges = newA(uintE, GA.n);
  // Keeps track of entries used in the wedges array
  uintE* used = newA(uintE, GA.n);
  const size_t eltsPerCacheLine = 64/sizeof(long);

  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();

  for(long i=0;i<GA.n;i++) { wedges[i] = 0; }

  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");
  t2.start();
#endif

  // Consider every vertex i
  for(long i=0; i < GA.n; ++i){
    intT used_idx = 0;
    intT u_offset  = GA.offsets[i];
    intT u_deg = GA.offsets[i+1]-u_offset;
    for (intT j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
      if (v <= i) break;
#endif
      // Iterate through all 2-hop neighbors of i (with appropriate rank)
      for (intT k=0; k < v_deg; ++k) { 
	uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
  if (u2_idx < i && u2_idx < v) {
#else
	if (u2_idx > i) {
#endif
	  // Count wedges on the second endpoint
	  wedges[u2_idx]++;
	  // Keep track of used second endpoints
	  if (wedges[u2_idx] == 1) used[used_idx++] = GA.edges[v_offset+k];
	}
	else break;
      }
    }

    // Iterate through wedges to store butterfly counts on the edges
    for (long j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
      if (v <= i) break;
#endif
      for (long k=0; k < v_deg; ++k) { 
	uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
  if (u2_idx < i && u2_idx < v) {
#else
	if (u2_idx > i) {
#endif
	  if (wedges[u2_idx] > 1) {
	    butterflies[eltsPerCacheLine*(v_offset+k)] += (long)(wedges[u2_idx]-1);
	    butterflies[eltsPerCacheLine*(u_offset+j)] += (long)(wedges[u2_idx]-1);
	  }
	}
	else break;
      }
    }

    // Clear wedges array for reuse (only need to clear on used second
    // endpoints)
    for(long j=0; j < used_idx; ++j) {
      uintE u2_idx = used[j] >> 1;
      wedges[u2_idx] = 0;
    }
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
}

/*
 *  Computes butterfly counts per edge, using by side ranking and a serial
 *  algorithm.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  butterflies   : Array to store butterfly counts per edge
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 */
void CountESerial(uintE* eti, long* butterflies, bipartiteCSR& GA, bool use_v) {
  timer t1,t2,t3;
  t1.start();

  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Stores number of wedges per endpoints
  uintE* wedges = newA(uintE, nu);
  // Keeps track of entries used in the wedges array
  uintE* used = newA(uintE, nu);
 
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();

  for(long i=0;i<nu;i++) { wedges[i] = 0; }
  const size_t eltsPerCacheLine = 64/sizeof(long);

  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");
  t2.start();
#endif

  // Consider every vertex i
  for(long i=0; i < nu; ++i){
    intT used_idx = 0;
    intT u_offset  = offsetsU[i];
    intT u_deg = offsetsU[i+1]-u_offset;
    for (long j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1]-offsetsV[v];
      // Iterate through all 2-hop neighbors of i
      for (long k=0; k < v_deg; ++k) { 
        uintE u2_idx = edgesV[v_offset+k];
        if (u2_idx < i) {
          // Count wedges on the second endpoint
          wedges[u2_idx]++;
          // Keep track of used second endpoints
          if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
        }
        else break;
      }
    }

    // Iterate through centers of wedges to store butterfly counts on the edges
    for(long j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      for(intT k=0; k < v_deg; ++k) {
        uintE u2_idx = edgesV[v_offset+k];
        // Store butterfly counts per edge
        if (u2_idx < i) {
          butterflies[eltsPerCacheLine*eti[u_offset + j]] += wedges[u2_idx] - 1;
          butterflies[eltsPerCacheLine*(v_offset + k)] += wedges[u2_idx] - 1;
        }
        else break;
      }
    }

    // Clear wedges array for reuse (only need to clear on used second
    // endpoints)
    for(long j=0; j < used_idx; ++j) { wedges[used[j]] = 0; }
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
}

/*
 *  Computes butterfly counts per edge, given ranked vertices.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges    : Number of wedges produced by the given ranking
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  max_array_size:
 *  type          : Wedge/butterfly aggregation type
 *  ranks         : Sorted U and V vertices by rank (where U indices all come before V indices)
 *  rankV         : Array mapping V indices to their rank
 *  rankU         : Array mapping U indices to their rank
 * 
 *  Returns an array of butterfly counts, indexed by edge.
 *  Note: This function takes ownership of and frees ranks, rankV, and rankU
 */
long* CountERank(uintE* eti, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size,
                 CountType type, uintE* ranks, uintE* rankV, uintE* rankU) {
#ifdef VERBOSE
  timer t_rank;
  t_rank.start();
#endif

  // Rank graph using ranks, rankV, and rankU
  // Vertices in g are indexed by rank
  // Also, obtain a converter that maps ranked edge indices to edge indices in
  // the original graph
  auto g_pair = rankGraphEdges(GA, use_v, ranks, rankV, rankU);
  auto g = g_pair.first;
  auto edge_converter = g_pair.second;
  free(ranks); free(rankU); free(rankV); 

#ifdef VERBOSE
  t_rank.reportTotal("ranking");
  timer t_time, t_time2;
  t_time.start();
#endif

  // Initialize array to hold butterfly counts per edge, as indexed by rank
  const size_t eltsPerCacheLine = 64/sizeof(long);
  long numEdges = g.offsets[g.n];
  long* rank_butterflies = newA(long,eltsPerCacheLine*numEdges);

#ifdef VERBOSE
  t_time.reportTotal("preprocess (malloc)");
  t_time2.start();
#endif
  granular_for(i,0,numEdges,numEdges > 10000, { rank_butterflies[eltsPerCacheLine*i] = 0; });

  // Compute wedge indices so that wedges can be stored in parallel (for any
  // aggregation type except simple batching)
  long* wedge_idxs = (type == BATCHS || type == SERIAL) ? nullptr : countWedgesScan(g);

  // Initialize structure that holds all space for counting algorithms, to be
  // reused with wedge batches
  CountESpace cs = CountESpace(type, numEdges, true);

#ifdef VERBOSE
  if (type == BATCHWA) t_time2.reportTotal("counting (scan)");
#endif

  // Choose wedge/butterfly aggregation type
  if (type == BATCHS) CountEWorkEfficientParallel(g, rank_butterflies, max_array_size);
  else if (type == SERIAL) CountEWorkEfficientSerial(g, rank_butterflies);
  else if (type == BATCHWA) CountEOrigCompactParallel_WedgeAware(g, rank_butterflies, max_array_size, wedge_idxs);
  else {
    if (max_wedges >= num_wedges) {
      if (type == ASORT) CountESort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == SORT) CountESortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == AHASH) CountEHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == HASH) CountEHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == AHIST) CountEHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == HIST) CountEHistCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
    }
    else {
      intT curr_idx = 0;
      while(curr_idx < g.n) {
	if (type == ASORT) curr_idx = CountESort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == SORT) curr_idx = CountESortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == AHASH) curr_idx = CountEHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == HASH) curr_idx = CountEHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == AHIST) curr_idx = CountEHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == HIST) curr_idx = CountEHistCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	cs.clear();
      }
    }
  }

  cs.del();
  if (type != BATCHS && type != SERIAL) free(wedge_idxs);

#ifdef VERBOSE
  if (type != BATCHS && type != BATCHWA && type != SERIAL) t_time2.reportTotal("counting");
  timer t_convert;
  t_convert.start();
#endif

  // Initialize array to hold butterfly counts per vertex, by original indices
  long* butterflies = newA(long, eltsPerCacheLine*GA.numEdges);
  granular_for(i,0,GA.numEdges,GA.numEdges > 10000, { butterflies[eltsPerCacheLine * i] = 0; });

  // Convert all butterfly counts to be indexed by original indices
  parallel_for(intT i=0; i < g.n; ++i) {
    intT v_offset = g.offsets[i];
    intT v_deg = g.offsets[i+1] - v_offset;
    granular_for(j, 0, v_deg, v_deg > 1000, {
	uintE u = g.edges[v_offset + j] >> 1;
	if (g.edges[v_offset + j] & 0b1) writeAdd(&butterflies[eltsPerCacheLine*(get<1>(edge_converter[v_offset + j]))],
						  rank_butterflies[eltsPerCacheLine*(v_offset + j)]);
	else writeAdd(&butterflies[eltsPerCacheLine*eti[(get<1>(edge_converter[v_offset + j]))]],
		      rank_butterflies[eltsPerCacheLine*(v_offset + j)]);
      });
  }
  free(rank_butterflies);
  free(edge_converter);
  g.del();

#ifdef VERBOSE
  t_convert.reportTotal("convert");
#endif

  return butterflies;
}

/*
 *  Computes butterfly counts per edge.
 * 
 *  eti          : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges    : Number of wedges produced by the bipartition given by use_v, assuming ranking by side
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  max_array_size:
 *  type          : Wedge/butterfly aggregation type
 *  tw            : Vertex ordering option
 * 
 *  Returns an array of butterfly counts, indexed by edge.
 */
long* CountE(uintE* eti, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size,
             CountType type, RankType tw) {
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Any ranking that is not by side
  if (tw != SIDE) {
#ifdef VERBOSE
    timer t_rank;
    t_rank.start();
#endif

    // Pick the correct ranking
    tuple<uintE*,uintE*,uintE*> rank_tup;
    if (tw == COCORE) rank_tup = getCoCoreRanks(GA);
    else if (tw == ACOCORE) rank_tup = getApproxCoCoreRanks(GA);
    else if (tw == DEG) rank_tup = getDegRanks(GA);
    else if (tw == ADEG) rank_tup = getApproxDegRanks(GA);

    // Compute the number of wedges that the ranking produces
    long num_ccwedges = sequence::reduce<long>(
					       (long) 0, GA.nu, addF<long>(), rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup), get<1>(rank_tup)));
    num_ccwedges += sequence::reduce<long>(
					   (long) 0, GA.nv, addF<long>(), rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup)));

#ifdef VERBOSE
    t_rank.reportTotal("ranking");
#endif

    // Compute butterfly counts per edge using the specified ranking
    return CountERank(eti, GA, use_v, num_ccwedges, max_wedges, max_array_size, type, get<0>(rank_tup),
                      get<1>(rank_tup), get<2>(rank_tup));
  }

  // Ranking by side
#ifdef VERBOSE
  timer t_malloc;
  t_malloc.start();
#endif

  long* butterflies = newA(long, eltsPerCacheLine*GA.numEdges);
  parallel_for(long i=0; i < GA.numEdges; ++i){
    butterflies[eltsPerCacheLine*i] = 0;
  }

  // For storing butterfly counts without translating edges accessed on
  // different bipartitions
  long* butterflies_u = (type == HASH || type == SORT || type == HIST || type == SERIAL) ? nullptr :
    newA(long, eltsPerCacheLine*GA.numEdges);
  if (type != HASH && type != SORT && type != HIST && type != SERIAL) {
    parallel_for(long i=0; i < GA.numEdges; ++i){
      butterflies_u[eltsPerCacheLine*i] = 0;
    }
  }

#ifdef VERBOSE
  t_malloc.reportTotal("preprocess (malloc)");
  timer t_time;
  t_time.start();
#endif

  // Compute wedge indices so that wedges can be stored in parallel (for any
  // aggregation type except simple batching)
  long* wedge_idxs = type == BATCHS ? nullptr : countWedgesScan(GA, use_v, true);

  // Initialize structure that holds all space for counting algorithms, to be
  // reused with wedge batches
  CountESpace cs = CountESpace(type, GA.numEdges, false);

#ifdef VERBOSE
  if (type == BATCHWA) t_time.reportTotal("counting (scan)");
#endif

  // Choose wedge/butterfly aggregation type
  if (type == BATCHS) CountEOrigCompactParallel(eti, butterflies, butterflies_u, GA, use_v, max_array_size);
  else if (type == SERIAL) CountESerial(eti, butterflies, GA, use_v);
  else if (type == BATCHWA) CountEOrigCompactParallel_WedgeAware(GA, butterflies, butterflies_u, use_v, max_array_size,
                                                                 wedge_idxs);
  else {
    if (max_wedges >= num_wedges) {
      if (type == ASORT) CountESort(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti);
      else if (type == SORT) CountESortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
      else if (type == AHASH) CountEHash(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs,
                                         eti);
      else if (type == HASH) CountEHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
      else if (type == AHIST) CountEHist(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs,
                                         eti);
      else if (type == HIST) CountEHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    }
    else {
      intT curr_idx = 0;
      while(curr_idx < nu) {
	if (type == ASORT)
	  curr_idx = CountESort(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == SORT)
	  curr_idx = CountESortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == AHASH)
	  curr_idx = CountEHash(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == HASH)
	  curr_idx = CountEHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == AHIST)
	  curr_idx = CountEHist(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == HIST)
	  curr_idx = CountEHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
	cs.clear();
      }
    }
  }
  cs.del();

  if (type != BATCHS && type != SERIAL) free(wedge_idxs);
#ifdef VERBOSE
  if (type != SERIAL && type != BATCHWA && type != SERIAL) t_time.reportTotal("counting");
#endif

  if (type != HASH && type != SORT && type != HIST && type != SERIAL) {
#ifdef VERBOSE
    timer t_convert;
    t_convert.start();
#endif
  
    // Collate butterflies if counts stored in two arrays
    parallel_for(long i=0; i < GA.numEdges; ++i) {
      butterflies[eltsPerCacheLine*eti[i]] += butterflies_u[eltsPerCacheLine*i];
    }
    free(butterflies_u);
  
#ifdef VERBOSE
    t_convert.reportTotal("convert");
#endif
  }

  return butterflies;
}

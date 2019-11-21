#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
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

#include "butterfly_utils.h"

using namespace std;

/*
 *  Computes butterfly counts per vertex, using by side ranking, sort
 *  aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *               wedges. If true, stores counts on U. If false, stores counts on V.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountSort(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
               long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair = getWedges<UVertexPair>(cs.wedges_seq_uvp, GA, use_v, UVertexPairCons(), max_wedges,
                                                        curr_idx, num_wedges, wedge_idxs);
  UVertexPair* wedges = cs.wedges_seq_uvp.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, sort wedges to aggregate
#ifdef OPENMP
  sampleSort(wedges, num_wedges_curr, UVertexPairCmp());
#else
  pbbs::sample_sort(wedges, num_wedges_curr, UVertexPairCmp());
#endif
  // Retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UVertexPairCmp(), UVertexPairEq(), LONG_MAX, nonMaxLongF());
  auto wedge_freqs_f = freq_pair.first;
  auto num_wedge_freqs_f = freq_pair.second;

  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Store butterfly counts per vertex
  parallel_for (long i = 0; i < num_wedge_freqs_f-1; ++i) {
    long num_butterflies = wedge_freqs_f[i+1] - wedge_freqs_f[i];
    long wedge_idx = wedge_freqs_f[i];
    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v1],num_butterflies); 
    writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v2],num_butterflies); 
  }

  free(freq_pair.first);
  return wedges_pair.second;
}

/*
 *  Computes butterfly counts per vertex, using by side ranking and sort
 *  aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *               wedges. If true, stores counts on U. If false, stores counts on V.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountSortCE(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
                 long* wedge_idxs, intT curr_idx=0) {
  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair = getWedges<UVertexPair>(cs.wedges_seq_uvp, GA, use_v, UVertexPairCons(), max_wedges,
                                                        curr_idx, num_wedges, wedge_idxs);
  UVertexPair* wedges = cs.wedges_seq_uvp.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, sort wedges to aggregate
#ifdef OPENMP
  sampleSort(wedges, num_wedges_f, UVertexPairCmp());
#else
  pbbs::sample_sort(wedges, num_wedges_f, UVertexPairCmp());
#endif
  // Retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UVertexPairCmp(), UVertexPairEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  using X = tuple<uintE, long>;
  if (cs.butterflies_seq_intt.n < 2*freq_pair.second-2) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*freq_pair.second-2);
    cs.butterflies_seq_intt.n = 2*freq_pair.second-2;
  }

  // Store butterfly counts per vertex in another array
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num = freq_arr[i+1] - freq_arr[i];
    num = num * (num-1) / 2;
    long j = freq_arr[i];
    cs.butterflies_seq_intt.A[2*i] = make_tuple(wedges[j].v1, num);
    cs.butterflies_seq_intt.A[2*i+1] = make_tuple(wedges[j].v2, num);
  }

  free(freq_arr);

  // Aggregate butterfly counts
  // Sort butterfly counts to aggregate
#ifdef OPENMP
  sampleSort(cs.butterflies_seq_intt.A, 2*freq_pair.second-2, tupleLt<uintE,long>());
#else
  pbbs::sample_sort(cs.butterflies_seq_intt.A, 2*freq_pair.second-2, tupleLt<uintE,long>());
#endif
  // Retrieve a list of indices where consecutive counts have different keys
  pair<long*, long> b_freq_pair = getFreqs<long>(cs.butterflies_seq_intt.A, 2*freq_pair.second-2, tupleLt<uintE,long>(),
                                                 tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Store butterfly counts per vertex
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
 *  Computes butterfly counts per vertex, given a ranked graph and using
 *  sort aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountSortCE(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs,
                 intT curr_idx=0) {
  const size_t eltsPerCacheLine = 64/sizeof(long);
  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx, num_wedges,
                                                   wedge_idxs);
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

  using X = tuple<uintE,long>;
  if (cs.butterflies_seq_intt.n < num_wedges_f) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, num_wedges_f);
    cs.butterflies_seq_intt.n = num_wedges_f;
  }

  // Store butterfly counts per vertex in another array
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];
    // Store butterfly counts per endpoint, if it is on the right bipartition
    if ((wedges[wedge_idx].v2 & 0b1) && num_butterflies > 1) {
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      cs.butterflies_seq_intt.A[freq_arr[i]] = make_tuple(wedges[wedge_idx].v1, num_butterflies);
      cs.butterflies_seq_intt.A[freq_arr[i]+1] = make_tuple(wedges[wedge_idx].v2 >> 1, num_butterflies);
      granular_for(j,freq_arr[i]+2,freq_arr[i+1],freq_arr[i+1]-(freq_arr[i]+2) > 10000, {
	  cs.butterflies_seq_intt.A[j] = make_tuple(UINT_E_MAX, 0);
	});
    }
    // Store butterfly counts per center, if it is on the right bipartition
    else if (!(wedges[wedge_idx].v2 & 0b1) && num_butterflies > 1){
      granular_for(j,freq_arr[i],freq_arr[i+1],freq_arr[i+1]-freq_arr[i] > 10000, {
	  cs.butterflies_seq_intt.A[j] = make_tuple(wedges[j].u, num_butterflies - 1);
	});
    }
    // Fill in the rest of the array
    else {
      granular_for(j,freq_arr[i],freq_arr[i+1],freq_arr[i+1]-freq_arr[i] > 10000, {
	  cs.butterflies_seq_intt.A[j] = make_tuple(UINT_E_MAX, 0);
	});
    }
  }
  free(freq_arr);

  // Aggregate butterfly counts
  // Sort butterfly counts to aggregate
#ifdef OPENMP
  sampleSort(cs.butterflies_seq_intt.A, num_wedges_f, tupleLt<uintE,long>());
#else
  pbbs::sample_sort(cs.butterflies_seq_intt.A, num_wedges_f, tupleLt<uintE,long>());
#endif
  // Retrieve a list of indices where consecutive counts have different keys
  auto b_freq_pair = getFreqs<long>(cs.butterflies_seq_intt.A, num_wedges_f, tupleLt<uintE,long>(),
                                    tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;

  // Store butterfly counts per vertex
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
 *  Computes butterfly counts per vertex, given a ranked graph and using
 *  sort aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountSort(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs,
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
  // Store butterfly counts per vertex
  parallel_for (long i = 0; i < freq_pair.second-1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];
    // Store butterfly counts per endpoint, if it is on the right bipartition
    if ((wedges[wedge_idx].v2 & 0b1) && (num_butterflies > 1)) {
      num_butterflies = (num_butterflies * (num_butterflies - 1))/2;
      writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v1],num_butterflies); 
      writeAdd(&butterflies[eltsPerCacheLine*(wedges[wedge_idx].v2 >> 1)],num_butterflies);
    }
    // Store butterfly counts per center, if it is on the right bipartition
    else if (!(wedges[wedge_idx].v2 & 0b1) && (num_butterflies > 1)) {
      granular_for(j,freq_pair.first[i],freq_pair.first[i+1],freq_pair.first[i+1]-freq_pair.first[i] > 10000, {
	  writeAdd(&butterflies[eltsPerCacheLine*wedges[j].u], num_butterflies - 1);
	});
    }
  }

  free(freq_pair.first);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per vertex, given a ranked graph and using
 *  hash aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHash(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
               intT curr_idx=0) {
  const size_t eltsPerCacheLine = 64/sizeof(long);
  // Construct a hash table to aggregate wedges
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges,
                                wedge_idxs);
  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  // Retrieve aggregated wedge counts and store butterfly count per endpoint, if it is on the right bipartition
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    if (num_butterflies > 1 && (wedge_freq_pair.first & 0b1)) {
      uintE u2 = ((wedge_freq_pair.first >> 1) % (GA.n));
      uintE u = ((wedge_freq_pair.first >> 1) / (GA.n));
      num_butterflies = (num_butterflies * (num_butterflies - 1))/2;
      writeAdd(&butterflies[eltsPerCacheLine*u], num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*u2], num_butterflies);
    }
  }

  // Retrieve aggregated wedge counts and store butterfly count per center, if it is on the right bipartition
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    if (u_deg > 0 && (GA.edges[u_offset] & 0b1)) {
      granular_for(j,0,u_deg,u_deg>1000,{
	  uintE v = GA.edges[u_offset+j] >> 1;
	  intT v_offset = GA.offsets[v];
	  intT v_deg = GA.offsets[v+1] - v_offset;
	  if (v > i) {
	    for (intT k=0; k < v_deg; ++k) { 
	      uintE u2 = GA.edges[v_offset+k] >> 1;
	      if (u2 > i) {
		// Retrieve the wedge count corresponding to the wedge
		long to_find = (((long) i) *GA.n + (long) u2) << 1;
		long num_butterflies = cs.wedges_hash.find(to_find).second;
		if (num_butterflies > 1) writeAdd(&butterflies[eltsPerCacheLine*v], num_butterflies - 1);
	      }
	      else break;
	    }
	  }
	});
    }
  }

  return next_idx;
}

/*
 *  Computes butterfly counts per vertex, given a ranked graph and using
 *  hash aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHashCE(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
                 intT curr_idx=0) {
  using T = pair<uintE,long>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Construct a hash table to aggregate wedges
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges,
                                wedge_idxs);
  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  // Retrieve aggregated wedge counts and store butterfly count per endpoint, if it is on the right bipartition
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    if (num_butterflies > 1 && (wedge_freq_pair.first & 0b1)) {
      uintE u2 = (wedge_freq_pair.first >> 1) % (GA.n);
      uintE u = (wedge_freq_pair.first >> 1) / (GA.n);
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      cs.butterflies_hash.insert(T(u, num_butterflies));
      cs.butterflies_hash.insert(T(u2, num_butterflies));
    } 
  }

  // Retrieve aggregated wedge counts and store butterfly count per center, if it is on the right bipartition
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    if (u_deg > 0 && (GA.edges[u_offset] & 0b1)) {
      granular_for(j,0,u_deg,u_deg > 10000, { 
	  uintE v = GA.edges[u_offset+j] >> 1;
	  intT v_offset = GA.offsets[v];
	  intT v_deg = GA.offsets[v+1] - v_offset;
	  if (v > i) {
	    for (intT k=0; k < v_deg; ++k) { 
	      uintE u2 = GA.edges[v_offset+k] >> 1;
	      if (u2 > i) {
		// Retrieve the wedge count corresponding to the wedge
		long to_find = ((long)i * GA.n + u2) << 1;
		long num_butterflies = cs.wedges_hash.find(to_find).second;
		if (num_butterflies > 1) cs.butterflies_hash.insert(T(v, num_butterflies - 1));
	      }
	      else break;
	    }
	  }
	});
    }
  }
  num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.butterflies_seq_intp);

  // Store butterfly counts per vertex
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterfly_pair = cs.butterflies_seq_intp.A[i];
    butterflies[eltsPerCacheLine * butterfly_pair.first] += butterfly_pair.second;
  }

  return next_idx;
}

/*
 *  Computes butterfly counts per vertex, using by side ranking and hash
 *  aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *               wedges. If true, stores counts on U. If false, stores counts on V.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHash(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
               long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Construct a hash table to aggregate wedges
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges,
                                wedge_idxs);
  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  // Retrieve aggregated wedge counts and store butterfly count per vertex
  parallel_for (long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;

    if (num_butterflies > 1){
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      writeAdd(&butterflies[eltsPerCacheLine*v1],num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*v2],num_butterflies);
    }
  }

  return next_idx;
}

/*
 *  Computes butterfly counts per vertex, using by side ranking and hash
 *  aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *               wedges. If true, stores counts on U. If false, stores counts on V.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHashCE(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
                 long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Construct a hash table to aggregate wedges
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges,
                                wedge_idxs);
  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);
  using T = pair<uintE, long>;

  // Retrieve aggregated wedge counts and store butterfly count per vertex in another hash table
  parallel_for (long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;
    if (num_butterflies > 1) {
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      cs.butterflies_hash.insert(T(v1, num_butterflies));
      cs.butterflies_hash.insert(T(v2, num_butterflies));
    }
  }
  num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.butterflies_seq_intp);

  // Store butterfly counts per vertex
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterfly_pair = cs.butterflies_seq_intp.A[i];
    butterflies[eltsPerCacheLine*butterfly_pair.first] += butterfly_pair.second;
  }
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per vertex, given a ranked graph and using
 *  histogram aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHistCE(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
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
  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp,
                                                                       cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  // Set up a hash table to hold wedge counts (for storing counts in centers)
  cs.wedges_hash.resize(wedge_freqs_n);
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }
  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  // Set up an array to hold butterfly counts, so that they can be aggregated again
  if (cs.butterflies_seq_intt.n < num_wedges_list) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(T, num_wedges_list);
    cs.butterflies_seq_intt.n = num_wedges_list;
  }
  parallel_for(long i=0; i < num_wedges_list; ++i) {
    cs.butterflies_seq_intt.A[i] = make_tuple(UINT_E_MAX, 0);
  }

  // Store butterfly counts per endpoint, if it is on the right bipartition
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    if (num_butterflies > 1 && (wedge_freq_pair.first & 0b1)) {
      uintE u2 = ((wedge_freq_pair.first >> 1) % (GA.n));
      uintE u = ((wedge_freq_pair.first >> 1) / (GA.n));
      long wedge_idx = wedge_idxs[u] - wedge_idxs[curr_idx];
      num_butterflies = (num_butterflies * (num_butterflies - 1))/2;
      cs.butterflies_seq_intt.A[wedge_idx] = make_tuple(u, num_butterflies);
      cs.butterflies_seq_intt.A[wedge_idx+1] = make_tuple(u2, num_butterflies);
    }
  }

  // Store butterfly counts per center, if it is on the right bipartition
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    long idx = 0;
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    if (u_deg > 0 && (GA.edges[u_offset] & 0b1)) {
      for(intT j=0;j<u_deg;j++){
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1] - v_offset;
	if (v > i) {
	  for (intT k=0; k < v_deg; ++k) { 
	    uintE u2 = GA.edges[v_offset+k] >> 1;
	    if (u2 > i) {
	      // Retrieve the wedge count corresponding to the wedge
	      long to_find = (((long) i) *GA.n + (long) u2) << 1;
	      long num_butterflies = cs.wedges_hash.find(to_find).second;
	      if (num_butterflies > 1) {
		cs.butterflies_seq_intt.A[wedge_idx+idx] = make_tuple(v, num_butterflies-1);
	      }
	      ++idx;
	    }
	    else break;
	  }
	}
      }
    }
  }

  // Construct a histogram to aggregate butterfly counts
  pbbsa::sequence<T> wedge_freqs_i_seq = pbbsa::sequence<T>(cs.butterflies_seq_intt.A,num_wedges_list);
  tuple<size_t, T*> butterflies_tuple =
    pbbsa::sparse_histogram_f<uintE,long>(wedge_freqs_i_seq, GA.n, getAdd<uintE,long>, getAddReduce<uintE,long>,
                                          cs.tmp_uint, cs.out_uint);
  T* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  // Store butterfly counts per vertex
  parallel_for (long i=0; i < butterflies_n; ++i) {
    if (get<0>(butterflies_l[i]) != UINT_E_MAX)
      butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }
  return next_idx;
#endif
}

/*
 *  Computes butterfly counts per vertex, given a ranked graph and using
 *  histogram aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Ranked graph in CSR format
 *  num_wedges : Number of wedges produced by the ranking
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHist(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
               intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<long,uintE>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_pair = getWedges<long>(cs.wedges_seq_int, GA, UWedgeIntRankCons(GA.n), max_wedges, curr_idx,
                                                 num_wedges, wedge_idxs);
  intT next_idx = wedges_pair.second;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(cs.wedges_seq_int.A, wedges_pair.first);

  // Construct a histogram to aggregate wedges
  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp,
                                                                       cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  // Set up a hash table to hold wedge counts (for storing counts in centers)
  cs.wedges_hash.resize(wedge_freqs_n);
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }
  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  // Store butterfly counts per endpoint, if it is on the right bipartition
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    if (num_butterflies > 1 && (wedge_freq_pair.first & 0b1)) {
      uintE u2 = ((wedge_freq_pair.first >> 1) % (GA.n));
      uintE u = ((wedge_freq_pair.first >> 1) / (GA.n));
      num_butterflies = (num_butterflies * (num_butterflies - 1))/2;
      writeAdd(&butterflies[eltsPerCacheLine*u], num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*u2], num_butterflies);
    }
  }

  // Store butterfly counts per center, if it is on the right bipartition
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    if (u_deg > 0 && (GA.edges[u_offset] & 0b1)) {
      parallel_for(intT j=0;j<u_deg;j++){
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1] - v_offset;
	if (v > i) {
	  for (intT k=0; k < v_deg; ++k) { 
	    uintE u2 = GA.edges[v_offset+k] >> 1;
	    if (u2 > i) {
	      // Retrieve the wedge count corresponding to the wedge
	      long to_find = (((long) i) *GA.n + (long) u2) << 1;
	      long num_butterflies = cs.wedges_hash.find(to_find).second;
	      if (num_butterflies > 1) writeAdd(&butterflies[eltsPerCacheLine*v], num_butterflies - 1);
	    }
	    else break;
	  }
	}
      }
    }
  }

  return next_idx;
#endif
}

/*
 *  Computes butterfly counts per vertex, using by side ranking, histogram
 *  aggregation (for wedge aggregation), and atomic CAS (for butterfly
 *  aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *               wedges. If true, stores counts on U. If false, stores counts on V.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHist(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
               long* wedge_idxs, intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  using T = tuple<long,uintE>;

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_list_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, UVertexPairIntCons(nu), max_wedges,
                                                      curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_curr = wedges_list_pair.first;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_curr);

  // Construct a histogram to aggregate wedges
  tuple<size_t,T*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  T* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  // Retrieve aggregated wedge counts and store butterfly count per vertex using CAS
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    auto wedge_freq_pair = wedge_freqs[i];
    long num_butterflies = get<1>(wedge_freq_pair);
    auto wedge_freq_pair_first = get<0>(wedge_freq_pair);

    uintE v1 = wedge_freq_pair_first / nu;
    uintE v2 = wedge_freq_pair_first % nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    if (num_butterflies > 0) {
      writeAdd(&butterflies[eltsPerCacheLine*v1],num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*v2],num_butterflies);
    }
  }

  return wedges_list_pair.second;
#endif
}

/*
 *  Computes butterfly counts per vertex, using by side ranking and histogram
 *  aggregation (for wedge and butterfly aggregation).
 * 
 *  cs         : Holds all array space needed, to be reused between batches
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *               wedges. If true, stores counts on U. If false, stores counts on V.
 *  num_wedges : Number of wedges produced by the specified bipartition
 *  butterflies: Array to store butterfly counts per vertex
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  wedges_idxs: Wedge indices to allow for wedge retrieval in parallel
 *  curr_idx   : Denotes the vertex to start processing the batch from
 * 
 *  Returns an index specifying the next vertex that still needs to be processed
 */
intT CountHistCE(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
                 long* wedge_idxs, intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  const long nu = use_v ? GA.nu : GA.nv;
  const long nv = use_v ? GA.nv : GA.nu;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  using T = tuple<long, uintE>;
  using X = tuple<uintE, long>;

  // Retrieve all wedges in this batch
  pair<long, intT> wedges_list_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, UVertexPairIntCons(nu), max_wedges,
                                                      curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_curr = wedges_list_pair.first;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list, num_wedges_curr);

  // Construct a histogram to aggregate wedges
  tuple<size_t,T*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nv *  nu + nu, cs.tmp, cs.out);
  T* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  if (cs.butterflies_seq_intt.n < 2*wedge_freqs_n) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*wedge_freqs_n);
    cs.butterflies_seq_intt.n = 2*wedge_freqs_n;
  }

  // Retrieve aggregated wedge counts and store butterfly count per vertex
  parallel_for(long i = 0; i < wedge_freqs_n; i++) {
    auto wedge_freq_pair = wedge_freqs[i];
    long num = get<1>(wedge_freq_pair);
    num = (num * (num-1)) / 2;
    auto wedge_num = get<0>(wedge_freq_pair);
    cs.butterflies_seq_intt.A[2*i] = make_tuple((uintE) ((long) wedge_num % nu), (long) num);
    cs.butterflies_seq_intt.A[2*i + 1] = make_tuple((uintE) ((long) wedge_num / nu), (long) num);
  }

  // Construct a histogram to aggregate butterfly counts
  pbbsa::sequence<X> wedge_freqs_i_seq = pbbsa::sequence<X>(cs.butterflies_seq_intt.A,2*wedge_freqs_n);
  tuple<size_t, X*> butterflies_tuple =  pbbsa::sparse_histogram_f<uintE,long>(wedge_freqs_i_seq,nu,getAdd<uintE,long>,
                                                                               getAddReduce<uintE,long>, cs.tmp_uint,
                                                                               cs.out_uint);
  X* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  // Store butterfly counts per vertex
  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine * get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  return wedges_list_pair.second;
#endif
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per vertex, using by side ranking and wedge-aware
 *  batching.
 * 
 *  GA             : Bipartite graph in CSR format
 *  use_v          : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *                   wedges. If true, stores counts on U. If false, stores counts on V.
 *  butterflies    : Array to store butterfly counts per vertex
 *  max_array_size :
 *  wedgesPrefixSum: Wedge indices to identify which vertices to process per batch
 */
void CountOrigCompactParallel_WedgeAware(bipartiteCSR& GA, bool use_v, long* butterflies, long max_array_size,
                                         long* wedgesPrefixSum) {
  timer t1,t2,t3,t4;
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu); // tunable parameter

  t1.start();

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
	  for (intT j=0; j < u_deg; ++j ) {
	    uintE v = edgesU[u_offset+j];
	    intT v_offset = offsetsV[v];
	    intT v_deg = offsetsV[v+1]-offsetsV[v];
	    // Iterate through all 2-hop neighbors of i
	    for (intT k=0; k < v_deg; ++k) {
	      uintE u2_idx = edgesV[v_offset+k];
	      if (u2_idx < i) {
		// Increment the number of wedges on the second endpoint
		wedges[shift+u2_idx]++;
		// Keep track of used second endpoints
		if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	      }
	      else break;
	    }
	  }
	  // Iterate through all second endpoints to store butterfly counts on
	  for(intT j=0; j < used_idx; ++j) {
	    uintE u2_idx = used[shift+j];
	    writeAdd(&butterflies[eltsPerCacheLine*i],  (long)((long) wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
	    writeAdd(&butterflies[eltsPerCacheLine*u2_idx], (long)((long) wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
	    // Clear wedges array for reuse
	    wedges[u2_idx+shift] = 0;
	  }
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
 *  Computes butterfly counts per vertex, using by side ranking and simple
 *  batching.
 * 
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *                  wedges. If true, stores counts on U. If false, stores counts on V.
 *  butterflies   : Array to store butterfly counts per vertex
 *  max_array_size:
 */
void CountOrigCompactParallel(bipartiteCSR& GA, bool use_v, long* butterflies, long max_array_size) {
  timer t1,t2,t3;
  t1.start();
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu); // tunable parameter

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

  // Consider vertices i in batches, as given by stepSize
  for(intT step = 0; step < (nu+stepSize-1)/stepSize; step++) {
    parallel_for_1(intT i=step*stepSize; i < min((step+1)*stepSize,nu); ++i){
      intT used_idx = 0;
      long shift = nu*(i-step*stepSize);
      intT u_offset  = offsetsU[i];
      intT u_deg = offsetsU[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1]-offsetsV[v];
	// Iterate through all 2-hop neighbors of i
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx < i) {
	    // Increment the number of wedges on the second endpoint
	    wedges[shift+u2_idx]++;
	    // Keep track of used second endpoints
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	  }
	  else break;
	}
      }

      // Iterate through all second endpoints to store butterfly counts on
      for(intT j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j];
        writeAdd(&butterflies[eltsPerCacheLine*i], (long)((long) wedges[shift+u2_idx] * (wedges[shift+u2_idx]-1) / 2));
        writeAdd(&butterflies[eltsPerCacheLine*u2_idx],
                 (long)((long) wedges[shift+u2_idx] * (wedges[shift+u2_idx]-1) / 2));
        // Clear wedges array for reuse
        wedges[u2_idx+shift] = 0;
      }
    }
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif

  free(wedges);
  free(used);
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per vertex, given a ranked graph and using
 *  wedge-aware batching.
 * 
 *  GA             : Ranked graph in CSR format
 *  butterflies    : Array to store butterfly counts per vertex
 *  max_array_size :
 *  wedgesPrefixSum: Wedge indices to identify which vertices to process per batch
 */
void CountOrigCompactParallel_WedgeAware(graphCSR& GA, long* butterflies, long max_array_size, long* wedgesPrefixSum) {
  const size_t eltsPerCacheLine = 64/sizeof(long);
  long stepSize = min<long>(getWorkers() * 40, max_array_size/GA.n); // tunable parameter

  timer t1,t2,t3,t4;
  t1.start();

  // Stores number of wedges per endpoints
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
		// Increment the number of wedges on the second endpoint
		wedges[shift+u2_idx]++;
		// Keep track of used second endpoints
		if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = GA.edges[v_offset+k];
	      }
	      else break;
	    }
	  }

	  // Iterate through centers of wedges to store butterfly counts on
	  for (long j=0; j < u_deg; ++j ) {
	    uintE v = GA.edges[u_offset+j] >> 1;
	    intT v_offset = GA.offsets[v];
	    intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
	    if (v <= i) break; 
#endif
	    if (!(GA.edges[u_offset+j] & 0b1)) continue;
	    for (long k=0; k < v_deg; ++k) { 
	      uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
        if (u2_idx < i && u2_idx < v) {
#else
	      if (u2_idx > i) {
#endif
		// Only store butterfly counts per center if it is on the right
		// bipartition
		if (wedges[shift+u2_idx] > 1) writeAdd(&butterflies[eltsPerCacheLine*v], (long)(wedges[shift+u2_idx]-1));
	      }
	      else break;
	    }
	  }
      
	  // Iterate through all second endpoints to store butterfly counts on
	  for(long j=0; j < used_idx; ++j) {
	    uintE u2_idx = used[shift+j] >> 1;
	    // Only store butterfly counts per center if it is on the right bipartition
	    if(used[shift+j] & 0b1) {
	      writeAdd(&butterflies[eltsPerCacheLine*i], (long)((long)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
	      writeAdd(&butterflies[eltsPerCacheLine*u2_idx],
                 (long)((long)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
	    }
	    // Clear wedges array for reuse
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
 *  Computes butterfly counts per vertex, given a ranked graph and using simple
 *  batching.
 * 
 *  GA            : Ranked graph in CSR format
 *  butterflies   : Array to store butterfly counts per vertex
 *  max_array_size:
 */
void CountWorkEfficientParallel(graphCSR& GA, long* butterflies, long max_array_size) {
  timer t1,t2,t3;
  t1.start();

  long stepSize = min<long>(getWorkers() * 60, max_array_size/GA.n);

  // Stores number of wedges per endpoints
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
  
  // Consider vertices i in batches, as given by stepSize
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
	    // Increment the number of wedges on the second endpoint
	    wedges[shift+u2_idx]++;
	    // Keep track of used second endpoints
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = GA.edges[v_offset+k];
	  }
	  else break;
	}
      }

      // Iterate through centers of wedges to store butterfly counts on
      for (long j=0; j < u_deg; ++j ) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
	if (v <= i) break;
#endif
	if (!(GA.edges[u_offset+j] & 0b1)) continue;
	for (long k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
    if (u2_idx < i && u2_idx < v) {
#else
	  if (u2_idx > i) {
#endif
	    // Only store butterfly counts per center if it is on the right
	    // bipartition
	    if (wedges[shift+u2_idx] > 1) writeAdd(&butterflies[eltsPerCacheLine*v], (long)(wedges[shift+u2_idx]-1));
	  }
	  else break;
	}
      }
      
      // Iterate through all second endpoints to store butterfly counts on
      for(long j=0; j < used_idx; ++j) {
	uintE u2_idx = used[shift+j] >> 1;
	// Only store butterfly counts per endpoint if it is on the right bipartition
	if(used[shift+j] & 0b1) {
	  writeAdd(&butterflies[eltsPerCacheLine*i],  (long)((long)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
	  writeAdd(&butterflies[eltsPerCacheLine*u2_idx], (long)((long)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
	}
	// Clear wedges array for reuse
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

//********************************************************************************************
//********************************************************************************************


/*
 *  Computes butterfly counts per vertex, using side ranking and a
 *  work-efficient algorithm.
 * 
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition to store counts on, based on which bipartition produces the fewest
 *               wedges. If true, stores counts on U. If false, stores counts on V.
 *  butterflies: Array to store butterfly counts per vertex
 */
void CountOrigCompactSerial(bipartiteCSR& GA, bool use_v, long* butterflies) {
  timer t1,t2;
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

  for(long i=0; i < nu; ++i) { wedges[i] = 0; }

#ifdef VERBOSE
  t1.reportTotal("preprocess");
#endif
  t2.start();

  // Consider every vertex i
  for(intT i=0; i < nu; ++i){
    intT used_idx = 0;
    intT u_offset  = offsetsU[i];
    intT u_deg = offsetsU[i+1]-u_offset;
    for (intT j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1]-offsetsV[v];
      // Iterate through all 2-hop neighbors of i
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2_idx = edgesV[v_offset+k];
        if (u2_idx < i) {
          // Store butterfly counts per endpoint
          butterflies[i] += wedges[u2_idx];
          butterflies[u2_idx] += wedges[u2_idx];
          // Increment number of wedges on second endpoint
          wedges[u2_idx]++;
          // Keep track of used second endpoints
          if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
        }
        else break;
      }
    }

    // Clear wedges array for reuse (only need to clear on used second
    // endpoints)
    for(intT j=0; j < used_idx; ++j) { wedges[used[j]] = 0; }
  }

#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
}

/*
 *  Computes butterfly counts per vertex, given a ranked graph and using a
 *  serial work-efficient algorithm.
 * 
 *  GA         : Ranked graph in CSR format
 *  butterflies: Array to store butterfly counts per vertex
 */
void CountWorkEfficientSerial(graphCSR& GA, long* butterflies) {
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
	  // Only store butterfly counts per endpoint if it is on the right
	  // bipartition
	  if (GA.edges[v_offset+k] & 0b1) {
	    butterflies[eltsPerCacheLine*i] += wedges[u2_idx];
	    butterflies[eltsPerCacheLine*u2_idx] += wedges[u2_idx];
	  }
	  // Keep a running count x of wedges on the second endpoint, and adding
	  // this incrementally gives us our desired x choose 2 number of
	  // butterflies per endpoint
	  wedges[u2_idx]++;
	  // Keep track of used second endpoints
	  if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
	}
	else break;
      }
    }

    // Iterate through centers of wedges to store butterfly counts on
    for (long j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1]-v_offset;
#ifndef INVERSE
      if (v <= i) break;
#endif
      // Only store butterfly counts per center if it is on the right
      // bipartition
      if (!(GA.edges[u_offset+j] & 0b1)) continue;
      for (long k=0; k < v_deg; ++k) { 
	uintE u2_idx = GA.edges[v_offset+k] >> 1;
#ifdef INVERSE
  if (u2_idx < i && u2_idx < v) {
#else
	if (u2_idx > i) {
#endif
	  if (wedges[u2_idx] > 1) butterflies[eltsPerCacheLine*v] += wedges[u2_idx]-1;
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

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes butterfly counts per vertex, given ranked vertices.
 * 
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition to store counts on, based on which bipartition produces the fewest wedges.
 *                  If true, stores counts on U. If false, stores counts on V.
 *  num_wedges    : Number of wedges produced by the given ranking
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  max_array_size:
 *  type          : Wedge/butterfly aggregation type
 *  ranks         : Sorted U and V vertices by rank (where U indices all come before V indices)
 *  rankV         : Array mapping V indices to their rank
 *  rankU         : Array mapping U indices to their rank
 * 
 *  Returns an array of butterfly counts, indexed by vertex.
 *  Note: This function takes ownership of and frees ranks, rankV, and rankU
 */
long* CountRank(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, CountType type, 
                uintE* ranks, uintE* rankV, uintE* rankU) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

#ifdef VERBOSE
  timer t_rank;
  t_rank.start();
#endif

  // Rank graph using ranks, rankV, and rankU
  // Vertices in g are indexed by rank
  auto g = rankGraph(GA, use_v, ranks, rankV, rankU);
  free(ranks);

#ifdef VERBOSE
  t_rank.reportTotal("ranking");
  timer t_malloc;
  t_malloc.start();
#endif

  // Initialize array to hold butterfly counts per vertex, as indexed by rank
  const size_t eltsPerCacheLine = 64/sizeof(long);
  long* rank_butterflies = newA(long,eltsPerCacheLine*g.n);

#ifdef VERBOSE
  t_malloc.reportTotal("preprocess (malloc)");
  timer t_time;
  t_time.start();
#endif
  granular_for(i,0,g.n,g.n > 10000, { rank_butterflies[eltsPerCacheLine*i] = 0; });

  // Compute wedge indices so that wedges can be stored in parallel (for any
  // aggregation type except simple batching)
  long* wedge_idxs = (type == BATCHS || type == SERIAL) ? nullptr : countWedgesScan(g);

#ifdef VERBOSE
  if (type == BATCHWA) t_time.reportTotal("counting (scan)");
#endif

  // Choose wedge/butterfly aggregation type
  if (type == BATCHS) CountWorkEfficientParallel(g, rank_butterflies, max_array_size);
  else if (type == SERIAL) CountWorkEfficientSerial(g, rank_butterflies);
  else if (type == BATCHWA) CountOrigCompactParallel_WedgeAware(g, rank_butterflies, max_array_size, wedge_idxs);
  else {
    // Initialize structure that holds all space for counting algorithms, to be
    // reused with wedge batches
    CountSpace cs = CountSpace(type, g.n, true);

    if (max_wedges >= num_wedges) {
      if (type == ASORT) CountSort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == SORT) CountSortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == AHASH) CountHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == HASH) CountHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == AHIST) CountHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == HIST) CountHistCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
    }
    else {
      intT curr_idx = 0;
      while(curr_idx < g.n) {
	if (type == ASORT) curr_idx = CountSort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == SORT) curr_idx = CountSortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == AHASH) curr_idx = CountHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == HASH) curr_idx = CountHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == AHIST) curr_idx = CountHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == HIST) curr_idx = CountHistCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	cs.clear();
      }
    }
    cs.del();
  }
  g.del();
  if (type != BATCHS && type != SERIAL) free(wedge_idxs);

#ifdef VERBOSE
  if (type != BATCHS && type != BATCHWA && type != SERIAL) t_time.reportTotal("counting");
  timer t_convert;
  t_convert.start();
#endif

  // Initialize array to hold butterfly counts per vertex, by original indices
  long* butterflies = newA(long, eltsPerCacheLine*nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[eltsPerCacheLine * i] = 0; });

  // Convert all butterfly counts to be indexed by original indices
  auto rank_converter = use_v ? rankU : rankV;
  granular_for(i,0,nu,nu > 10000, {
      butterflies[eltsPerCacheLine*i] = rank_butterflies[eltsPerCacheLine*rank_converter[i]];
    });
  free(rank_butterflies);
  free(rankU); free(rankV);

#ifdef VERBOSE
  t_convert.reportTotal("convert");
#endif

  return butterflies;
}

/*
 *  Computes butterfly counts per vertex.
 * 
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition to store counts on, based on which bipartition produces the fewest wedges.
 *                  If true, stores counts on U. If false, stores counts on V.
 *  num_wedges    : Number of wedges produced by the bipartition given by use_v, assuming ranking by side
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  max_array_size:
 *  type          : Wedge/butterfly aggregation type
 *  tw            : Vertex ordering option
 * 
 *  Returns an array of butterfly counts, indexed by vertex.
 */
long* Count(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, CountType type,
            RankType tw) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

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
					   (long) 0, GA.nv, addF<long>(), rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup))
					   );
#ifdef VERBOSE
    t_rank.reportTotal("ranking");
#endif
    // Compute butterfly counts per vertex using the specified ranking
    return CountRank(GA, use_v, num_ccwedges, max_wedges, max_array_size, type, get<0>(rank_tup), get<1>(rank_tup),
                     get<2>(rank_tup));
  }

  // Ranking by side

#ifdef VERBOSE
  timer t_malloc;
  t_malloc.start();
#endif

  const size_t eltsPerCacheLine = 64/sizeof(long);
  long* butterflies = newA(long, eltsPerCacheLine*nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[eltsPerCacheLine * i] = 0; });

#ifdef VERBOSE
  t_malloc.reportTotal("preprocess (malloc)");
  timer t_time;
  t_time.start();
#endif

  // Compute wedge indices so that wedges can be stored in parallel (for any
  // aggregation type except simple batching)
  long* wedge_idxs = (type == BATCHS || type == SERIAL) ? nullptr : countWedgesScan(GA, use_v, true);

#ifdef VERBOSE
  if (type == BATCHWA) t_time.reportTotal("counting (scan)");
#endif

  // Choose wedge/butterfly aggregation type
  if (type == BATCHS) CountOrigCompactParallel(GA, use_v, butterflies, max_array_size);
  else if (type == BATCHWA) CountOrigCompactParallel_WedgeAware(GA, use_v, butterflies, max_array_size, wedge_idxs);
  else if (type == SERIAL) CountOrigCompactSerial(GA, use_v, butterflies);
  else {
    // Initialize structure that holds all space for counting algorithms, to be
    // reused with wedge batches
    CountSpace cs = CountSpace(type, nu, false);

    if (max_wedges >= num_wedges) {
      if (type == ASORT) CountSort(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
      else if (type == SORT) CountSortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
      else if (type == AHASH) CountHash(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
      else if (type == HASH) CountHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
      else if (type == AHIST) CountHist(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
      else if (type == HIST) CountHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    }
    else {
      intT curr_idx = 0;
      while(curr_idx < nu) {
	if (type == ASORT) curr_idx = CountSort(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == SORT) curr_idx = CountSortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == AHASH) curr_idx = CountHash(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == HASH) curr_idx = CountHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == AHIST) curr_idx = CountHist(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == HIST) curr_idx = CountHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
	cs.clear();
      }
    }
    cs.del();
  }

  if (type != BATCHS && type != SERIAL) free(wedge_idxs);

#ifdef VERBOSE
  if (type != BATCHS && type != BATCHWA && type != SERIAL) t_time.reportTotal("counting");
#endif

  return butterflies;
}

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#include "hygra.h"
#include "index_map.h"
#include "bucket.h"
#include "edgeMapReduce.h"

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "hypergraphIO.h"
#include "parseCommandLine.h"
#include "vertex.h"
#include "sequence.h"
#include "sparseSet.h"
#include "sampleSort.h"
#ifndef OPENMP
#include "../../lib/histogram.h"
#include "../../lib/sample_sort.h"
#endif
#include "../../radixsort/RadixSort/radixSort.h"

#include "butterfly_putils.h"
#include "butterfly_utils.h"

/*
 *  Peels butterflies given bucket of active edges, using histogram aggregation.
 * 
 *  eti         : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite         : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current     : Array to keep track of active edges
 *  ps          : Holds all array space needed, to be reused between buckets
 *  active      : Set of active edges
 *  butterflies : Butterfly counts per edge
 *  update_dense: Keeps track of edges with updated butterfly counts
 *  GA          : Bipartite graph in CSR format
 *  use_v       : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges  : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                memory
 *  curr_idx    : Denotes the starting vertex in this batch
 * 
 *  Returns the next vertex to process for the next batch, and the number of
 *  edges with updated butterfly counts.
 */
pair<intT, long> PeelESort(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active,
                           long* butterflies, bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx=0) {
  using X = tuple<uintE,long>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve number of butterflies peeled on corresponding edges, after peeling active subset
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto ret = getIntersectWedges(eti, ite, current, ps, active_map, active.size(), GA, use_v, max_wedges, curr_idx);

  // Aggregate butterfly count updates
#ifdef OPENMP
  sampleSort(ps.wedges_seq_tup_fil.A, ret.first, tupleLt<uintE,long>());
#else
  pbbs::sample_sort(ps.wedges_seq_tup_fil.A, ret.first, tupleLt<uintE,long>());
#endif
  auto b_freq_pair = getFreqs<long>(ps.wedges_seq_tup_fil.A, ret.first, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;
  ps.resize_update(b_freq_pair.second-1);
  uintE* update = ps.update_seq_int.A;

  // Update butterfly counts for those edges
  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    long num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(ps.wedges_seq_tup_fil.A[b_freq_arr[i-1]]), num_freq, tupleAdd<uintE, long>());
    // These are our butterfly counts
    butterflies[eltsPerCacheLine*get<0>(reduce)] -= get<1>(reduce);
    update[i-1] = get<0>(reduce);
  }

  free(b_freq_arr);

  return make_pair(ret.second, b_freq_pair.second-1);
}

/*
 *  Peels butterflies given bucket of active edges, using histogram aggregation.
 * 
 *  eti         : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite         : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current     : Array to keep track of active edges
 *  ps          : Holds all array space needed, to be reused between buckets
 *  active      : Set of active edges
 *  butterflies : Butterfly counts per edge
 *  update_dense: Keeps track of edges with updated butterfly counts
 *  GA          : Bipartite graph in CSR format
 *  use_v       : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges  : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                memory
 *  curr_idx    : Denotes the starting vertex in this batch
 * 
 *  Returns the next vertex to process for the next batch, and the number of
 *  edges with updated butterfly counts.
 */
pair<intT, long> PeelEHist(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active,
                           long* butterflies, bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx=0) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<uintE,long>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Retrieve number of butterflies peeled on corresponding edges, after peeling active subset
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto ret = getIntersectWedges(eti, ite, current, ps, active_map, active.size(), GA, use_v, max_wedges, curr_idx);
  pbbsa::sequence<X> wedge_seq = pbbsa::sequence<X>(ps.wedges_seq_tup_fil.A, ret.first);

  // Aggregate butterfly count updates
  pbbsa::sequence<tuple<uintE, long>> tmp = pbbsa::sequence<tuple<uintE, long>>();
  pbbsa::sequence<tuple<uintE, long>> out = pbbsa::sequence<tuple<uintE, long>>();
  tuple<size_t, X*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE, long>(wedge_seq, GA.numEdges, getAdd<uintE,long>, getAddReduce<uintE,long>, tmp,
                                           out);
  X* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  ps.resize_update(butterflies_n);
  uintE* update = ps.update_seq_int.A;

  // Update butterfly counts for those edges
  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] -= get<1>(butterflies_l[i]);
    update[i] = get<0>(butterflies_l[i]);
  }

  return make_pair(ret.second, butterflies_n);
#endif
}

/*
 *  Peels butterflies given bucket of active edges, using hash aggregation.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite           : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current       : Array to keep track of active edges
 *  ps            : Holds all array space needed, to be reused between buckets
 *  active        : Set of active edges
 *  butterflies   : Butterfly counts per edge
 *  update_dense  : Keeps track of edges with updated butterfly counts
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 * 
 *  Returns the number of edges with updated butterfly counts.
 */
long PeelEHash(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, long* butterflies,
               bipartiteCSR& GA, bool use_v) {
  // Retrieve number of butterflies peeled on corresponding edges, after peeling active subset
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  getIntersectWedgesHash(eti, ite, current, ps, active_map, active.size(), GA, use_v);
 
  auto update_seq = ps.update_hash.entries();
  long num_updates = update_seq.n;
  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;

  const size_t eltsPerCacheLine = 64/sizeof(long);
  // Update butterfly counts for those edges
  granular_for(i, 0, num_updates, (num_updates>1000), {
      auto update_pair = update_seq.A[i];
      uintE idx = update_pair.first;
      butterflies[eltsPerCacheLine*idx] -= update_pair.second;
      update[i] = idx;
    });

  update_seq.del();

  return num_updates;
}

/*
 *  Peels butterflies given bucket of active edges, using wedge-aware batching.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite           : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current       : Array to keep track of active edges
 *  ps            : Holds all array space needed, to be reused between buckets
 *  active        : Set of active edges
 *  butterflies   : Butterfly counts per edge
 *  update_dense  : Keeps track of edges with updated butterfly counts
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 * 
 *  Returns an array of edges with updated butterfly counts and the size of that
 *  array.
 */
pair<uintE*, long> PeelEOrigParallel_WedgeAware(uintE* eti, uintE* ite, bool* current, PeelESpace& ps,
                                                vertexSubset& active, long* butterflies, bool* update_dense,
                                                bipartiteCSR& GA, bool use_v) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = active.size() < MAX_STEP_SIZE ? active.size() : MAX_STEP_SIZE; // tunable parameter

  
  auto wedges_seq = ps.wedges_seq_int;
  auto used_seq = ps.used_seq_int;

  // Stores number of wedges per endpoints
  auto wedges = wedges_seq.A;
  // Keeps track of entries used in the wedges array
  auto used = used_seq.A;

  const size_t eltsPerCacheLine = 64/sizeof(long);
  // Bookkeeping
  granular_for(i,0,GA.numEdges,GA.numEdges > 10000, { update_dense[eltsPerCacheLine*i] = false; });
  parallel_for(intT i=0; i < active.size(); ++i){ current[active.vtx(i)] = 1; }

  // Determine wedge indices so that they can be processed in proper-sized batches
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  long* wedgesPrefixSum = getIntersectWedgeIdxs(ite, active_map, active.size(), GA, use_v);

  // Consider active vertices i in batches, as given by wedgesPrefixSum
  for(intT step = 0; step < (active.size()+stepSize-1)/stepSize; step++) {
    std::function<void(intT,intT)> recursive_lambda =
      [&]
      (intT start, intT end){
      if ((start == end-1) || (wedgesPrefixSum[end]-wedgesPrefixSum[start] < 1000)){ 
	for (intT i = start; i < end; i++){
	  intT used_idx = 0;
	  intT shift = nu*(i-step*stepSize);
	  intT idx_vu = active.vtx(i);

	  // Retrieve active edge (u, v)
	  uintE u = edgesV[idx_vu];
	  intT u_offset = offsetsU[u];
	  intT u_deg = offsetsU[u+1] - u_offset;

	  uintE v = edgesU[ite[idx_vu]];
	  intT v_offset = offsetsV[v];
	  intT v_deg = offsetsV[v+1] - v_offset;

	  // Iterate through all neighbors of v
	  for (intT k=0; k < v_deg; ++k) {
	    uintE u2 = edgesV[v_offset + k];
	    // Check that the neighbor u2 hasn't been peeled before
	    if (u2 != u && u2 != UINT_E_MAX && (!current[v_offset+k] || idx_vu < v_offset + k)) {
	      intT u2_offset = offsetsU[u2];
	      intT u2_deg = offsetsU[u2+1] - u2_offset;
	      intT int_offset = 0;

	      // Intersect N(u) with N(u2)
	      // Iterate through all neighbors of u2
	      for(intT j = 0; j < u2_deg; ++j) {
		uintE v2 = edgesU[u2_offset + j];
		uintE idx_v2u2 = eti[u2_offset + j];
		// Check that the neighbor v2 hasn't been peeled before
		if(v2 != UINT_E_MAX && v2 != v && (!current[idx_v2u2] || idx_v2u2 >= idx_vu)) {
		  // Check if v2 is also a a neighbor of u
		  while(int_offset < u_deg &&
			(edgesU[u_offset + int_offset] == UINT_E_MAX || edgesU[u_offset + int_offset] < v2)) { int_offset++; }
		  if (int_offset < u_deg && edgesU[u_offset+int_offset] == v2 &&
		      (!current[eti[u_offset + int_offset]] || idx_vu < eti[u_offset + int_offset])) {
		    // If v2 is a neighbor of u and u2, then increment our count
		    wedges[shift+k]++;
		    // Subtract the butterfly from (v2, u2) and (v2, u)
		    writeAdd(&butterflies[eltsPerCacheLine*eti[u2_offset + j]], (long) -1);
		    writeAdd(&butterflies[eltsPerCacheLine*eti[u_offset + int_offset]], (long) -1);
		    // Ensure that  (v2, u2), (v2, u), and (v, u2) have been marked with changed butterfly counts
		    if(!update_dense[eltsPerCacheLine*eti[u2_offset + j]])
		      CAS(&update_dense[eltsPerCacheLine*eti[u2_offset + j]], false, true);
		    if(!update_dense[eltsPerCacheLine*eti[u_offset + int_offset]])
		      CAS(&update_dense[eltsPerCacheLine*eti[u_offset + int_offset]], false, true);
		    if (wedges[shift+k] == 1) {
		      used[shift+used_idx++] = k;
		      if(!update_dense[eltsPerCacheLine*(v_offset+k)]) CAS(&update_dense[eltsPerCacheLine*(v_offset+k)],false,true);
		    }
		  }
		  else if(int_offset >= u_deg) break;
		}
	      }
	    }
	  }

	  // Update the butterfly count on (v, u2) and clear the wedges array
	  granular_for(j,0,used_idx,used_idx > 10000, { 
	      intT k = used[shift+j];
	      writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)],(long) -1* ((long)wedges[shift+k]));
	      wedges[shift+k] = 0;
	    });
	}
      } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,active.size()));
  }
  free(wedgesPrefixSum);

  // Update array that keeps track of active edges
  parallel_for(intT i=0; i < active.size(); ++i){
    edgesV[active.vtx(i)] = UINT_E_MAX; edgesU[ite[active.vtx(i)]] = UINT_E_MAX;
    current[active.vtx(i)] = 0;
  }

  // Collate edges with updated butterfly counts
  auto f = [&] (size_t i) { return update_dense[eltsPerCacheLine*i]; };
  auto f_in = make_in_imap<bool>(GA.numEdges, f);
  auto out = pbbs::pack_index<uintE>(f_in);
  out.alloc = false;

  return make_pair(out.s, out.size());
}

/*
 *  Peels butterflies given bucket of active edges, using simple batching.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite           : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current       : Array to keep track of active edges
 *  ps            : Holds all array space needed, to be reused between buckets
 *  active        : Set of active edges
 *  butterflies   : Butterfly counts per edge
 *  update_dense  : Keeps track of edges with updated butterfly counts
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 * 
 *  Returns an array of edges with updated butterfly counts and the size of that
 *  array.
 */
pair<uintE*, long> PeelEOrigParallel(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active,
                                     long* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = active.size() < MAX_STEP_SIZE ? active.size() : MAX_STEP_SIZE; // tunable parameter

  auto wedges_seq = ps.wedges_seq_int;
  auto used_seq = ps.used_seq_int;

  // Stores number of wedges per endpoints
  uintE* wedges = wedges_seq.A;
  // Keeps track of entries used in the wedges array
  uintE* used = used_seq.A;

  const size_t eltsPerCacheLine = 64/sizeof(long);
  // Bookkeeping
  granular_for(i,0,GA.numEdges,GA.numEdges > 10000, { update_dense[eltsPerCacheLine*i] = false; });
  parallel_for(intT i=0; i < active.size(); ++i){ current[active.vtx(i)] = 1; }

  // Consider active vertices i in batches
  for(intT step = 0; step < (active.size()+stepSize-1)/stepSize; step++) {
    parallel_for_1(intT i=step*stepSize; i < min((step+1)*stepSize,active.size()); ++i){
      intT used_idx = 0;
      intT shift = nu*(i-step*stepSize);
      intT idx_vu = active.vtx(i);
  
      // Retrieve active edge (u, v)
      uintE u = edgesV[idx_vu];
      intT u_offset = offsetsU[u];
      intT u_deg = offsetsU[u+1] - u_offset;

      uintE v = edgesU[ite[idx_vu]];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
  
      // Iterate through all neighbors of v
      for (intT k=0; k < v_deg; ++k) { 
	uintE u2 = edgesV[v_offset + k];
	// Check that the neighbor u2 hasn't been peeled before
	if (u2 != u && u2 != UINT_E_MAX && (!current[v_offset+k] || idx_vu < v_offset + k)) {
	  intT u2_offset = offsetsU[u2];
	  intT u2_deg = offsetsU[u2+1] - u2_offset;
	  intT int_offset = 0;
	  // Intersect N(u) with N(u2)
	  // Iterate through all neighbors of u2
	  for(intT j = 0; j < u2_deg; ++j) {
	    uintE v2 = edgesU[u2_offset + j];
	    uintE idx_v2u2 = eti[u2_offset + j];
	    // Check that the neighbor v2 hasn't been peeled before
	    if(v2 != UINT_E_MAX && v2 != v && (!current[idx_v2u2] || idx_v2u2 >= idx_vu)) {
	      // Check if v2 is also a a neighbor of u
	      while(int_offset < u_deg &&
		    (edgesU[u_offset + int_offset] == UINT_E_MAX || edgesU[u_offset + int_offset] < v2)) {
		int_offset++;
	      }
	      if (int_offset < u_deg && edgesU[u_offset+int_offset] == v2 &&
		  (!current[eti[u_offset + int_offset]] || idx_vu < eti[u_offset + int_offset])) {
		// If v2 is a neighbor of u and u2, then increment our count
		wedges[shift+k]++;
		// Subtract the butterfly from (v2, u2) and (v2, u)
		writeAdd(&butterflies[eltsPerCacheLine*eti[u2_offset + j]], (long) -1);
		writeAdd(&butterflies[eltsPerCacheLine*eti[u_offset + int_offset]], (long) -1);
		// Ensure that  (v2, u2), (v2, u), and (v, u2) have been marked with changed butterfly counts
		if(!update_dense[eltsPerCacheLine*eti[u2_offset + j]])
		  CAS(&update_dense[eltsPerCacheLine*eti[u2_offset + j]], false, true);
		if(!update_dense[eltsPerCacheLine*eti[u_offset + int_offset]])
		  CAS(&update_dense[eltsPerCacheLine*eti[u_offset + int_offset]], false, true);
		if (wedges[shift+k] == 1) {
		  used[shift+used_idx++] = k;
		  if(!update_dense[eltsPerCacheLine*(v_offset+k)]) CAS(&update_dense[eltsPerCacheLine*(v_offset+k)],false,true);
		}
	      }
	      else if(int_offset >= u_deg) break;
	    }
	  }
	}
      }

      // Update the butterfly count on (v, u2) and clear the wedges array
      granular_for(j,0,used_idx,used_idx > 10000, { 
	  intT k = used[shift+j];
	  writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], (long) -1 * ((long)wedges[shift+k]));
	  wedges[shift+k] = 0;
	});
    }
  }

  // Update array that keeps track of active edges
  parallel_for(intT i=0; i < active.size(); ++i){
    edgesV[active.vtx(i)] = UINT_E_MAX; edgesU[ite[active.vtx(i)]] = UINT_E_MAX;
    current[active.vtx(i)] = 0;
  }

  // Collate edges with updated butterfly counts
  auto f = [&] (size_t i) { return update_dense[eltsPerCacheLine*i]; };
  auto f_in = make_in_imap<bool>(GA.numEdges, f);
  auto out = pbbs::pack_index<uintE>(f_in);
  out.alloc = false;
  return make_pair(out.s, out.size());
}

/*
 *  Peels butterflies given bucket of active edges, and updates the bucketing
 *  structure appropriately.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite           : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current       : Array to keep track of active edges
 *  ps            : Holds all array space needed, to be reused between buckets
 *  active        : Set of active edges
 *  butterflies   : Butterfly counts per edge
 *  update_dense  : Depending on aggregation type, keeps track of edges with updated butterfly counts
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  type          : Wedge aggregation type
 *  D             : Map of edges to the buckets that they're in
 *  b             : Bucketing structure
 *  k             : ID of current bucket
 */
void PeelE_helper(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, long* butterflies,
                  bool* update_dense, bipartiteCSR& GA, bool use_v, long max_wedges, PeelType type, array_imap<long>& D,
                  buckets<array_imap<long>>& b, long k) {
  const long nu = use_v ? GA.nu : GA.nv;
  bool is_seq = (active.size() < 1000);

  // Choose wedge aggregation type
  if (type == PBATCHS) {
    auto ret = PeelEOrigParallel(eti, ite, current, ps, active, butterflies, update_dense, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ret.first, ret.second);
    free(ret.first);
    return;
  }
  else if (type == PBATCHWA) {
    auto ret = PeelEOrigParallel_WedgeAware(eti, ite, current, ps, active, butterflies, update_dense, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ret.first, ret.second);
    free(ret.first);
    return;
  }
  else if (type == PHASH) {
    auto num_updates = PeelEHash(eti, ite, current, ps, active, butterflies, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ps.update_seq_int.A, num_updates);
    ps.clear();
    return;
  }

  pair<intT, long> ret;
  intT curr_idx = 0;
  while(curr_idx < active.size()) {
    if (type == PSORT) ret = PeelESort(eti, ite, current, ps, active, butterflies, GA, use_v, max_wedges, curr_idx);
    else if (type == PHIST)
      ret = PeelEHist(eti, ite, current, ps, active, butterflies, GA, use_v, max_wedges, curr_idx);
    curr_idx = ret.first;
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ps.update_seq_int.A, ret.second);
  }
}

/*
 *  Peels butterflies per edge.
 * 
 *  eti           : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite           : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  butterflies   : Butterfly counts per edge
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  type          : Wedge aggregation type
 *  num_buckets   : Number of buckets to initialize bucketing structure with
 * 
 *  Returns a map of edges to the buckets they were peeled in.
 */
array_imap<long> PeelE(uintE* eti, uintE* ite, bipartiteCSR& GA, bool use_v, long* butterflies, long max_wedges,
                       PeelType type, size_t num_buckets=128) {
  // Butterflies are stored on GA.U if use_v is true, and on GA.V otherwise
  const long nu = use_v ? GA.nu : GA.nv;

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  
  using X = tuple<uintE,uintE>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Construct initial buckets based on butterfly counts
  auto D = array_imap<long>(GA.numEdges, [&] (size_t i) { return butterflies[eltsPerCacheLine*i]; });
  auto b = make_buckets(GA.numEdges, D, increasing, num_buckets);

  // Set up space for peeling algorithms
  bool* update_dense = nullptr;
  if (type == PBATCHS || type == PBATCHWA) update_dense = newA(bool, eltsPerCacheLine*GA.numEdges);
  PeelESpace ps = PeelESpace(type, GA.numEdges, MAX_STEP_SIZE, nu);

  // Set up array to keep track of current active vertices
  bool* current = newA(bool, GA.numEdges);
  parallel_for(size_t i=0; i < GA.numEdges; ++i) { current[i] = false; }

  size_t finished = 0;

  // Peel each bucket
  while (finished != GA.numEdges) {
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    if (active.size() == 0) {active.del(); continue;}
    long k = bkt.id;
    finished += active.size();
    bool is_seq = (active.size() < 1000);

    // Peel butterflies
    PeelE_helper(eti, ite, current, ps, active, butterflies, update_dense, GA, use_v, max_wedges, type, D, b, k);

    active.del();
  }

  ps.del();
  b.del();
  if (type == PBATCHS || type == PBATCHWA) free(update_dense);
  free(current);
  
  return D;
}

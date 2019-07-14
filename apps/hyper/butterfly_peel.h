#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#include "hypergraphIO.h"
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
#include "../../lib/histogram.h"
//#include "../../lib/gbbs-histogram.h"
#include "../../lib/sample_sort.h"
#include "../../radixsort/RadixSort/radixSort.h"

#include "butterfly_putils.h"
#include "butterfly_utils.h"

template <class E>
struct seagullSumHelper { 
  uintE u;
  uintT* offsetsU;
  uintT* offsetsV;
  uintE* edgesU;
seagullSumHelper(uintE _u, uintT* _offsetsU, uintT* _offsetsV, uintE* _edgesU) : u(_u), offsetsU(_offsetsU), offsetsV(_offsetsV), edgesU(_edgesU) {}
  inline E operator() (const E& i) const {
    intT u_offset = offsetsU[u];
    uintE v = edgesU[u_offset + i];
    return (E) (offsetsV[v+1] - offsetsV[v] - 1);
  }
};

template <class E>
struct seagullSum { 
  uintT* offsetsU;
  uintT* offsetsV;
  uintE* edgesU;
  uintE* active;
seagullSum(uintT* _offsetsU, uintT* _offsetsV, uintE* _edgesU, uintE* _active) : offsetsU(_offsetsU), offsetsV(_offsetsV), edgesU(_edgesU), active(_active) {}
  inline E operator() (const E& i) const {
    /*uintE u = active[i];
      intT u_offset = offsetsU[u];
      intT u_deg = offsetsU[active[i]+1] - offsetsU[active[i]];
      E ret=0;
      for (long k=0; k < u_deg; ++k) {
      uintE v = edgesU[u_offset + k];
      ret += (offsetsV[v+1] - offsetsV[v] - 1);
      }
      return ret;*/
    E u_deg = offsetsU[active[i]+1] - offsetsU[active[i]];
    return sequence::reduce<E>((E) 0, u_deg, addF<E>(), seagullSumHelper<E>(active[i], offsetsU, offsetsV, edgesU));
  }
};

/*
 *  Counts the number of seagulls on vertices in active. The vertices in active are assumed to
 *  be a subset of U, and V represents the other bipartition of the bipartite graph.
 * 
 *  V     : One bipartition of vertices
 *  U     : The other bipartition of vertices
 *  active: Set of active vertices (subset of U)
 * 
 *  Returns: Number of seagulls on vertices in active
 */
long countSeagulls_seq(bipartiteCSR& GA, bool use_v, vertexSubset active) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long num_sg = 0;
  for (uintT i=0; i < active.size(); ++i) {
    uintE u = active.vtx(i);
    uintT u_offset = offsetsU[u];
    long u_deg = offsetsU[u+1] - u_offset;
    for (uintT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      num_sg += (offsetsV[v+1] - offsetsV[v] - 1);
    }
  }
  return num_sg;
}

long countSeagulls(bipartiteCSR& GA, bool use_v, vertexSubset active) {
  if (active.size() < 1000) return countSeagulls_seq(GA, use_v, active);

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  return sequence::reduce<long>((long) 0, (long) active.size(), addF<long>(), seagullSum<long>(offsetsU, offsetsV, edgesU, active.s));
}


//***************************************************************************************************
//***************************************************************************************************

/*
 *  Retrieves the number of butterflies remaining on each non-active vertex (where the number of 
 *  butterflies has changed), given an array of seagulls (constructed from some set of active vertices).
 *  Uses sorting to obtain these counts.
 * 
 *  nu         : The number of vertices in U, which represents one bipartition of vertices, of which 
 *               the active vertices are taken to be a subset of
 *  seagulls   : Array of seagulls constructed from some set of active vertices. The active vertex that
 *               forms the endpoint of a seagull is assumed to be the first vertex in the unordered
 *               vertex pair, and the non-active vertex is assumed to be the second vertex.
 *  num_sgs    : The number of seagulls
 *  butterflies: An updated array containing the number of butterflies on each vertex in U (where
 *               the vertex is given by the index in butterflies)
 * 
 *  Returns: Array of butterfly count changes, including both vertex indices and the changed values
 */
pair<tuple<uintE,long>*, long> getSeagullFreqs(const long nu, UVertexPair* seagulls, long num_sgs) {
  using X = tuple<uintE,long>;
  // Sort seagulls (considering both active + non-active endpoints), and retrieve frequency counts
  //radix::parallelIntegerSort<uintE>(seagulls, num_sgs, UVPFirst());
  //radix::parallelIntegerSort<uintE>(seagulls, num_sgs, UVPSecond());
  pbbs::sample_sort(seagulls, num_sgs, UVertexPairCmp2());
  pair<long*, long> freq_pair = getFreqs<long>(seagulls, num_sgs, UVertexPairCmp2(), UVertexPairEq(), LONG_MAX, nonMaxLongF());
  long num_sg_freqs = freq_pair.second - 1;
  X* sg_freqs = newA(X, num_sg_freqs);
  // When retrieving frequency counts, store the frequency choose 2 with the non-active endpoint
  // This gives us the number of butterflies to be removed on the non-active endpoint
  parallel_for(long i=1; i < freq_pair.second; ++i) {
    long num = freq_pair.first[i] - freq_pair.first[i-1];
    uintE idx = seagulls[freq_pair.first[i-1]].v2;
    sg_freqs[i-1] = make_tuple(idx, (num * (num-1)) / 2);
  }
  free(freq_pair.first);

  // Filter out any entries that have 0 butterflies to be removed (bucketing cannot handle these)
  //X* sg_freqs_filter = newA(X,num_sg_freqs);
  //num_sg_freqs = sequence::filter(sg_freqs,sg_freqs_filter,num_sg_freqs,nonZeroF());

  //free(sg_freqs);
  //return make_pair(sg_freqs_filter,num_sg_freqs);
  return make_pair(sg_freqs,num_sg_freqs);
}

pair<tuple<uintE,long>*, long> getSeagullFreqs_seq(const long nu, UVertexPair* seagulls, long num_sgs) {
  //radix::parallelIntegerSort<uintE>(seagulls, num_sgs, UVPFirst());
  //radix::parallelIntegerSort<uintE>(seagulls, num_sgs, UVPSecond());
  pbbs::sample_sort(seagulls, num_sgs, UVertexPairCmp2());
  return getFreqs_seq<long,uintE>(seagulls, num_sgs, UVertexPairCmp2(), UVertexPairEq(), true,
				  UVertexPairV2(), choose2(), reflCount<UVertexPair>());
}

long getUpdates_seq(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, long* butterflies) {
  using X = tuple<uintE,long>;
  //radix::parallelIntegerSort<uintE>(sg_freqs, num_sg_freqs, tupleFirst<uintE,long>());
  pbbs::sample_sort(sg_freqs, num_sg_freqs, tupleLt<uintE,long>());
  pair<X*, long> b_freq_pair = getFreqs_seq<long,uintE>(sg_freqs, num_sg_freqs, tupleLt<uintE,long>(), tupleEq<uintE,long>(), false,
							uintETupleGet0(), refl<long>(), uintECount());
  long num_updates = b_freq_pair.second;
  //uintE* update = newA(uintE, num_updates);
  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;
  X* b_updates = b_freq_pair.first;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  for(long i=0; i<num_updates; ++i) {
    uintE u_idx = get<0>(b_updates[i]);
    butterflies[eltsPerCacheLine*u_idx] -= get<1>(b_updates[i]);
    //update_dense[eltsPerCacheLine*u_idx] = get<1>(b_updates[i]) > 0 ? true : false; //this is b/c seq doesn't filter
    update[i] = u_idx;
  }
  free(b_updates);
  return num_updates;
}


/*
 *  Collates butterfly count updates on each non-active vertex, and updates butterflies with
 *  these values.
 *  Uses sorting to collate these updated counts.
 * 
 *  sg_freqs    : Array of butterfly counts to be removed from corresponding vertex indices.
 *                Vertex indices are not necessarily distinct amongst entries.
 *  num_sg_freqs: Length of sg_freqs
 *  butterflies : An updated array containing the number of butterflies on each vertex in U (where
 *                the vertex is given by the index in butterflies)
 * 
 *  Returns: Array of vertex indices on which the butterfly counts change
 */
long getUpdates(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, long* butterflies) {
  using X = tuple<uintE,long>;
  // Now, collate all butterflies to be removed with the same non-active endpoint
  // Do this by sorting on the non-active endpoint, and summing the frequencies
  //radix::parallelIntegerSort<uintE>(sg_freqs, num_sg_freqs, tupleFirst<uintE,long>());
  pbbs::sample_sort(sg_freqs, num_sg_freqs, tupleLt<uintE,long>());
  auto b_freq_pair = getFreqs<long>(sg_freqs, num_sg_freqs, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF(), false);
  long num_updates = b_freq_pair.second - 1;
  //uintE* update = newA(uintE, num_updates);
  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  //uintE* update = newA(uintE, num_updates);
  parallel_for(long i=1; i < num_updates + 1; ++i) {
    long num_freq = b_freq_pair.first[i] - b_freq_pair.first[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(sg_freqs[b_freq_pair.first[i-1]]), num_freq, tupleAdd<uintE,long>());
    uintE u_idx = get<0>(reduce);
    // Remove these butterflies from our array of butterfly counts
    butterflies[eltsPerCacheLine*u_idx] -= get<1>(reduce);
    // Add the updated index to an update array, to be returned
    //update_dense[eltsPerCacheLine*u_idx] = get<1>(reduce) > 0 ? true : false;
    update[i-1] = u_idx;
  }
  free(b_freq_pair.first);

  return num_updates;
}

/*
 *  Precisely getSeagullFreqs, but using histograms instead of sorting.
 */
pair<tuple<uintE,long>*, long> getSeagullFreqsHist(PeelSpace& ps, const long nu, long* seagulls, long num_sgs) {
  using X = tuple<long,uintE>;
  using T = tuple<uintE,long>;
  // TODO integrate sequence into histogram code (so we don't have to convert?)
  pbbsa::sequence<long> sg_seq = pbbsa::sequence<long>(seagulls,num_sgs);

  // Place seagulls into a histogram to retrieve frequency counts (considering both active + non-active endpoints)
  tuple<size_t,X*> sg_hist_tuple = pbbsa::sparse_histogram<long,uintE>(sg_seq,nu*nu+nu,ps.tmp,ps.out);
  X* sg_freqs = get<1>(sg_hist_tuple);

  // Filter out any frequency count <= 1, since these won't contribute towards a butterfly
  X* sg_freqs_filter = newA(X,get<0>(sg_hist_tuple));
  size_t num_sg_freqs = sequence::filter(sg_freqs, sg_freqs_filter, get<0>(sg_hist_tuple), greaterOneLongF());
  //free(sg_freqs);
  T* sg_freqs_ret = newA(T, num_sg_freqs);

  // Update our frequency counts to represent the number of butterflies to remove on each
  // non-active endpoint
  // Do this by storing the frequency choose 2 on the non-active endpoint
  parallel_for(long i=0; i < num_sg_freqs; ++i) {
    long num = get<1>(sg_freqs_filter[i]);
    // The non-active endpoint is given by our seagull id % nu, because it is always the second endpoint
    sg_freqs_ret[i] = make_tuple(get<0>(sg_freqs_filter[i]) % nu, num * (num-1) / 2);
  }
  free(sg_freqs_filter);
  return make_pair(sg_freqs_ret, num_sg_freqs);
}

pair<tuple<uintE,long>*, long> getSeagullFreqsHist_seq(PeelSpace& ps, const long nu, long* seagulls, long num_sgs) {
  using X = tuple<long,uintE>;
  using T = tuple<uintE,long>;
  pbbsa::sequence<long> sg_seq = pbbsa::sequence<long>(seagulls,num_sgs);
  tuple<size_t,X*> sg_hist_tuple = pbbsa::sparse_histogram<long,uintE>(sg_seq,nu*nu+nu,ps.tmp,ps.out);
  X* sg_freqs = get<1>(sg_hist_tuple);
  size_t num_sg_freqs = get<0>(sg_hist_tuple);
  T* sg_freqs_ret = newA(T, num_sg_freqs);
  long idx = 0;
  for(long i=0; i < num_sg_freqs; ++i) {
    long num = get<1>(sg_freqs[i]);
    if (num > 1) {
      sg_freqs_ret[idx] = make_tuple(get<0>(sg_freqs[i]) % nu, num * (num-1)/2);
      idx++;
    }
  }
  //free(sg_freqs);
  return make_pair(sg_freqs_ret, idx);
}

// TODO just use granular for on getUpdatesHis?
long getUpdatesHist_seq(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, const long nu, long* butterflies) {
  using X = tuple<uintE, long>;
  using T = tuple<long, uintE>;
  pbbsa::sequence<X> sg_freqs_seq = pbbsa::sequence<X>(sg_freqs,num_sg_freqs);
  tuple<size_t,X*> b_hist_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(sg_freqs_seq,nu,getAdd<uintE,long>,getAddReduce<uintE,long>);
  X* b_freqs = get<1>(b_hist_tuple);
  const size_t eltsPerCacheLine = 64/sizeof(long);
  size_t num_updates = get<0>(b_hist_tuple);
  //uintE* update = newA(uintE, num_updates);
  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;
  for(long i=0; i < num_updates; ++i) {
    uintE u_idx = get<0>(b_freqs[i]);
    butterflies[eltsPerCacheLine*u_idx] -=get<1>(b_freqs[i]);
    //update_dense[eltsPerCacheLine*u_idx] = true;
    update[i] = u_idx;
  }
  free(b_freqs);
  return (long) num_updates;
}

/*
 *  Precisely getUpdates, but using histograms instead of sorting.
 */
long getUpdatesHist(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, const long nu, long* butterflies) {
  using X = tuple<uintE,long>;
  // Now, collate all butterflies to be removed with the same non-active endpoint
  // Do this by using a histogram to sum frequencies on the non-active endpoint
  pbbsa::sequence<X> sg_freqs_seq = pbbsa::sequence<X>(sg_freqs,num_sg_freqs);
  tuple<size_t,X*> b_hist_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(sg_freqs_seq,nu,getAdd<uintE,long>,getAddReduce<uintE,long>);
  X* b_freqs = get<1>(b_hist_tuple);
  size_t num_updates = get<0>(b_hist_tuple);
  //uintE* update = newA(uintE, num_updates);
  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Remove these butterflies from our array of butterfly counts, and change our frequency
  // array to be an update array on the new number of butterflies on our non-active endpoint
  parallel_for(long i=0; i < num_updates; ++i) {
    uintE u_idx = get<0>(b_freqs[i]);
    butterflies[eltsPerCacheLine*u_idx] -=get<1>(b_freqs[i]);
    //update_dense[eltsPerCacheLine*u_idx] = true;
    update[i] = u_idx;
  }
  free(b_freqs);

  return (long) num_updates;
}

//***************************************************************************************************
//***************************************************************************************************

/*
 *  Computes updated butterfly counts after active vertices are removed. The vertices in active are
 *  assumed to be a subset of U, and V represents the other bipartition of the bipartite graph.
 *  Uses sorting to obtain counts.
 * 
 *  active     : Set of active vertices (subset of U)
 *  butterflies: An updated array containing the number of butterflies on each vertex in U (where
 *               the vertex is given by the index in butterflies)
 *  V          : One bipartition of vertices
 *  U          : The other bipartition of vertices
 *  nu         : The number of vertices in U
 * 
 *  Returns: Array of butterfly counts to be changed with corresponding vertex, and the number of
 *           such changes
 */
pair<intT, long> PeelSort(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long,intT> sg_pair = getActiveWedges<UVertexPair>(ps.wedges_seq_uvp, active_map, active.size(), GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  auto sg_freqs_pair = getSeagullFreqs(nu, ps.wedges_seq_uvp.A, sg_pair.first);

  // Compute updated butterfly counts
  auto ret = getUpdates(ps, sg_freqs_pair.first, sg_freqs_pair.second, butterflies);

  free(sg_freqs_pair.first);
  return make_pair(next_idx, ret);
}

//TODO these are all CE versions -- do version where we don't use getUpdates_seq, and we writeadd (ish) to update
pair<intT, long> PeelSort_seq(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long, intT> sg_pair = getActiveWedges<UVertexPair>(ps.wedges_seq_uvp, active_map, active.size(), GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  auto sg_freqs_pair = getSeagullFreqs_seq(nu, ps.wedges_seq_uvp.A, sg_pair.first);

  // Compute updated butterfly counts
  auto ret = getUpdates_seq(ps, sg_freqs_pair.first, sg_freqs_pair.second, butterflies);

  free(sg_freqs_pair.first);
  return make_pair(next_idx, ret);
}

/*
 *  Precisely PeelSort, but using histograms instead of repeated sortings.
 */
pair<intT, long> PeelHist(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long, intT> sg_pair = getActiveWedges<long>(ps.wedges_seq_long, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  auto sg_freqs_pair = getSeagullFreqsHist(ps, nu, ps.wedges_seq_long.A, sg_pair.first);
  // Compute updated butterfly counts
  auto ret = getUpdatesHist(ps, sg_freqs_pair.first, sg_freqs_pair.second, nu, butterflies);

  free(sg_freqs_pair.first);
  return make_pair(next_idx, ret);
}

pair<intT, long> PeelHist_seq(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long, intT> sg_pair = getActiveWedges<long>(ps.wedges_seq_long, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  auto sg_freqs_pair = getSeagullFreqsHist_seq(ps, nu, ps.wedges_seq_long.A, sg_pair.first);
  // Compute updated butterfly counts
  auto ret = getUpdatesHist_seq(ps, sg_freqs_pair.first, sg_freqs_pair.second, nu, butterflies);

  free(sg_freqs_pair.first);
  return make_pair(next_idx, ret);
}

/*
 *  Precisely PeelSort, but using hash tables instead of repeated sortings.
 */
pair<intT, long> PeelHash(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);
  // Retrieve all seagulls
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  intT next_idx = getActiveWedgesHash(ps, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges);

  auto wedges_seq_num = ps.wedges_hash->entries_no_init(ps.wedges_seq_intp);

  using T = pair<long,long>;
  granular_for ( j, 0, wedges_seq_num, (wedges_seq_num > 1000), {
      auto sg_freq_pair = ps.wedges_seq_intp.A[j];
      long num = sg_freq_pair.second;
      if (num > 1) {ps.update_hash.insert(T(((long) sg_freq_pair.first % nu), (num * (num - 1)/2)));}
    });

  // Compute updated butterfly counts
  auto num_updates = ps.update_hash.entries_no_init(ps.butterflies_seq_intp);
  ps.resize_update(num_updates);

  granular_for(i, 0, num_updates, (num_updates>1000), {
      auto update_pair = ps.butterflies_seq_intp.A[i];
      butterflies[eltsPerCacheLine*update_pair.first] -= update_pair.second;
      ps.update_seq_int.A[i] = update_pair.first;
    });

  return make_pair(next_idx, num_updates);
}

pair<uintE*, long> PeelOrigParallel(PeelSpace& ps, vertexSubset& active, long* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v,
array_imap<long>& D, buckets<array_imap<long>>& b, uintE k2) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = active.size() < ps.stepSize ? active.size() : ps.stepSize; //tunable parameter

  auto wedges_seq = ps.wedges_seq_int;
  auto used_seq = ps.used_seq_int;

  auto wedges = wedges_seq.A;
  auto used = used_seq.A;
  //long* update_idx = ps.update_idx_seq_int.A;

  const size_t eltsPerCacheLine = 64/sizeof(long);
  granular_for(i,0,nu,nu > 10000, { update_dense[eltsPerCacheLine*i] = false; });

  for(long step = 0; step < (active.size()+stepSize-1)/stepSize; step++) {
    parallel_for_1(long i=step*stepSize; i < min((step+1)*stepSize,active.size()); ++i){
      intT used_idx = 0;
      long shift = nu*(i-step*stepSize);
      intT u_idx = active.vtx(i);
      intT u_offset  = offsetsU[u_idx];
      intT u_deg = offsetsU[u_idx+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1]-offsetsV[v];
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx != u_idx) {
	    if (wedges[shift+u2_idx] > 0) writeAdd(&butterflies[eltsPerCacheLine*u2_idx], (long) -1* ((long)wedges[shift+u2_idx]));
	    // alternatively, at end before clear do this
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) {
        used[shift+used_idx++] = u2_idx;
        if(!update_dense[eltsPerCacheLine*u2_idx]) CAS(&update_dense[eltsPerCacheLine*u2_idx],false,true);
      }
    }
	}
      }
      granular_for(j,0,used_idx,used_idx > 10000, { wedges[shift+used[shift+j]] = 0; });
      //update_idx[i-step*stepSize] = used_idx;
    }
    // used[shift+j] for relevant j as indicated by update contain all of the indices to be updated; must call update here
    /*long use_stepSize = active.size() < (step+1)*stepSize ? active.size() - step*stepSize : stepSize;
    update_idx[use_stepSize] = 0;
    sequence::plusScan(update_idx,update_idx,use_stepSize+1);
    ps.resize_update(update_idx[use_stepSize]);
    auto update_seq = ps.update_seq_int.A;
    parallel_for(long i=step*stepSize; i < min((step+1)*stepSize,active.size()); ++i){
      long shift = nu*(i-step*stepSize);
      granular_for(j,0,update_idx[i+1-step*stepSize] - update_idx[i-step*stepSize],update_idx[i+1-step*stepSize] - update_idx[i-step*stepSize] > 10000, { 
      //parallel_for(long j=0; j < update_idx[i+1-step*stepSize] - update_idx[i-step*stepSize]; ++j) {
        update_seq[update_idx[i-step*stepSize]+j] = used[shift+j];
      });
    }
    updateBuckets(D, b, k2, butterflies, (stepSize < 1000), nu, update_seq, update_idx[use_stepSize]);*/
  }
  auto f = [&] (size_t i) { return update_dense[eltsPerCacheLine*i]; };
  auto f_in = make_in_imap<bool>(nu, f);
  auto out = pbbs::pack_index<uintE>(f_in);
  out.alloc = false;
  return make_pair(out.s, out.size());
}

pair<uintE*, long> PeelOrigParallel_WedgeAware(PeelSpace& ps, vertexSubset& active, long* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v,
array_imap<long>& D, buckets<array_imap<long>>& b, uintE k2) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = active.size() < ps.stepSize ? active.size() : ps.stepSize; //tunable parameter

  auto wedges_seq = ps.wedges_seq_int;
  auto used_seq = ps.used_seq_int;

  auto wedges = wedges_seq.A;
  uintE* used = used_seq.A;
  //long* update_idx = ps.update_idx_seq_int.A;

  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  long* wedgesPrefixSum = getActiveWedgeIdxs(active_map,active.size(),GA,use_v);

  const size_t eltsPerCacheLine = 64/sizeof(long);
  granular_for(i,0,nu,nu > 10000, { update_dense[eltsPerCacheLine*i] = false; });

  for(intT step = 0; step < (active.size()+stepSize-1)/stepSize; step++) {
    std::function<void(intT,intT)> recursive_lambda =
      [&]
      (intT start, intT end){
      if ((start == end-1) || (wedgesPrefixSum[end]-wedgesPrefixSum[start] < 1000)){ 
	for (intT i = start; i < end; i++){
      intT used_idx = 0;
      long shift = nu*(i-step*stepSize);
      intT u_idx = active.vtx(i);
      intT u_offset  = offsetsU[u_idx];
      intT u_deg = offsetsU[u_idx+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1]-offsetsV[v];
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx != u_idx) {
	    if (wedges[shift+u2_idx] > 0) writeAdd(&butterflies[eltsPerCacheLine*u2_idx], (long) -1* ((long)wedges[shift+u2_idx]));
	    // alternatively, at end before clear do this
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) {
        used[shift+used_idx++] = u2_idx; 
        if(!update_dense[eltsPerCacheLine*u2_idx]) CAS(&update_dense[eltsPerCacheLine*u2_idx],false,true);
      }
    }
	}
      }
      granular_for(j,0,used_idx,used_idx > 10000, { wedges[shift+used[shift+j]] = 0; });
      //update_idx[i-step*stepSize] = used_idx;
    }

    } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,active.size()));
    
// used[shift+j] for relevant j as indicated by update contain all of the indices to be updated; must call update here
    /* long use_stepSize = active.size() < (step+1)*stepSize ? active.size() - step*stepSize : stepSize;
    update_idx[use_stepSize] = 0;
    sequence::plusScan(update_idx,update_idx,use_stepSize+1);
    ps.resize_update(update_idx[use_stepSize]);
    auto update_seq = ps.update_seq_int.A;
    parallel_for(long i=step*stepSize; i < min((step+1)*stepSize,active.size()); ++i){
      long shift = nu*(i-step*stepSize);
      granular_for(j,0,update_idx[i+1-step*stepSize] - update_idx[i-step*stepSize],update_idx[i+1-step*stepSize] - update_idx[i-step*stepSize] > 10000, { 
      //parallel_for(long j=0; j < update_idx[i+1-step*stepSize] - update_idx[i-step*stepSize]; ++j) {
        update_seq[update_idx[i-step*stepSize]+j] = used[shift+j];
      });
    }
    updateBuckets(D, b, k2, butterflies, (stepSize < 1000), nu, update_seq, update_idx[use_stepSize]);*/
  }
  free(wedgesPrefixSum);
  auto f = [&] (size_t i) { return update_dense[eltsPerCacheLine*i]; };
  auto f_in = make_in_imap<bool>(nu, f);
  auto out = pbbs::pack_index<uintE>(f_in);
  out.alloc = false;
  return make_pair(out.s, out.size());
}

//***************************************************************************************************
//***************************************************************************************************

void Peel_helper (PeelSpace& ps, vertexSubset& active, long* butterflies, bool* update_dense,
		  bipartiteCSR& GA, bool use_v, long max_wedges, long type, array_imap<long>& D, buckets<array_imap<long>>& b, long k) {
  const long nu = use_v ? GA.nu : GA.nv;
  bool is_seq = (active.size() < 1000);
  if (type == 3) {
    auto ret = PeelOrigParallel(ps, active, butterflies, update_dense, GA, use_v, D, b, k);
    updateBuckets(D, b, k, butterflies, is_seq, nu, ret.first, ret.second);
    free(ret.first);
    return;
  }
  else if (type == 5) {
    auto ret = PeelOrigParallel_WedgeAware(ps, active, butterflies, update_dense, GA, use_v, D, b, k);
    updateBuckets(D, b, k, butterflies, is_seq, nu, ret.first, ret.second);
    free(ret.first);
    return;
  }
  
  long num_wedges = countSeagulls(GA, use_v, active);
  pair<intT, long> ret;
  //TODO remove update dense from all of these
  if (max_wedges >= num_wedges) {
    if (type == 0) ret = PeelHash(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else if (type == 1 && is_seq) ret = PeelSort_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else if (type == 1) ret = PeelSort(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else if (type == 2 && is_seq) ret = PeelHist_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else ret = PeelHist(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);

    updateBuckets(D, b, k, butterflies, is_seq, nu, ps.update_seq_int.A, ret.second);
    ps.clear();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < active.size()) {
    if (type == 0) ret = PeelHash(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == 1 && is_seq) ret = PeelSort_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == 1) ret = PeelSort(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == 2 && is_seq) ret = PeelHist_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else ret = PeelHist(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);

    curr_idx = ret.first;
    updateBuckets(D, b, k, butterflies, is_seq, nu, ps.update_seq_int.A, ret.second);
    ps.clear();
  }
}

array_imap<long> Peel(bipartiteCSR& GA, bool use_v, long* butterflies, long max_wedges, long type, long max_array_size, size_t num_buckets=128) {
  // Butterflies are assumed to be stored on U
  const long nu = use_v ? GA.nu : GA.nv;
  
  const size_t eltsPerCacheLine = 64/sizeof(long);

  auto D = array_imap<long>(nu, [&] (size_t i) { return butterflies[eltsPerCacheLine*i]; });

  auto b = make_buckets(nu, D, increasing, num_buckets); // TODO may also have to fix buckets to use long

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu);
  PeelSpace ps = PeelSpace(type, nu, stepSize);

  bool* update_dense = nullptr;
  if (type == 3 || type == 5) update_dense = newA(bool, eltsPerCacheLine*nu);

  size_t finished = 0;
  //long nonZeroRounds = 0;
  //long totalRounds = 0;
  while (finished != nu) {
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    if (active.size() == 0) {active.del(); continue;}
    long k = bkt.id;
    finished += active.size();
    //totalRounds++;
    //if(active.size() > 0) {nonZeroRounds++; }
    bool is_seq = (active.size() < 1000);

    Peel_helper(ps, active, butterflies, update_dense, GA, use_v, max_wedges, type, D, b, k);

    active.del();
  }
  ps.del();
  b.del();
  if (type == 3 || type == 5) free(update_dense);
  //cout << "totalRounds = " << totalRounds << endl;
  //cout << "nonZeroRounds = " << nonZeroRounds << endl;
  
  return D;
}

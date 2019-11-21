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
#ifndef OPENMP
#include "../../lib/histogram.h"
#include "../../lib/sample_sort.h"
#endif
#include "../../radixsort/RadixSort/radixSort.h"

#include "butterfly_putils.h"
#include "butterfly_utils.h"

// Helper function for seagullSum that sums two-hop degrees
template <class E>
struct seagullSumHelper { 
  uintE u;
  uintT* offsetsU;
  uintT* offsetsV;
  uintE* edgesU;
seagullSumHelper(uintE _u, uintT* _offsetsU, uintT* _offsetsV, uintE* _edgesU) :
  u(_u), offsetsU(_offsetsU), offsetsV(_offsetsV), edgesU(_edgesU) {}
  inline E operator() (const E& i) const {
    intT u_offset = offsetsU[u];
    uintE v = edgesU[u_offset + i];
    return (E) (offsetsV[v+1] - offsetsV[v] - 1);
  }
};

// Computes the number of wedges that each vertex in active contributes
template <class E>
struct seagullSum { 
  uintT* offsetsU;
  uintT* offsetsV;
  uintE* edgesU;
  uintE* active;
seagullSum(uintT* _offsetsU, uintT* _offsetsV, uintE* _edgesU, uintE* _active) :
  offsetsU(_offsetsU), offsetsV(_offsetsV), edgesU(_edgesU), active(_active) {}
  inline E operator() (const E& i) const {
    E u_deg = offsetsU[active[i]+1] - offsetsU[active[i]];
    return sequence::reduce<E>((E) 0, u_deg, addF<E>(), seagullSumHelper<E>(active[i], offsetsU, offsetsV, edgesU));
  }
};

/*
 *  Counts the number of wedges on vertices in active, sequentially.
 * 
 *  GA    : Bipartite graph in CSR format
 *  use_v : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *          produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  active: Set of active vertices
 * 
 *  Returns: Number of wedges on vertices in active
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

/*
 *  Counts the number of wedges on vertices in active, in parallel.
 * 
 *  GA    : Bipartite graph in CSR format
 *  use_v : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *          produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  active: Set of active vertices
 * 
 *  Returns: Number of wedges on vertices in active
 */
long countSeagulls(bipartiteCSR& GA, bool use_v, vertexSubset active) {
  if (active.size() < 1000) return countSeagulls_seq(GA, use_v, active);

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  return sequence::reduce<long>((long) 0, (long) active.size(), addF<long>(),
                                seagullSum<long>(offsetsU, offsetsV, edgesU, active.s));
}


//***************************************************************************************************
//***************************************************************************************************

/*
 *  Retrieve updated butterfly counts given wedges, using sort aggregation.
 * 
 *  nu      : Number of vertices on the bipartition we're peeling from
 *  seagulls: Array of wedges
 *  num_sgs : Number of wedges
 * 
 *  Returns an array of updated butterfly counts associated with the
 *  corresponding vertices, and the number of elements in this array.
 */
pair<tuple<uintE,long>*, long> getSeagullFreqs(const long nu, UVertexPair* seagulls, long num_sgs) {
  using X = tuple<uintE,long>;
  // Sort wedges (considering both active + non-active endpoints), and retrieve frequency counts
#ifdef OPENMP
  sampleSort(seagulls, num_sgs, UVertexPairCmp2());
#else
  pbbs::sample_sort(seagulls, num_sgs, UVertexPairCmp2());
#endif
  pair<long*, long> freq_pair = getFreqs<long>(seagulls, num_sgs, UVertexPairCmp2(), UVertexPairEq(), LONG_MAX,
                                               nonMaxLongF());
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

  return make_pair(sg_freqs,num_sg_freqs);
}

/*
 *  Retrieve updated butterfly counts sequentially given wedges, using sort
 *  aggregation.
 * 
 *  nu      : Number of vertices on the bipartition we're peeling from
 *  seagulls: Array of wedges
 *  num_sgs : Number of wedges
 * 
 *  Returns an array of updated butterfly counts associated with the
 *  corresponding vertices, and the number of elements in this array.
 */
pair<tuple<uintE,long>*, long> getSeagullFreqs_seq(const long nu, UVertexPair* seagulls, long num_sgs) {
#ifdef OPENMP
  sampleSort(seagulls, num_sgs, UVertexPairCmp2());
#else
  pbbs::sample_sort(seagulls, num_sgs, UVertexPairCmp2());
#endif
  return getFreqs_seq<long,uintE>(seagulls, num_sgs, UVertexPairCmp2(), UVertexPairEq(), true,
                                  UVertexPairV2(), choose2(), reflCount<UVertexPair>());
}

/*
 *  Aggregate updated butterfly counts sequentially, using sort aggregation.
 * 
 *  ps          : Holds all array space needed, to be reused between buckets
 *  sg_freqs    : Array of butterfly counts to be removed from corresponding vertex indices. Vertex indices are not
 *                necessarily distinct amongst entries.
 *  num_sg_freqs: Size of sg_freqs
 *  butterflies : Butterfly counts per vertex
 * 
 *  Returns the number of vertices with updated butterfly counts.
 */
long getUpdates_seq(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, long* butterflies) {
  using X = tuple<uintE,long>;
  // Collate all butterflies to be removed with the same non-active endpoint
  // Do this by sorting on the non-active endpoint, and summing the frequencies
#ifdef OPENMP
  sampleSort(sg_freqs, num_sg_freqs, tupleLt<uintE,long>());
#else
  pbbs::sample_sort(sg_freqs, num_sg_freqs, tupleLt<uintE,long>());
#endif
  pair<X*, long> b_freq_pair =
    getFreqs_seq<long,uintE>(sg_freqs, num_sg_freqs, tupleLt<uintE,long>(), tupleEq<uintE,long>(), false,
                             uintETupleGet0(), refl<long>(), uintECount());
  long num_updates = b_freq_pair.second;

  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;
  X* b_updates = b_freq_pair.first;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Remove these butterflies from our array of butterfly counts, and record
  // the vertex whose count has been updated
  for(long i=0; i<num_updates; ++i) {
    uintE u_idx = get<0>(b_updates[i]);
    butterflies[eltsPerCacheLine*u_idx] -= get<1>(b_updates[i]);
    update[i] = u_idx;
  }
  free(b_updates);

  return num_updates;
}

/*
 *  Aggregate updated butterfly counts, using sort aggregation.
 * 
 *  ps          : Holds all array space needed, to be reused between buckets
 *  sg_freqs    : Array of butterfly counts to be removed from corresponding vertex indices. Vertex indices are not
 *                necessarily distinct amongst entries.
 *  num_sg_freqs: Size of sg_freqs
 *  butterflies : Butterfly counts per vertex
 * 
 *  Returns the number of vertices with updated butterfly counts.
 */
long getUpdates(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, long* butterflies) {
  using X = tuple<uintE,long>;
  // Collate all butterflies to be removed with the same non-active endpoint
  // Do this by sorting on the non-active endpoint, and summing the frequencies
#ifdef OPENMP
  sampleSort(sg_freqs, num_sg_freqs, tupleLt<uintE,long>());
#else
  pbbs::sample_sort(sg_freqs, num_sg_freqs, tupleLt<uintE,long>());
#endif
  auto b_freq_pair = getFreqs<long>(sg_freqs, num_sg_freqs, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX,
                                    nonMaxLongF(), false);
  long num_updates = b_freq_pair.second - 1;

  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Remove these butterflies from our array of butterfly counts, and record
  // the vertex whose count has been updated
  parallel_for(long i=1; i < num_updates + 1; ++i) {
    long num_freq = b_freq_pair.first[i] - b_freq_pair.first[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(sg_freqs[b_freq_pair.first[i-1]]), num_freq, tupleAdd<uintE,long>());
    uintE u_idx = get<0>(reduce);
    // Remove these butterflies from our array of butterfly counts
    butterflies[eltsPerCacheLine*u_idx] -= get<1>(reduce);
    // Add the updated index to an update array
    update[i-1] = u_idx;
  }
  free(b_freq_pair.first);

  return num_updates;
}

/*
 *  Retrieve updated butterfly counts given wedges, using histogram aggregation.
 * 
 *  ps      : Holds all array space needed, to be reused between buckets
 *  nu      : Number of vertices on the bipartition we're peeling from
 *  seagulls: Array of wedges
 *  num_sgs : Number of wedges
 * 
 *  Returns an array of updated butterfly counts associated with the
 *  corresponding vertices, and the number of elements in this array.
 */
pair<tuple<uintE,long>*, long> getSeagullFreqsHist(PeelSpace& ps, const long nu, long* seagulls, long num_sgs) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<long,uintE>;
  using T = tuple<uintE,long>;
  pbbsa::sequence<long> sg_seq = pbbsa::sequence<long>(seagulls,num_sgs);

  // Place wedges into a histogram to retrieve frequency counts (considering both active + non-active endpoints)
  tuple<size_t,X*> sg_hist_tuple = pbbsa::sparse_histogram<long,uintE>(sg_seq,nu*nu+nu,ps.tmp,ps.out);
  X* sg_freqs = get<1>(sg_hist_tuple);

  // Filter out any frequency count <= 1, since these won't contribute towards a butterfly
  X* sg_freqs_filter = newA(X,get<0>(sg_hist_tuple));
  size_t num_sg_freqs = sequence::filter(sg_freqs, sg_freqs_filter, get<0>(sg_hist_tuple), greaterOneLongF());

  T* sg_freqs_ret = newA(T, num_sg_freqs);

  // Update our frequency counts to represent the number of butterflies to remove on each
  // non-active endpoint
  // Do this by storing the frequency choose 2 on the non-active endpoint
  parallel_for(long i=0; i < num_sg_freqs; ++i) {
    long num = get<1>(sg_freqs_filter[i]);
    // The non-active endpoint is given by our wedge id % nu, because it is always the second endpoint
    sg_freqs_ret[i] = make_tuple(get<0>(sg_freqs_filter[i]) % nu, num * (num-1) / 2);
  }
  free(sg_freqs_filter);

  return make_pair(sg_freqs_ret, num_sg_freqs);
#endif
}

/*
 *  Retrieve updated butterfly counts sequentially given wedges, using histogram
 *  aggregation.
 * 
 *  ps      : Holds all array space needed, to be reused between buckets
 *  nu      : Number of vertices on the bipartition we're peeling from
 *  seagulls: Array of wedges
 *  num_sgs : Number of wedges
 * 
 *  Returns an array of updated butterfly counts associated with the
 *  corresponding vertices, and the number of elements in this array.
 */
pair<tuple<uintE,long>*, long> getSeagullFreqsHist_seq(PeelSpace& ps, const long nu, long* seagulls, long num_sgs) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<long,uintE>;
  using T = tuple<uintE,long>;
  // Aggregate wedges using a histogram
  pbbsa::sequence<long> sg_seq = pbbsa::sequence<long>(seagulls,num_sgs);
  tuple<size_t,X*> sg_hist_tuple = pbbsa::sparse_histogram<long,uintE>(sg_seq,nu*nu+nu,ps.tmp,ps.out);
  X* sg_freqs = get<1>(sg_hist_tuple);
  size_t num_sg_freqs = get<0>(sg_hist_tuple);

  T* sg_freqs_ret = newA(T, num_sg_freqs);
  long idx = 0;

  // Retrive the butterfly contributions of these wedges
  for(long i=0; i < num_sg_freqs; ++i) {
    long num = get<1>(sg_freqs[i]);
    if (num > 1) {
      sg_freqs_ret[idx] = make_tuple(get<0>(sg_freqs[i]) % nu, num * (num-1)/2);
      idx++;
    }
  }

  return make_pair(sg_freqs_ret, idx);
#endif
}

/*
 *  Aggregate updated butterfly counts sequentially, using histogram aggregation.
 * 
 *  ps          : Holds all array space needed, to be reused between buckets
 *  sg_freqs    : Array of butterfly counts to be removed from corresponding vertex indices. Vertex indices are not
 *                necessarily distinct amongst entries.
 *  num_sg_freqs: Size of sg_freqs
 *  nu          : Number of vertices on the bipartition we're peeling from
 *  butterflies : Butterfly counts per vertex
 * 
 *  Returns the number of vertices with updated butterfly counts.
 */
long getUpdatesHist_seq(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, const long nu,
                        long* butterflies) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<uintE, long>;
  using T = tuple<long, uintE>;

  // Collate all butterflies to be removed with the same non-active endpoint
  // Do this by using a histogram to sum frequencies on the non-active endpoint
  pbbsa::sequence<X> sg_freqs_seq = pbbsa::sequence<X>(sg_freqs,num_sg_freqs);
  tuple<size_t,X*> b_hist_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(sg_freqs_seq,nu,getAdd<uintE,long>,getAddReduce<uintE,long>);
  X* b_freqs = get<1>(b_hist_tuple);
  const size_t eltsPerCacheLine = 64/sizeof(long);
  size_t num_updates = get<0>(b_hist_tuple);

  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;

  // Remove these butterflies from our array of butterfly counts, and record
  // the vertex whose count has been updated
  for(long i=0; i < num_updates; ++i) {
    uintE u_idx = get<0>(b_freqs[i]);
    butterflies[eltsPerCacheLine*u_idx] -=get<1>(b_freqs[i]);
    update[i] = u_idx;
  }
  free(b_freqs);

  return (long) num_updates;
#endif
}

/*
 *  Aggregate updated butterfly counts, using histogram aggregation.
 * 
 *  ps          : Holds all array space needed, to be reused between buckets
 *  sg_freqs    : Array of butterfly counts to be removed from corresponding vertex indices. Vertex indices are not
 *                necessarily distinct amongst entries.
 *  num_sg_freqs: Size of sg_freqs
 *  nu          : Number of vertices on the bipartition we're peeling from
 *  butterflies : Butterfly counts per vertex
 * 
 *  Returns the number of vertices with updated butterfly counts.
 */
long getUpdatesHist(PeelSpace& ps, tuple<uintE,long>* sg_freqs, long num_sg_freqs, const long nu, long* butterflies) {
#ifdef OPENMP
  cout << "Histogram on OPENMP not supported\n";
  exit(0);
#else
  using X = tuple<uintE,long>;
  // Collate all butterflies to be removed with the same non-active endpoint
  // Do this by using a histogram to sum frequencies on the non-active endpoint
  pbbsa::sequence<X> sg_freqs_seq = pbbsa::sequence<X>(sg_freqs,num_sg_freqs);
  tuple<size_t,X*> b_hist_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(sg_freqs_seq,nu,getAdd<uintE,long>,getAddReduce<uintE,long>);
  X* b_freqs = get<1>(b_hist_tuple);
  size_t num_updates = get<0>(b_hist_tuple);

  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Remove these butterflies from our array of butterfly counts, and record
  // the vertex whose count has been updated
  parallel_for(long i=0; i < num_updates; ++i) {
    uintE u_idx = get<0>(b_freqs[i]);
    butterflies[eltsPerCacheLine*u_idx] -= get<1>(b_freqs[i]);
    update[i] = u_idx;
  }
  free(b_freqs);

  return (long) num_updates;
#endif
}

//***************************************************************************************************
//***************************************************************************************************

/*
 *  Peels butterflies given bucket of active vertices, using sort
 *  aggregation.
 * 
 *  ps         : Holds all array space needed, to be reused between buckets
 *  active     : Set of active vertices
 *  butterflies: Butterfly counts per vertex
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *               produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges : Number of wedges produced by the active set
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  curr_idx   : Denotes the starting vertex in this batch
 * 
 *  Returns the next vertex to process for the next batch, and the number of
 *  vertices with updated butterfly counts.
 */
pair<long, long> PeelSort(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v,
                          long num_wedges, long max_wedges, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  // Retrieve all wedges on active vertices
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto sg_pair = getActiveWedges<UVertexPair>(ps.wedges_seq_uvp, active_map, active.size(), GA, use_v,
                                              UVertexPairCons(), max_wedges, curr_idx, num_wedges);
  long next_idx = sg_pair.second;

  // Aggregate wedges by endpoints by sorting
  auto sg_freqs_pair = getSeagullFreqs(nu, ps.wedges_seq_uvp.A, sg_pair.first);

  // Compute updated butterfly counts
  auto ret = getUpdates(ps, sg_freqs_pair.first, sg_freqs_pair.second, butterflies);

  free(sg_freqs_pair.first);
  return make_pair(next_idx, ret);
}

/*
 *  Peels butterflies sequentially given bucket of active vertices, using
 *  sort aggregation.
 * 
 *  ps         : Holds all array space needed, to be reused between buckets
 *  active     : Set of active vertices
 *  butterflies: Butterfly counts per vertex
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *               produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges : Number of wedges produced by the active set
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  curr_idx   : Denotes the starting vertex in this batch
 * 
 *  Returns the next vertex to process for the next batch, and the number of
 *  vertices with updated butterfly counts.
 */
pair<long, long> PeelSort_seq(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v,
                              long num_wedges, long max_wedges, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  // Retrieve all wedges on active vertices
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto sg_pair = getActiveWedges<UVertexPair>(ps.wedges_seq_uvp, active_map, active.size(), GA, use_v,
                                              UVertexPairCons(), max_wedges, curr_idx, num_wedges);
  long next_idx = sg_pair.second;

  // Aggregate wedges by endpoints by sorting
  auto sg_freqs_pair = getSeagullFreqs_seq(nu, ps.wedges_seq_uvp.A, sg_pair.first);

  // Compute updated butterfly counts
  auto ret = getUpdates_seq(ps, sg_freqs_pair.first, sg_freqs_pair.second, butterflies);

  free(sg_freqs_pair.first);
  return make_pair(next_idx, ret);
}

/*
 *  Peels butterflies given bucket of active vertices, using histogram
 *  aggregation.
 * 
 *  ps         : Holds all array space needed, to be reused between buckets
 *  active     : Set of active vertices
 *  butterflies: Butterfly counts per vertex
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *               produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges : Number of wedges produced by the active set
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  curr_idx   : Denotes the starting vertex in this batch
 * 
 *  Returns the next vertex to process for the next batch, and the number of
 *  vertices with updated butterfly counts.
 */
pair<long, long> PeelHist(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v,
                          long num_wedges, long max_wedges, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  // Retrieve all wedges on active vertices
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto sg_pair = getActiveWedges<long>(ps.wedges_seq_long, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu),
                                       max_wedges, curr_idx, num_wedges);
  long next_idx = sg_pair.second;

  // Aggregate wedges by endpoints using a histogram
  auto sg_freqs_pair = getSeagullFreqsHist(ps, nu, ps.wedges_seq_long.A, sg_pair.first);

  // Compute updated butterfly counts
  auto ret = getUpdatesHist(ps, sg_freqs_pair.first, sg_freqs_pair.second, nu, butterflies);

  free(sg_freqs_pair.first);
  return make_pair(next_idx, ret);
}

/*
 *  Peels butterflies sequentially given bucket of active vertices, using
 *  histogram aggregation.
 * 
 *  ps         : Holds all array space needed, to be reused between buckets
 *  active     : Set of active vertices
 *  butterflies: Butterfly counts per vertex
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *               produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges : Number of wedges produced by the active set
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  curr_idx   : Denotes the starting vertex in this batch
 * 
 *  Returns the next vertex to process for the next batch, and the number of
 *  vertices with updated butterfly counts.
 */
pair<long, long> PeelHist_seq(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v,
                              long num_wedges, long max_wedges, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  // Retrieve all wedges on active vertices
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto sg_pair = getActiveWedges<long>(ps.wedges_seq_long, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu),
                                       max_wedges, curr_idx, num_wedges);
  long next_idx = sg_pair.second;

  // Aggregate wedges by endpoints using a histogram
  auto sg_freqs_pair = getSeagullFreqsHist_seq(ps, nu, ps.wedges_seq_long.A, sg_pair.first);

  // Compute updated butterfly counts
  auto ret = getUpdatesHist_seq(ps, sg_freqs_pair.first, sg_freqs_pair.second, nu, butterflies);

  free(sg_freqs_pair.first);

  return make_pair(next_idx, ret);
}

/*
 *  Peels butterflies given bucket of active vertices, using hash aggregation.
 * 
 *  ps         : Holds all array space needed, to be reused between buckets
 *  active     : Set of active vertices
 *  butterflies: Butterfly counts per vertex
 *  GA         : Bipartite graph in CSR format
 *  use_v      : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *               produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges : Number of wedges produced by the active set
 *  max_wedges : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *               memory
 *  curr_idx   : Denotes the starting vertex in this batch
 * 
 *  Returns the next vertex to process for the next batch, and the number of
 *  vertices with updated butterfly counts.
 */
pair<long, long> PeelHash(PeelSpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, bool use_v,
                          long num_wedges, long max_wedges, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Hash wedges on active vertices, by endpoints
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  long next_idx = getActiveWedgesHash(ps, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu), max_wedges,
                                      curr_idx, num_wedges);

  auto wedges_seq_num = ps.wedges_hash->entries_no_init(ps.wedges_seq_intp);

  using T = pair<long,long>;
  // Compute butterflies on active vertices and store in another hash table
  granular_for ( j, 0, wedges_seq_num, (wedges_seq_num > 1000), {
      auto sg_freq_pair = ps.wedges_seq_intp.A[j];
      long num = sg_freq_pair.second;
      if (num > 1) {ps.update_hash.insert(T(((long) sg_freq_pair.first % nu), (num * (num - 1)/2)));}
    });
  auto num_updates = ps.update_hash.entries_no_init(ps.butterflies_seq_intp);
  ps.resize_update(num_updates);

  // Update butterfly counts for those vertices
  granular_for(i, 0, num_updates, (num_updates>1000), {
      auto update_pair = ps.butterflies_seq_intp.A[i];
      butterflies[eltsPerCacheLine*update_pair.first] -= update_pair.second;
      ps.update_seq_int.A[i] = update_pair.first;
    });

  return make_pair(next_idx, num_updates);
}

/*
 *  Peels butterflies given bucket of active vertices, using simple batching.
 * 
 *  ps            : Holds all array space needed, to be reused between buckets
 *  active        : Set of active vertices
 *  butterflies   : Butterfly counts per vertex
 *  update_dense  : Keeps track of vertices with updated butterfly counts
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 * 
 *  Returns an array of vertices with updated butterfly counts and the size of
 *  that array.
 */
pair<uintE*, long> PeelOrigParallel(PeelSpace& ps, vertexSubset& active, long* butterflies, bool* update_dense,
                                    bipartiteCSR& GA, bool use_v) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = active.size() < ps.stepSize ? active.size() : ps.stepSize; // tunable parameter

  auto wedges_seq = ps.wedges_seq_int;
  auto used_seq = ps.used_seq_int;

  // Stores number of wedges per endpoints
  auto wedges = wedges_seq.A;
  // Keeps track of entries used in the wedges array
  auto used = used_seq.A;

  const size_t eltsPerCacheLine = 64/sizeof(long);
  granular_for(i,0,nu,nu > 10000, { update_dense[eltsPerCacheLine*i] = false; });

  // Consider active vertices i in batches
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
	// Iterate through all two-hop neighbors of u_idx
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx != u_idx) {
	    // Subtract the contributing butterfly count from u2_idx
	    if (wedges[shift+u2_idx] > 0)
	      writeAdd(&butterflies[eltsPerCacheLine*u2_idx], (long) -1* ((long)wedges[shift+u2_idx]));
	    // Increment number of wedges on second endpoint
	    wedges[shift+u2_idx]++;
	    // Keep track of used second endpoints
	    if (wedges[shift+u2_idx] == 1) {
	      used[shift+used_idx++] = u2_idx;
	      // Keep track that u2_idx has an updated butterfly count
	      if(!update_dense[eltsPerCacheLine*u2_idx]) CAS(&update_dense[eltsPerCacheLine*u2_idx], false, true);
	    }
	  }
	}
      }
      // Clear wedges array for reuse
      granular_for(j,0,used_idx,used_idx > 10000, { wedges[shift+used[shift+j]] = 0; });
    }
  }

  // Collate edges with updated butterfly counts
  auto f = [&] (size_t i) { return update_dense[eltsPerCacheLine*i]; };
  auto f_in = make_in_imap<bool>(nu, f);
  auto out = pbbs::pack_index<uintE>(f_in);
  out.alloc = false;

  return make_pair(out.s, out.size());
}

/*
 *  Peels butterflies given bucket of active vertices, using wedge-aware batching.
 * 
 *  ps            : Holds all array space needed, to be reused between buckets
 *  active        : Set of active vertices
 *  butterflies   : Butterfly counts per vertex
 *  update_dense  : Keeps track of vertices with updated butterfly counts
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 * 
 *  Returns an array of vertices with updated butterfly counts and the size of
 *  that array.
 */
pair<uintE*, long> PeelOrigParallel_WedgeAware(PeelSpace& ps, vertexSubset& active, long* butterflies,
                                               bool* update_dense, bipartiteCSR& GA, bool use_v) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = active.size() < ps.stepSize ? active.size() : ps.stepSize; // tunable parameter

  auto wedges_seq = ps.wedges_seq_int;
  auto used_seq = ps.used_seq_int;

  // Stores number of wedges per endpoints
  auto wedges = wedges_seq.A;
  // Keeps track of entries used in the wedges array
  uintE* used = used_seq.A;

  // Compute the indices of wedges given by the active subset, so that batches can be set accordingly
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  long* wedgesPrefixSum = getActiveWedgeIdxs(active_map,active.size(),GA,use_v);

  const size_t eltsPerCacheLine = 64/sizeof(long);
  granular_for(i,0,nu,nu > 10000, { update_dense[eltsPerCacheLine*i] = false; });

  // Consider active vertices i in batches, as given by wedgesPrefixSum
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
	    // Iterate through all two-hop neighbors of u_idx
	    for (intT k=0; k < v_deg; ++k) { 
	      uintE u2_idx = edgesV[v_offset+k];
	      if (u2_idx != u_idx) {
		// Subtract the contributing butterfly count from u2_idx
		if (wedges[shift+u2_idx] > 0)
		  writeAdd(&butterflies[eltsPerCacheLine*u2_idx], (long) -1* ((long)wedges[shift+u2_idx]));
		// Increment number of wedges on second endpoint
		wedges[shift+u2_idx]++;
		// Keep track of used second endpoints
		if (wedges[shift+u2_idx] == 1) {
		  used[shift+used_idx++] = u2_idx;
		  // Keep track that u2_idx has an updated butterfly count
		  if(!update_dense[eltsPerCacheLine*u2_idx]) CAS(&update_dense[eltsPerCacheLine*u2_idx], false, true);
		}
	      }
	    }
	  }
	  // Clear wedges array for reuse
	  granular_for(j,0,used_idx,used_idx > 10000, { wedges[shift+used[shift+j]] = 0; });
	}
      } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,active.size()));

  }
  free(wedgesPrefixSum);

  // Collate edges with updated butterfly counts
  auto f = [&] (size_t i) { return update_dense[eltsPerCacheLine*i]; };
  auto f_in = make_in_imap<bool>(nu, f);
  auto out = pbbs::pack_index<uintE>(f_in);
  out.alloc = false;

  return make_pair(out.s, out.size());
}

//***************************************************************************************************
//***************************************************************************************************

/*
 *  Peels butterflies given bucket of active vertices, and updates the bucketing
 *  structure appropriately.
 * 
 *  ps          : Holds all array space needed, to be reused between buckets
 *  active      : Set of active vertices
 *  butterflies : Butterfly counts per edge
 *  update_dense: Depending on aggregation type, keeps track of vertices with updated butterfly counts
 *  GA          : Bipartite graph in CSR format
 *  use_v       : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges  : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                memory
 *  type        : Wedge aggregation type
 *  D           : Map of vertices to the buckets that they're in
 *  b           : Bucketing structure
 *  k           : ID of current bucket
 */
void Peel_helper (PeelSpace& ps, vertexSubset& active, long* butterflies, bool* update_dense,
                  bipartiteCSR& GA, bool use_v, long max_wedges, PeelType type, array_imap<long>& D,
                  buckets<array_imap<long>>& b, long k) {
  const long nu = use_v ? GA.nu : GA.nv;
  bool is_seq = (active.size() < 1000);

  // Choose wedge aggregation type
  if (type == PBATCHS) {
    auto ret = PeelOrigParallel(ps, active, butterflies, update_dense, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, nu, ret.first, ret.second);
    free(ret.first);
    return;
  }
  else if (type == PBATCHWA) {
    auto ret = PeelOrigParallel_WedgeAware(ps, active, butterflies, update_dense, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, nu, ret.first, ret.second);
    free(ret.first);
    return;
  }
  
  // Compute the number of wedges given by the active subset
  long num_wedges = countSeagulls(GA, use_v, active);
  pair<long, long> ret;

  if (max_wedges >= num_wedges) {
    if (type == PHASH) ret = PeelHash(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else if (type == PSORT && is_seq) ret = PeelSort_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else if (type == PSORT) ret = PeelSort(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else if (type == PHIST && is_seq) ret = PeelHist_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else ret = PeelHist(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);

    updateBuckets(D, b, k, butterflies, is_seq, nu, ps.update_seq_int.A, ret.second);
    ps.clear();
    return;
  }
  long curr_idx = 0;
  while(curr_idx < active.size()) {
    if (type == PHASH) ret = PeelHash(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == PSORT && is_seq)
      ret = PeelSort_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == PSORT) ret = PeelSort(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == PHIST && is_seq)
      ret = PeelHist_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else ret = PeelHist(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);

    curr_idx = ret.first;
    updateBuckets(D, b, k, butterflies, is_seq, nu, ps.update_seq_int.A, ret.second);
    ps.clear();
  }
}

/*
 *  Peels butterflies per vertex.
 * 
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  butterflies   : Butterfly counts per vertex
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  type          : Wedge aggregation type
 *  max_array_size:
 *  num_buckets   : Number of buckets to initialize bucketing structure with
 * 
 *  Returns a map of vertices to the buckets they were peeled in.
 */
array_imap<long> Peel(bipartiteCSR& GA, bool use_v, long* butterflies, long max_wedges, PeelType type,
                      long max_array_size, size_t num_buckets=128) {
  // Butterflies are stored on GA.U if use_v is true, and on GA.V otherwise
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Construct initial buckets based on butterfly counts
  auto D = array_imap<long>(nu, [&] (size_t i) { return butterflies[eltsPerCacheLine*i]; });
  auto b = make_buckets(nu, D, increasing, num_buckets); // TODO may also have to fix buckets to use long
  
  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu);
  // Set up space for peeling algorithms
  PeelSpace ps = PeelSpace(type, nu, stepSize);
  bool* update_dense = nullptr;
  if (type == PBATCHS || type == PBATCHWA) update_dense = newA(bool, eltsPerCacheLine*nu);

  size_t finished = 0;
  // Peel each bucket
  while (finished != nu) {
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    if (active.size() == 0) {active.del(); continue;}
    long k = bkt.id;
    finished += active.size();
    bool is_seq = (active.size() < 1000);

    // Peel butterflies
    Peel_helper(ps, active, butterflies, update_dense, GA, use_v, max_wedges, type, D, b, k);

    active.del();
  }

  ps.del();
  b.del();
  if (type == PBATCHS || type == PBATCHWA) free(update_dense);

  return D;
}

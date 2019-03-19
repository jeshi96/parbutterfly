#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#define HYPER 1
#define LONG 1
#define EDGELONG 1
#define MCX16 1

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
#include "../../lib/gbbs-histogram.h"

#include "butterfly_count.h"
#include "butterfly_ecount.h"
#include "butterfly_peel.h"

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }

#include <vector>

using namespace std;

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
      uintT v = edgesU[u_offset + j];
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
pair<tuple<uintE,uintE>*, long> getSeagullFreqs(const long nu, UVertexPair* seagulls, long num_sgs) {
  using X = tuple<uintE,uintE>;
  // Sort seagulls (considering both active + non-active endpoints), and retrieve frequency counts
  pair<uintE*, long> freq_pair = getFreqs(seagulls, num_sgs, UVertexPairCmp2(), UVertexPairEq());
  long num_sg_freqs = freq_pair.second - 1;
  X* sg_freqs = newA(X, num_sg_freqs);
  // When retrieving frequency counts, store the frequency choose 2 with the non-active endpoint
  // This gives us the number of butterflies to be removed on the non-active endpoint
  parallel_for(long i=1; i < freq_pair.second; ++i) {
    uintE num = freq_pair.first[i] - freq_pair.first[i-1];
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

pair<tuple<uintE,uintE>*, long> getSeagullFreqs_seq(const long nu, UVertexPair* seagulls, long num_sgs) {
  return getFreqs_seq<uintE>(seagulls, num_sgs, UVertexPairCmp2(), UVertexPairEq(), true,
    UVertexPairV2(), choose2(), reflCount<UVertexPair>());
}

long getUpdates_seq(tuple<uintE,uintE>* sg_freqs, long num_sg_freqs, uintE* butterflies, bool* update_dense) {
  using X = tuple<uintE,uintE>;

  pair<X*, long> b_freq_pair = getFreqs_seq<uintE>(sg_freqs, num_sg_freqs, uintETupleLt(), uintETupleEq(), false,
    uintETupleGet0(), refl<uintE>(), uintECount());
  long num_updates = b_freq_pair.second;
  X* b_updates = b_freq_pair.first;
  const intT eltsPerCacheLine = 64/sizeof(long);

  for(long i=0; i<num_updates; ++i) {
    uintE u_idx = get<0>(b_updates[i]);
    butterflies[eltsPerCacheLine*u_idx] -= get<1>(b_updates[i]);
    update_dense[u_idx] = get<1>(b_updates[i]) > 0 ? true : false; //this is b/c seq doesn't filter
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
long getUpdates(tuple<uintE,uintE>* sg_freqs, long num_sg_freqs, uintE* butterflies, bool* update_dense) {
  using X = tuple<uintE,uintE>;
  // Now, collate all butterflies to be removed with the same non-active endpoint
  // Do this by sorting on the non-active endpoint, and summing the frequencies
  // TODO MAKE SURE THIS IS RIGHT LOL
  pair<uintE*, long> b_freq_pair = getFreqs(sg_freqs, num_sg_freqs, uintETupleLt(), uintETupleEq(), false);
  long num_updates = b_freq_pair.second - 1;
  const intT eltsPerCacheLine = 64/sizeof(long);

  //uintE* update = newA(uintE, num_updates);
  parallel_for(long i=1; i < num_updates + 1; ++i) {
    uintE num_freq = b_freq_pair.first[i] - b_freq_pair.first[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(sg_freqs[b_freq_pair.first[i-1]]), num_freq, uintETupleAdd());
    uintE u_idx = get<0>(reduce);
    // Remove these butterflies from our array of butterfly counts
    butterflies[eltsPerCacheLine*u_idx] -= get<1>(reduce);
    // Add the updated index to an update array, to be returned
    update_dense[u_idx] = get<1>(reduce) > 0 ? true : false;
  }
  free(b_freq_pair.first);

  return num_updates;
}

/*
 *  Precisely getSeagullFreqs, but using histograms instead of sorting.
 */
pair<tuple<uintE,uintE>*, long> getSeagullFreqsHist(const long nu, uintE* seagulls, long num_sgs) {
  using X = tuple<uintE,uintE>;
  // TODO integrate sequence into histogram code (so we don't have to convert?)
  pbbsa::sequence<uintE> sg_seq = pbbsa::sequence<uintE>(seagulls,num_sgs);

  // Place seagulls into a histogram to retrieve frequency counts (considering both active + non-active endpoints)
  tuple<size_t,X*> sg_hist_tuple = pbbsa::sparse_histogram<uintE,uintE>(sg_seq,nu*nu+nu);
  X* sg_freqs = get<1>(sg_hist_tuple);
  // Filter out any frequency count <= 1, since these won't contribute towards a butterfly
  X* sg_freqs_filter = newA(X,get<0>(sg_hist_tuple));
  size_t num_sg_freqs = sequence::filter(sg_freqs, sg_freqs_filter, get<0>(sg_hist_tuple), greaterOneF());
  //free(sg_freqs);

  // Update our frequency counts to represent the number of butterflies to remove on each
  // non-active endpoint
  // Do this by storing the frequency choose 2 on the non-active endpoint
  parallel_for(long i=0; i < num_sg_freqs; ++i) {
    uintE num = get<1>(sg_freqs_filter[i]);
    // The non-active endpoint is given by our seagull id % nu, because it is always the second endpoint
    sg_freqs_filter[i] = make_tuple(get<0>(sg_freqs_filter[i]) % nu, num * (num-1) / 2);
  }
  return make_pair(sg_freqs_filter, num_sg_freqs);
}

pair<tuple<uintE,uintE>*, long> getSeagullFreqsHist_seq(const long nu, uintE* seagulls, long num_sgs) {
  using X = tuple<uintE,uintE>;
  pbbsa::sequence<uintE> sg_seq = pbbsa::sequence<uintE>(seagulls,num_sgs);
  tuple<size_t,X*> sg_hist_tuple = pbbsa::sparse_histogram<uintE,uintE>(sg_seq,nu*nu+nu);
  X* sg_freqs = get<1>(sg_hist_tuple);
  size_t num_sg_freqs = get<0>(sg_hist_tuple);
  long idx = 0;
  for(long i=0; i < num_sg_freqs; ++i) {
    uintE num = get<1>(sg_freqs[i]);
    if (num > 1) {
      sg_freqs[idx] = make_tuple(get<0>(sg_freqs[i]) % nu, num * (num-1)/2);
      idx++;
    }
  }
  return make_pair(sg_freqs, idx);
}

// TODO just use granular for on getUpdatesHis?
long getUpdatesHist_seq(tuple<uintE,uintE>* sg_freqs, long num_sg_freqs, const long nu, uintE* butterflies, bool* update_dense) {
  using X = tuple<uintE, uintE>;
  pbbsa::sequence<X> sg_freqs_seq = pbbsa::sequence<X>(sg_freqs,num_sg_freqs);
  tuple<size_t,X*> b_hist_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(sg_freqs_seq,nu,getAdd<uintE,uintE>,getAddReduce<uintE,uintE>);
  X* b_freqs = get<1>(b_hist_tuple);
  const intT eltsPerCacheLine = 64/sizeof(long);
  size_t num_updates = get<0>(b_hist_tuple);
  for(long i=0; i < num_updates; ++i) {
    uintE u_idx = get<0>(b_freqs[i]);
    butterflies[eltsPerCacheLine*u_idx] -=get<1>(b_freqs[i]);
    update_dense[u_idx] = true;
  }
  free(b_freqs);
  return (long)num_updates;
}

/*
 *  Precisely getUpdates, but using histograms instead of sorting.
 */
long getUpdatesHist(tuple<uintE,uintE>* sg_freqs, long num_sg_freqs, const long nu, uintE* butterflies, bool* update_dense) {
  using X = tuple<uintE,uintE>;
  // Now, collate all butterflies to be removed with the same non-active endpoint
  // Do this by using a histogram to sum frequencies on the non-active endpoint
  pbbsa::sequence<X> sg_freqs_seq = pbbsa::sequence<X>(sg_freqs,num_sg_freqs);
  tuple<size_t,X*> b_hist_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(sg_freqs_seq,nu,getAdd<uintE,uintE>,getAddReduce<uintE,uintE>);
  X* b_freqs = get<1>(b_hist_tuple);
  size_t num_updates = get<0>(b_hist_tuple);
  const intT eltsPerCacheLine = 64/sizeof(long);

  // Remove these butterflies from our array of butterfly counts, and change our frequency
  // array to be an update array on the new number of butterflies on our non-active endpoint
  parallel_for(long i=0; i < num_updates; ++i) {
    uintE u_idx = get<0>(b_freqs[i]);
    butterflies[eltsPerCacheLine*u_idx] -=get<1>(b_freqs[i]);
    update_dense[u_idx] = true;
  }
  free(b_freqs);

  return (long)num_updates;
}


/*
 *  Precisely getSeagullFreqs, but using hash tables instead of sorting.
 */
sparseAdditiveSet<uintE> getSeagullFreqsHash(sparseAdditiveSet<uintE>* seagulls, const long nu) {
  sparseAdditiveSet<uintE> update_hash = sparseAdditiveSet<uintE>(nu, (float) 1, UINT_E_MAX);
  _seq<pair<uintE,uintE>> sg_freqs = seagulls->entries();

  using T = pair<uintE,uintE>;

  //parallel_for (long j=0; j < sg_freqs.n; ++j) {
  granular_for ( j, 0, sg_freqs.n, (sg_freqs.n > 1000), {
    T sg_freq_pair = sg_freqs.A[j];
    uintE num = sg_freq_pair.second;
    if (num > 1) {update_hash.insert(T(sg_freq_pair.first % nu, num * (num - 1)/2));}
  });
  sg_freqs.del();
  return update_hash;
}

/*
 *  Precisely getUpdates, but using hash tables instead of sorting.
 */
long getUpdatesHash(sparseAdditiveSet<uintE>& update_hash, uintE* butterflies, bool* update_dense) {

  _seq<pair<uintE,uintE>> update_seq = update_hash.entries();
  long num_updates = update_seq.n;
  //uintE* update = newA(uintE, num_updates);

  using T = pair<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  //parallel_for (long i=0; i < num_updates; ++i) {
  granular_for(i, 0, num_updates, (num_updates>1000), {
    T update_pair = update_seq.A[i];
    uintE u_idx = update_pair.first;
    butterflies[eltsPerCacheLine*u_idx] -= update_pair.second;
    //update[i] = u_idx;
    update_dense[u_idx] = true;
  });

  update_seq.del();
  return num_updates;
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
intT PeelSort(_seq<UVertexPair>& wedges_seq, vertexSubset& active, uintE* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long,intT> sg_pair = getActiveWedges<UVertexPair>(wedges_seq, active_map, active.size(), GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  pair<tuple<uintE,uintE>*, long> sg_freqs_pair = getSeagullFreqs(nu, wedges_seq.A, sg_pair.first);

  // Compute updated butterfly counts
  getUpdates(sg_freqs_pair.first, sg_freqs_pair.second, butterflies, update_dense);

  free(sg_freqs_pair.first);
  return next_idx;
}

//TODO these are all CE versions -- do version where we don't use getUpdates_seq, and we writeadd (ish) to update
intT PeelSort_seq(_seq<UVertexPair>& wedges_seq, vertexSubset& active, uintE* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long, intT> sg_pair = getActiveWedges<UVertexPair>(wedges_seq, active_map, active.size(), GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  pair<tuple<uintE,uintE>*, long> sg_freqs_pair = getSeagullFreqs_seq(nu, wedges_seq.A, sg_pair.first);

  // Compute updated butterfly counts
  getUpdates_seq(sg_freqs_pair.first, sg_freqs_pair.second, butterflies, update_dense);

  free(sg_freqs_pair.first);
  return next_idx;
}

/*
 *  Precisely PeelSort, but using histograms instead of repeated sortings.
 */
intT PeelHist(_seq<uintE>& wedges_seq, vertexSubset& active, uintE* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long, intT> sg_pair = getActiveWedges<uintE>(wedges_seq, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  pair<tuple<uintE,uintE>*, long> sg_freqs_pair = getSeagullFreqsHist(nu, wedges_seq.A, sg_pair.first);
  // Compute updated butterfly counts
  getUpdatesHist(sg_freqs_pair.first, sg_freqs_pair.second, nu, butterflies, update_dense);

  free(sg_freqs_pair.first);
  return next_idx;
}

intT PeelHist_seq(_seq<uintE>& wedges_seq, vertexSubset& active, uintE* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  // Retrieve all seagulls
  const long nu = use_v ? GA.nu : GA.nv;
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  pair<long, intT> sg_pair = getActiveWedges<uintE>(wedges_seq, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges);
  intT next_idx = sg_pair.second;

  pair<tuple<uintE,uintE>*, long> sg_freqs_pair = getSeagullFreqsHist_seq(nu, wedges_seq.A, sg_pair.first);
  // Compute updated butterfly counts
  getUpdatesHist_seq(sg_freqs_pair.first, sg_freqs_pair.second, nu, butterflies, update_dense);

  free(sg_freqs_pair.first);
  return next_idx;
}

/*
 *  Precisely PeelSort, but using hash tables instead of repeated sortings.
 */
intT PeelHash(PeelSpace& ps, vertexSubset& active, uintE* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  // Retrieve all seagulls
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  intT next_idx = getActiveWedgesHash(ps, active_map, active.size(), GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges);

  sparseAdditiveSet<uintE> update_hash = getSeagullFreqsHash(ps.wedges_hash, nu);

  // Compute updated butterfly counts]
  getUpdatesHash(update_hash, butterflies, update_dense);

  update_hash.del();
  return next_idx;
}

//***************************************************************************************************
//***************************************************************************************************

void Peel_helper (PeelSpace& ps, vertexSubset& active, uintE* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v, long max_wedges, long type) {
  const long nu = use_v ? GA.nu : GA.nv;
  bool is_seq = (active.size() < 1000);
  long num_wedges = countSeagulls(GA, use_v, active);

  if (max_wedges >= num_wedges) {
    if (type == 0) PeelHash(ps, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges);
    else if (type == 1 && is_seq) PeelSort_seq(ps.wedges_seq_uvp, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges);
    else if (type == 1) PeelSort(ps.wedges_seq_uvp, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges);
    else if (type == 2 && is_seq) PeelHist_seq(ps.wedges_seq_int, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges);
    else PeelHist(ps.wedges_seq_int, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges);
    ps.clear();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < active.size()) {
    if (type == 0) curr_idx = PeelHash(ps, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == 1 && is_seq) curr_idx = PeelSort_seq(ps.wedges_seq_uvp, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == 1) curr_idx = PeelSort(ps.wedges_seq_uvp, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges, curr_idx);
  	else if (type == 2 && is_seq) curr_idx = PeelHist_seq(ps.wedges_seq_int, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges, curr_idx);
    else curr_idx = PeelHist(ps.wedges_seq_int, active, butterflies, update_dense, GA, use_v, num_wedges, max_wedges, curr_idx);
    ps.clear();
  }
}

array_imap<uintE> Peel(bipartiteCSR& GA, bool use_v, uintE* butterflies, long max_wedges, long type=0, size_t num_buckets=128) {
  // Butterflies are assumed to be stored on U
  const long nu = use_v ? GA.nu : GA.nv;
  
  using X = tuple<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  auto D = array_imap<uintE>(nu, [&] (size_t i) { return butterflies[eltsPerCacheLine*i]; });

  auto b = make_buckets(nu, D, increasing, num_buckets);

  bool* update_dense = newA(bool, nu);
  PeelSpace ps = PeelSpace(type, nu);

  size_t finished = 0;
  while (finished != nu) {
  //cout <<"get1, ";fflush(stdout);
    auto bkt = b.next_bucket();
  //cout <<"get2, ";fflush(stdout);
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();
    bool is_seq = (active.size() < 1000);

    parallel_for(intT i=0; i < nu; ++i) { update_dense[i] = false; }
    Peel_helper(ps, active, butterflies, update_dense, GA, use_v, max_wedges, type);

    auto f = [&] (size_t i) { return update_dense[i]; };
    auto f_in = make_in_imap<bool>(nu, f);
    auto out = pbbs::pack_index<uintE>(f_in);
    out.alloc = false;

    pair<tuple<uintE,uintE>*,long> bucket_pair;
    if (is_seq) bucket_pair = updateBuckets_seq(out.s, out.size(), butterflies, D, b, k);
    else bucket_pair = updateBuckets(out.s, out.size(), butterflies, D, b, k);

    free(out.s);

    vertexSubsetData<uintE> moved = vertexSubsetData<uintE>(nu, bucket_pair.second, bucket_pair.first);
    b.update_buckets(moved.get_fn_repr(), moved.size());

    moved.del(); active.del();
  }
  ps.del();
  return D;
}

void CountOrigCompact(bipartiteCSR& GA, bool use_v) {
  timer t1,t2;
  t1.start();
  //cout << GA.nv << " " << GA.nu << " " << GA.numEdges << endl;
  cout << "Original Serial (make sure running with CILK_NWORKERS=1)" << endl;  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long results = 0;
  uintE* wedges = newA(uintE, nu);
  uintE* used = newA(uintE, nu);

  for(intT i=0; i < nu; ++i) { wedges[i] = 0; }

  uintE* butterflies = newA(uintE,nu);
  for(intT i=0; i < nu; ++i) { butterflies[i] = 0; }

  t1.reportTotal("preprocess");
  t2.start();

  for(intT i=0; i < nu; ++i){
    intT used_idx = 0;
    intT u_offset  = offsetsU[i];
    intT u_deg = offsetsU[i+1]-u_offset;
    for (intT j=0; j < u_deg; ++j ) {
      intT v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1]-offsetsV[v];
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2_idx = edgesV[v_offset+k];
        if (u2_idx < i) {
          //butterflies[i] += wedges[u2_idx];
          //butterflies[u2_idx] += wedges[u2_idx];
          results += wedges[u2_idx];
          wedges[u2_idx]++;
          if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
        }
        else break;
      }
    }
    for(intT j=0; j < used_idx; ++j) { wedges[used[j]] = 0; }
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);

  for(long i=0;i<nu;i++) results += butterflies[i];
  free(butterflies);
  cout << "num: " << results << "\n";
}

void CountOrigCompactParallel(bipartiteCSR& GA, bool use_v) {
  timer t1,t2;
  t1.start();
  cout << "Original Parallel" << endl;
  //cout << GA.nv << " " << GA.nu << " " << GA.numEdges << endl;
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = 1000; //tunable parameter

  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);

  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);
  
  long* results = newA(long,eltsPerCacheLine*stepSize); //one entry per cache line
  granular_for(i,0,stepSize,stepSize > 10000, {results[eltsPerCacheLine*i] = 0;});
  
  uintE* butterflies = newA(uintE,nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[i] = 0; });
  t1.reportTotal("preprocess");

  t2.start();

  for(intT step = 0; step < (nu+stepSize-1)/stepSize; step++) {
    parallel_for_1(intT i=step*stepSize; i < min((step+1)*stepSize,nu); ++i){
      intT used_idx = 0;
      intT shift = nu*(i-step*stepSize);
      intT u_offset  = offsetsU[i];
      intT u_deg = offsetsU[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	intT v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1]-offsetsV[v];
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx < i) {
	    //butterflies[i] += wedges[u2_idx];
	    //butterflies[u2_idx] += wedges[u2_idx];
	    results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = shift+u2_idx;
	  }
	  else break;
	}
      }
      for(intT j=0; j < used_idx; ++j) { wedges[used[shift+j]] = 0; }
    }
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);
  long total = 0;
  
  for(long i=0;i<nu;i++) total += butterflies[i];
  for(long i=0;i<stepSize;i++) total += results[i*eltsPerCacheLine];
  free(butterflies);
  free(results);
  cout << "num: " << total << "\n";
}

void CountOrigCompactParallel_WedgeAware(bipartiteCSR& GA, bool use_v, long* wedgesPrefixSum) {
  timer t1,t2;
  t1.start();
  cout << "Original Parallel Wedge-Aware" << endl;
  //cout << GA.nv << " " << GA.nu << " " << GA.numEdges << endl;
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = 1000; //tunable parameter

  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);

  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);
  
  long* results = newA(long,eltsPerCacheLine*stepSize); //one entry per cache line
  granular_for(i,0,stepSize,stepSize > 10000, {results[eltsPerCacheLine*i] = 0;});
  
  uintE* butterflies = newA(uintE,nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[i] = 0; });
  t1.reportTotal("preprocess");

  t2.start();
  
  //JS: try wedge-aware parallelism using wedge counts per vertex
  for(intT step = 0; step < (nu+stepSize-1)/stepSize; step++) {
      std::function<void(intT,intT)> recursive_lambda =
	[&]
	(intT start, intT end){
	if ((start == end-1) || (wedgesPrefixSum[end]-wedgesPrefixSum[start] < 1000)){ 
	  for (intT i = start; i < end; i++){
	    intT used_idx = 0;
	    intT shift = nu*(i-step*stepSize);
	    intT u_offset  = offsetsU[i];
	    intT u_deg = offsetsU[i+1]-u_offset;
	    for (intT j=0; j < u_deg; ++j ) {
	      intT v = edgesU[u_offset+j];
	      intT v_offset = offsetsV[v];
	      intT v_deg = offsetsV[v+1]-offsetsV[v];
	      for (intT k=0; k < v_deg; ++k) { 
		uintE u2_idx = edgesV[v_offset+k];
		if (u2_idx < i) {
		  //butterflies[i] += wedges[u2_idx];
		  //butterflies[u2_idx] += wedges[u2_idx];
		  results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
		  wedges[shift+u2_idx]++;
		  if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = shift+u2_idx;
		}
		else break;
	      }
	    }
	    for(intT j=0; j < used_idx; ++j) { wedges[used[shift+j]] = 0; }
	  }
	} else {
	  cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	  recursive_lambda(start + ((end-start)>>1), end);
	}
      }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,nu));
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);
  long total = 0;
  
  for(long i=0;i<nu;i++) total += butterflies[i];
  for(long i=0;i<stepSize;i++) total += results[i*eltsPerCacheLine];
  free(butterflies);
  free(results);
  cout << "num: " << total << "\n";
}

// Note: must be invoked with symmetricVertex
void Compute(bipartiteCSR& GA, commandLine P) {
  // Method type for counting + peeling
  long ty = P.getOptionLongValue("-t",0);
  long tp = P.getOptionLongValue("-tp",0);

  // # of max wedges
  long max_wedges = P.getOptionLongValue("-m",2577500000);

  //TODO wedgesPrefixSum not needed except in CountOrigCompactParallel_WedgeAware
  timer t1;
  t1.start();
  tuple<bool,long,long*> use_v_tuple = cmpWedgeCounts(GA);
  bool use_v = get<0>(use_v_tuple);
  long num_wedges = get<1>(use_v_tuple);
  t1.reportTotal("compute wedge counts + work prefix sum");
  
  
  //TODO seq code integrate w/count
  /*if (ty == 7) CountOrigCompactParallel(GA,use_v);
  else if (ty == 8) {
  	long* workPrefixSum = get<2>(use_v_tuple);
  	CountOrigCompactParallel_WedgeAware(GA,use_v,workPrefixSum);
  	free(workPrefixSum);
  	return;
  }
  else if (ty == 9) CountOrigCompact(GA,use_v);

  if (ty > 6) return;

  timer t;
  t.start();
  uintE* butterflies = Count(GA, use_v, num_wedges, max_wedges, ty);
  t.stop();

  if (ty==0) t.reportTotal("Sort:");
  else if (ty==1) t.reportTotal("SortCE:");
  else if (ty==2) t.reportTotal("Hash:");
  else if (ty==3) t.reportTotal("HashCE:");
  else if (ty==4) t.reportTotal("Hist:");
  else if (ty==6) t.reportTotal("HistCE:");

  const intT eltsPerCacheLine = 64/sizeof(long);
  long num_idxs = use_v ? GA.nu : GA.nv;
  long b = 0;
  for (long i=0; i < num_idxs; ++i) {b += butterflies[eltsPerCacheLine*i];}
  b = b / 2;
  cout << "number of butterflies: " << b << "\n";
  
  //uintE* butterflies2 = Count(GA, use_v, num_wedges, max_wedges, 2);
  //for (long i=0; i < num_idxs; ++i) { assertf(butterflies[eltsPerCacheLine*i] == butterflies2[eltsPerCacheLine*i], "%d, %d, %d", i, butterflies[eltsPerCacheLine*i], butterflies2[eltsPerCacheLine*i]); }

  timer t2;
  t2.start();
  auto cores = Peel(GA, use_v, butterflies, max_wedges, tp);
  t2.stop();
  if (tp ==0) t2.reportTotal("Hash Peel:");
  else if (tp==1) t2.reportTotal("Sort Peel:");
  else t2.reportTotal("Hist Peel:");

  uintE mc = 0;
  for (size_t i=0; i < num_idxs; i++) { mc = std::max(mc, cores[i]); }
  cout << "### Max core: " << mc << endl;

  free(butterflies);*/

 timer t3;
 
 auto eti = edgeToIdxs(GA, use_v);
 t3.start();
 uintE* ebutterflies = CountE(eti, GA, use_v, num_wedges, max_wedges, ty);
 t3.stop();
 if(ty==2) t3.reportTotal("E Hash:");
 else if (ty == 3) t3.reportTotal("E HashCE:");
 else if (ty == 0) t3.reportTotal("E Sort:");
 else if (ty==1) t3.reportTotal("E SortCE:");
 else if (ty==4) t3.reportTotal("E Hist:");

 const intT eltsPerCacheLine = 64/sizeof(long);
 long b=0;
 
 for (long i=0; i < GA.numEdges; ++i) {b += ebutterflies[eltsPerCacheLine*i];}
 cout << "number of edge butterflies: " << b/4 << "\n";


// timer t2;
// t2.start();
// auto cores = PeelE(G, use_v, ebutterflies, tp);
// t2.stop();
// if (tp ==0) t2.reportTotal("Hash Peel:");
// else if (tp==1) t2.reportTotal("Sort Peel:");
// else t2.reportTotal("Hist Peel:");

  free(eti);
  free(ebutterflies);
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," <inFile>");
  char* iFile = P.getArgument(0);
  //long rounds = P.getOptionLongValue("-rounds",3);

  bipartiteCSR G = readBipartite(iFile);

  Compute(G,P);
  // for(int r=0;r<rounds;r++) {
  //   startTime();
  //   Compute(G,P);
  //   nextTime("Running time");
  // }
  G.del();
}


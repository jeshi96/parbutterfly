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
#include "../../lib/histogram.h"
#include "../../lib/gbbs-histogram.h"

#include "butterfly_utils.h"

pair<intT, long> PeelESort(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, uintE* butterflies, bipartiteCSR& GA, 
  bool use_v, long max_wedges, intT curr_idx=0) {
  using X = tuple<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto ret = getIntersectWedges(eti, ite, current, ps, active_map, active.size(), GA, use_v, max_wedges, curr_idx);
  //ps.wedges_seq_tup_fil contains edge idx, val pairs; ret.first is size of this, ret.second is next idx
  // just need to hist up, update butterflies, return num updates

  pair<uintE*, long> b_freq_pair = getFreqs(ps.wedges_seq_tup_fil.A, ret.first, uintETupleLt(), uintETupleEq());
  uintE* b_freq_arr = b_freq_pair.first;
  ps.resize_update(b_freq_pair.second-1);
  uintE* update = ps.update_seq_int.A;

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    uintE num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(ps.wedges_seq_tup_fil.A[b_freq_arr[i-1]]), num_freq, uintETupleAdd());
    // These are our butterfly counts
    butterflies[get<0>(reduce)*eltsPerCacheLine] -= get<1>(reduce);
    update[i-1] = get<0>(reduce);
  }

  free(b_freq_arr);

  return make_pair(ret.second, b_freq_pair.second-1);
}

pair<intT, long> PeelEHist(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, uintE* butterflies, bipartiteCSR& GA, 
  bool use_v, long max_wedges, intT curr_idx=0) {
  using X = tuple<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto ret = getIntersectWedges(eti, ite, current, ps, active_map, active.size(), GA, use_v, max_wedges, curr_idx);
  //ps.wedges_seq_tup_fil contains edge idx, val pairs; ret.first is size of this, ret.second is next idx
  // just need to hist up, update butterflies, return num updates
  pbbsa::sequence<X> wedge_seq = pbbsa::sequence<X>(ps.wedges_seq_tup_fil.A, ret.first);
  // TODO fix these
  pbbsa::sequence<tuple<uintE, uintE>> tmp = pbbsa::sequence<tuple<uintE, uintE>>();
  pbbsa::sequence<tuple<uintE, uintE>> out = pbbsa::sequence<tuple<uintE, uintE>>();
  tuple<size_t, X*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(wedge_seq, GA.numEdges, getAdd<uintE,uintE>, getAddReduce<uintE,uintE>, tmp, out);
  X* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  ps.resize_update(butterflies_n);
  uintE* update = ps.update_seq_int.A;

  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] -= get<1>(butterflies_l[i]);
    update[i] = get<0>(butterflies_l[i]);
  }
  return make_pair(ret.second, butterflies_n);
}

long PeelEHash(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, uintE* butterflies, bipartiteCSR& GA, 
  bool use_v) {
  // Retrieve all seagulls
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  getIntersectWedgesHash(eti, ite, current, ps, active_map, active.size(), GA, use_v);
 
  _seq<pair<uintE,uintE>> update_seq = ps.update_hash.entries();
  long num_updates = update_seq.n;
  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;

  using T = pair<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);
  granular_for(i, 0, num_updates, (num_updates>1000), {
    T update_pair = update_seq.A[i];
    uintE idx = update_pair.first;
    butterflies[eltsPerCacheLine*idx] -= update_pair.second;
    update[i] = idx;
  });

  update_seq.del();

  return num_updates;
}

void PeelE_helper (uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, uintE* butterflies,
  bipartiteCSR& GA, bool use_v, long max_wedges, long type, array_imap<uintE>& D, buckets<array_imap<uintE>>& b, uintE k) {
  const long nu = use_v ? GA.nu : GA.nv;
  bool is_seq = (active.size() < 1000);
  /*if (type == 3) {
  	auto ret = PeelOrigParallel(ps, active, butterflies, update_dense, GA, use_v);
  	updateBuckets(D, b, k, butterflies, is_seq, nu, ret.first, ret.second);
  	free(ret.first);
  	return;
  }*/
  if (type == 0) {
    auto num_updates = PeelEHash(eti, ite, current, ps, active, butterflies, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ps.update_seq_int.A, num_updates);
    ps.clear();
    return;
  }

  pair<intT, long> ret;
  intT curr_idx = 0;
  while(curr_idx < active.size()) {
    if (type == 1) ret = PeelESort(eti, ite, current, ps, active, butterflies, GA, use_v, max_wedges, curr_idx);
  	else if (type == 2) ret = PeelEHist(eti, ite, current, ps, active, butterflies, GA, use_v, max_wedges, curr_idx);
    curr_idx = ret.first;
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ps.update_seq_int.A, ret.second);
    //ps.clear();
  }
}

array_imap<uintE> PeelE(uintE* eti, uintE* ite, bipartiteCSR& GA, bool use_v, uintE* butterflies, long max_wedges, long type=0, size_t num_buckets=128) {
  // Butterflies are assumed to be stored on U
  const long nu = use_v ? GA.nu : GA.nv;

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  
  using X = tuple<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  auto D = array_imap<uintE>(GA.numEdges, [&] (size_t i) { return butterflies[eltsPerCacheLine*i]; });

  auto b = make_buckets(GA.numEdges, D, increasing, num_buckets);

  //bool* update_dense = newA(bool, eltsPerCacheLine*GA.numEdges);
  PeelESpace ps = PeelESpace(type, GA.numEdges, MAX_STEP_SIZE);

  bool* current = newA(bool, GA.numEdges);
  parallel_for(size_t i=0; i < GA.numEdges; ++i) { current[i] = false; }

  size_t finished = 0;
  while (finished != GA.numEdges) {

    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();
    bool is_seq = (active.size() < 1000);

    //parallel_for(intT i=0; i < nu; ++i) { update_dense[eltsPerCacheLine*i] = false; }
    PeelE_helper(eti, ite, current, ps, active, butterflies, GA, use_v, max_wedges, type, D, b, k);

    active.del();
  }
  ps.del();
  free(current);
  
  return D;
}
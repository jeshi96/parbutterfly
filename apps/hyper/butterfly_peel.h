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

long PeelEHash(uintE* eti, uintE* ite, PeelSpace& ps, bool* used, vertexSubset& active, uintE* butterflies, bipartiteCSR& GA, bool use_v) {
  // Retrieve all seagulls
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  getIntersectWedgesHash(eti, ite, ps, used, active_map, active.size(), GA, use_v);
 
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

void PeelE_helper (uintE* eti, uintE* ite, PeelSpace& ps, bool* used, vertexSubset& active, uintE* butterflies,
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
    auto num_updates = PeelEHash(eti, ite, ps, used, active, butterflies, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ps.update_seq_int.A, num_updates);
    ps.clear();
    return;
  }
/*
  long num_wedges = countSeagulls(GA, use_v, active);
  pair<intT, long> ret;
//TODO remove update dense from all of these
  if (max_wedges >= num_wedges) {
    //else if (type == 1 && is_seq) ret = PeelSort_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else if (type == 1) ret = PeelESort(eti, ite, ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    //else if (type == 2 && is_seq) ret = PeelHist_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    else ret = PeelEHist(eti, ite, ps, active, butterflies, GA, use_v, num_wedges, max_wedges);
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ps.update_seq_int.A, ret.second);
    ps.clear();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < active.size()) {
    //else if (type == 1 && is_seq) ret = PeelSort_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else if (type == 1) ret = PeelESort(eti, ite, ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
  	//else if (type == 2 && is_seq) ret = PeelHist_seq(ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    else ret = PeelEHist(eti, ite, ps, active, butterflies, GA, use_v, num_wedges, max_wedges, curr_idx);
    curr_idx = ret.first;
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ps.update_seq_int.A, ret.second);
    ps.clear();
  }*/
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
  PeelSpace ps = PeelSpace(type, GA.numEdges, MAX_STEP_SIZE);

  bool* used = newA(bool, GA.numEdges);
  parallel_for(size_t i=0; i < GA.numEdges; ++i) { used[i] = false; }

  size_t finished = 0;
  while (finished != GA.numEdges) {

    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();
    bool is_seq = (active.size() < 1000);

    //parallel_for(intT i=0; i < nu; ++i) { update_dense[eltsPerCacheLine*i] = false; }
    PeelE_helper(eti, ite, ps, used, active, butterflies, GA, use_v, max_wedges, type, D, b, k);

    active.del();
  }
  ps.del();
  free(used);
  
  return D;
}
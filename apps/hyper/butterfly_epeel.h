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
//#include "../../lib/gbbs-histogram.h"

#include "butterfly_putils.h"
#include "butterfly_utils.h"

pair<intT, long> PeelESort(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, 
			   bool use_v, long max_wedges, intT curr_idx=0) {
  using X = tuple<uintE,long>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto ret = getIntersectWedges(eti, ite, current, ps, active_map, active.size(), GA, use_v, max_wedges, curr_idx);
  //ps.wedges_seq_tup_fil contains edge idx, val pairs; ret.first is size of this, ret.second is next idx
  // just need to hist up, update butterflies, return num updates

  auto b_freq_pair = getFreqs<long>(ps.wedges_seq_tup_fil.A, ret.first, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;
  ps.resize_update(b_freq_pair.second-1);
  uintE* update = ps.update_seq_int.A;

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    long num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(ps.wedges_seq_tup_fil.A[b_freq_arr[i-1]]), num_freq, tupleAdd<uintE, long>());
    // These are our butterfly counts
    butterflies[get<0>(reduce)*eltsPerCacheLine] -= get<1>(reduce);
    update[i-1] = get<0>(reduce);
  }

  free(b_freq_arr);

  return make_pair(ret.second, b_freq_pair.second-1);
}

pair<intT, long> PeelEHist(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, 
			   bool use_v, long max_wedges, intT curr_idx=0) {
  using X = tuple<uintE,long>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  auto ret = getIntersectWedges(eti, ite, current, ps, active_map, active.size(), GA, use_v, max_wedges, curr_idx);
  //ps.wedges_seq_tup_fil contains edge idx, val pairs; ret.first is size of this, ret.second is next idx
  // just need to hist up, update butterflies, return num updates
  pbbsa::sequence<X> wedge_seq = pbbsa::sequence<X>(ps.wedges_seq_tup_fil.A, ret.first);
  // TODO fix these
  pbbsa::sequence<tuple<uintE, long>> tmp = pbbsa::sequence<tuple<uintE, long>>();
  pbbsa::sequence<tuple<uintE, long>> out = pbbsa::sequence<tuple<uintE, long>>();
  tuple<size_t, X*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(wedge_seq, GA.numEdges, getAdd<uintE,long>, getAddReduce<uintE,long>, tmp, out);
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

long PeelEHash(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, long* butterflies, bipartiteCSR& GA, 
	       bool use_v) {
  // Retrieve all seagulls
  auto active_map = make_in_imap<uintT>(active.size(), [&] (size_t i) { return active.vtx(i); });
  getIntersectWedgesHash(eti, ite, current, ps, active_map, active.size(), GA, use_v);
 
  auto update_seq = ps.update_hash.entries();
  long num_updates = update_seq.n;
  ps.resize_update(num_updates);
  uintE* update = ps.update_seq_int.A;

  const intT eltsPerCacheLine = 64/sizeof(long);
  granular_for(i, 0, num_updates, (num_updates>1000), {
      auto update_pair = update_seq.A[i];
      uintE idx = update_pair.first;
      butterflies[eltsPerCacheLine*idx] -= update_pair.second;
      update[i] = idx;
    });

  update_seq.del();

  return num_updates;
}

pair<uintE*, long> PeelEOrigParallel(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, long* butterflies, bool* update_dense, bipartiteCSR& GA, bool use_v) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = active.size() < MAX_STEP_SIZE ? active.size() : MAX_STEP_SIZE; //tunable parameter

  auto wedges_seq = ps.wedges_seq_int;
  auto used_seq = ps.used_seq_int;

  uintE* wedges = wedges_seq.A;
  uintE* used = used_seq.A;

  //granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });

  const intT eltsPerCacheLine = 64/sizeof(long);
  granular_for(i,0,GA.numEdges,GA.numEdges > 10000, { update_dense[eltsPerCacheLine*i] = false; });

  parallel_for(intT i=0; i < active.size(); ++i){ current[active.vtx(i)] = 1; }

  for(intT step = 0; step < (active.size()+stepSize-1)/stepSize; step++) {
    parallel_for_1(intT i=step*stepSize; i < min((step+1)*stepSize,active.size()); ++i){//parallel_for_1
      intT used_idx = 0;
      intT shift = nu*(i-step*stepSize);
      intT idx_vu = active.vtx(i);
  
      uintE u = edgesV[idx_vu];
      intT u_offset = offsetsU[u];
      intT u_deg = offsetsU[u+1] - u_offset;

      uintE v = edgesU[ite[idx_vu]];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // we want to go through all neighbors of v
      // intersect each N(u') with N(u)
      // get all u, v, u2
      // intersection of N(u) and N(u2) - 1 is the num butterflies to subtract from v, u2
      for (intT k=0; k < v_deg; ++k) { 
	uintE u2 = edgesV[v_offset + k];
	if (u2 != u && u2 != UINT_E_MAX && (!current[v_offset+k] || idx_vu < v_offset + k)) {//u2 != UINT_E_MAX - 1

	  intT u2_offset = offsetsU[u2];
	  intT u2_deg = offsetsU[u2+1] - u2_offset;

	  //bool* same = newA(bool, u2_deg);
	  intT int_offset = 0;
	  //parallel_for(intT j = 0; j < u2_deg; ++j) { same[j] = 0; }
	  for(intT j = 0; j < u2_deg; ++j) {
	    uintE v2 = edgesU[u2_offset + j];
	    uintE idx_v2u2 = eti[u2_offset + j];
	    if(v2 != UINT_E_MAX && v2 != v && (!current[idx_v2u2] || idx_v2u2 >= idx_vu)) {
	      while(int_offset < u_deg && (edgesU[u_offset + int_offset] == UINT_E_MAX || edgesU[u_offset + int_offset] < v2)) { int_offset++; }
	      if (int_offset < u_deg && edgesU[u_offset+int_offset] == v2 && (!current[eti[u_offset + int_offset]] || idx_vu < eti[u_offset + int_offset])) {
		//same[j] = 1;
		wedges[shift+k]++;
		writeAdd(&butterflies[eltsPerCacheLine*eti[u2_offset + j]], (long) -1);
		writeAdd(&butterflies[eltsPerCacheLine*eti[u_offset + int_offset]], (long) -1);
		if(!update_dense[eltsPerCacheLine*eti[u2_offset + j]]) CAS(&update_dense[eltsPerCacheLine*eti[u2_offset + j]],false,true);
		if(!update_dense[eltsPerCacheLine*eti[u_offset + int_offset]]) CAS(&update_dense[eltsPerCacheLine*eti[u_offset + int_offset]],false,true);
		if (wedges[shift+k] == 1) {
		  used[shift+used_idx++] = k;
		  if(!update_dense[eltsPerCacheLine*(v_offset+k)]) CAS(&update_dense[eltsPerCacheLine*(v_offset+k)],false,true);
		}
	      }
	      else if(int_offset >= u_deg) break;
	    }
	  }
	  //long num_same = sequence::sum(same, u2_deg);
	  //if (num_same > 0) { 
	  //  writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k), -1*num_same);
	  //}
	  //free(same);
	}
      }

      granular_for(j,0,used_idx,used_idx > 10000, { 
	  intT k = used[shift+j];
	  writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)],(long) -1* ((long)wedges[shift+k]));
	  wedges[shift+k] = 0;
	});
    }
  }

  parallel_for(intT i=0; i < active.size(); ++i){
    edgesV[active.vtx(i)] = UINT_E_MAX; edgesU[ite[active.vtx(i)]] = UINT_E_MAX;
    current[active.vtx(i)] = 0;
  }
  auto f = [&] (size_t i) { return update_dense[eltsPerCacheLine*i]; };
  auto f_in = make_in_imap<bool>(GA.numEdges, f);
  auto out = pbbs::pack_index<uintE>(f_in);
  out.alloc = false;
  return make_pair(out.s, out.size());
}

void PeelE_helper (uintE* eti, uintE* ite, bool* current, PeelESpace& ps, vertexSubset& active, long* butterflies,  bool* update_dense,
		   bipartiteCSR& GA, bool use_v, long max_wedges, long type, array_imap<long>& D, buckets<array_imap<long>>& b, long k) {
  const long nu = use_v ? GA.nu : GA.nv;
  bool is_seq = (active.size() < 1000);
  if (type == 3) {
    auto ret = PeelEOrigParallel(eti, ite, current, ps, active, butterflies, update_dense, GA, use_v);
    updateBuckets(D, b, k, butterflies, is_seq, GA.numEdges, ret.first, ret.second);
    free(ret.first);
    return;
  }
  else if (type == 0) {
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

array_imap<long> PeelE(uintE* eti, uintE* ite, bipartiteCSR& GA, bool use_v, long* butterflies, long max_wedges, long type=0, size_t num_buckets=128) {
  // Butterflies are assumed to be stored on U
  const long nu = use_v ? GA.nu : GA.nv;

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  
  using X = tuple<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  auto D = array_imap<long>(GA.numEdges, [&] (size_t i) { return butterflies[eltsPerCacheLine*i]; });

  auto b = make_buckets(GA.numEdges, D, increasing, num_buckets);

  bool* update_dense = nullptr;
  if (type == 3) update_dense = newA(bool, eltsPerCacheLine*GA.numEdges);
  PeelESpace ps = PeelESpace(type, GA.numEdges, MAX_STEP_SIZE, nu);

  bool* current = newA(bool, GA.numEdges);
  parallel_for(size_t i=0; i < GA.numEdges; ++i) { current[i] = false; }

  size_t finished = 0;
  while (finished != GA.numEdges) {

    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    long k = bkt.id;
    finished += active.size();
    bool is_seq = (active.size() < 1000);

    PeelE_helper(eti, ite, current, ps, active, butterflies, update_dense, GA, use_v, max_wedges, type, D, b, k);

    active.del();
  }
  ps.del();
  b.del();
  if (type == 3) free(update_dense);
  free(current);
  
  return D;
}

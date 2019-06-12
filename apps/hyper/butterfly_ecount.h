#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "math.h"
#include "graph.h"
#include "vertex.h"
#include "sequence.h"
#include "sparseSet.h"
#include "sampleSort.h"
#include "../../lib/histogram.h"
#include "../../lib/gbbs-histogram.h"

#include "butterfly_utils.h"

using namespace std;

timer rehashWedgesTimer, getWedgesTimer, retrieveCountsTimer;

intT CountEHistCE(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
		intT curr_idx=0) {
  using X = tuple<long,uintE>;
  using T = tuple<uintE, long>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  pair<long, intT> wedges_pair  = getWedges<long>(cs.wedges_seq_int, GA, UWedgeIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);
  intT next_idx = wedges_pair.second;
  long num_wedges_list = wedges_pair.first;

  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(cs.wedges_seq_int.A, wedges_pair.first);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  cs.wedges_hash.resize(wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

  if (cs.butterflies_seq_intt.n < 2*num_wedges_list) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(T, 2*num_wedges_list);
    cs.butterflies_seq_intt.n = 2*num_wedges_list;
  }

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    long idx = 0;
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    for(intT j=0;j<u_deg;j++) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1] - v_offset;
	if (v > i) { 
	  for (intT k=0; k < v_deg; ++k) { 
	  uintE u2 = GA.edges[v_offset+k] >> 1;
	  if (u2 > i) {
	    long to_find = ((((long) i) *GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	    long num_butterflies = cs.wedges_hash.find(to_find).second;
	    if (num_butterflies > 1) {
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

  pbbsa::sequence<T> wedge_freqs_i_seq = pbbsa::sequence<T>(cs.butterflies_seq_intt.A,2*num_wedges_list);
  tuple<size_t, T*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(wedge_freqs_i_seq,GA.numEdges,getAdd<uintE,long>, getAddReduce<uintE,long>, cs.tmp_uint, cs.out_uint);
  T* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  parallel_for (long i=0; i < butterflies_n; ++i) {
    if (get<0>(butterflies_l[i]) != UINT_E_MAX) butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  return next_idx;
}

intT CountEHistCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  using X = tuple<long,uintE>;
  using T = tuple<uintE, long>;
  auto cons = UVertexPairIntCons(nu);

  pair<long, intT> wedges_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_list = wedges_pair.first;
  intT next_idx = wedges_pair.second;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_list);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  cs.wedges_hash.resize(wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }
  if (cs.butterflies_seq_intt.n < 2*num_wedges_list) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(T, 2*num_wedges_list);
    cs.butterflies_seq_intt.n = 2*num_wedges_list;
  }

  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    long idx = 0;
    //granular_for (j, 0, u_deg, (u_deg > 1000), {
    for(intT j=0;j<u_deg;j++) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1] - v_offset;
	// Find all seagulls with center v and endpoint u
	for (intT k=0; k < v_deg; ++k) {
	  uintE u2 = edgesV[v_offset+k];
	  if (u2 < i) {
	    // TODO should take out -1, put in if to store only if > 0
	    long num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
      cs.butterflies_seq_intt.A[2*(wedge_idx+idx)] = make_tuple(eti[u_offset + j], num_butterflies);
      cs.butterflies_seq_intt.A[2*(wedge_idx+idx)+1] = make_tuple((v_offset + k), num_butterflies);
      ++idx;
	  }
	  else break;
	}
    }
  }

  pbbsa::sequence<T> wedge_freqs_i_seq = pbbsa::sequence<T>(cs.butterflies_seq_intt.A,2*num_wedges_list);
  tuple<size_t, T*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(wedge_freqs_i_seq,GA.numEdges,getAdd<uintE,long>, getAddReduce<uintE,long>, cs.tmp_uint, cs.out_uint);
  T* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  return next_idx;
}

intT CountEHist(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long* butterflies_u, long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  using X = tuple<long,uintE>;
  auto cons = UVertexPairIntCons(nu);
  //getWedgesTimer.start();
  pair<long, intT> wedges_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_list = wedges_pair.first;
  intT next_idx = wedges_pair.second;

  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_list);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);
  //getWedgesTimer.stop();
  //rehashWedgesTimer.start();
  //TODO save this space too
  cs.wedges_hash.resize(wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }
  //rehashWedgesTimer.stop();
  //retrieveCountsTimer.start();

  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    //granular_for (j, 0, u_deg, (u_deg > 1000), {
    parallel_for(intT j=0;j<u_deg;j++) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1] - v_offset;
	// Find all seagulls with center v and endpoint u
	for (intT k=0; k < v_deg; ++k) {
	  uintE u2 = edgesV[v_offset+k];
	  if (u2 < i) {
	    // TODO should take out -1, put in if to store only if > 0
	    long num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
      writeAdd(&butterflies_u[eltsPerCacheLine*(u_offset + j)],num_butterflies); 
	    writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)],num_butterflies); 
	  }
	  else break;
	}
    }
  }
  //retrieveCountsTimer.stop();
  return next_idx;
}

intT CountEHist(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
		intT curr_idx=0) {
  using X = tuple<long,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  pair<long, intT> wedges_pair  = getWedges<long>(cs.wedges_seq_int, GA, UWedgeIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);
  intT next_idx = wedges_pair.second;

  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(cs.wedges_seq_int.A, wedges_pair.first);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  cs.wedges_hash.resize(wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    //granular_for (j, 0, u_deg, (u_deg > 1000), {
    parallel_for(intT j=0;j<u_deg;j++) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1] - v_offset;
	if (v > i) { 
	  for (intT k=0; k < v_deg; ++k) { 
	  uintE u2 = GA.edges[v_offset+k] >> 1;
	  if (u2 > i) {
	    long to_find = ((((long) i) *GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	    long num_butterflies = cs.wedges_hash.find(to_find).second;
	    if (num_butterflies > 1) {
	      writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], num_butterflies - 1);
	      writeAdd(&butterflies[eltsPerCacheLine*(u_offset+j)], num_butterflies - 1);
	    } // store on edges -- note must sum eti and ite versions in the end
	  }
	  else break;
	  }
	}
    }//);
  }

  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

intT CountESortCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  //getWedgesTimer.start();
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, use_v, UWedgeCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;
  //getWedgesTimer.stop();
  //rehashWedgesTimer.start();
  using X = tuple<uintE,long>;
  if (cs.butterflies_seq_intt.n < 2*num_wedges_f) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*num_wedges_f);
    cs.butterflies_seq_intt.n = 2*num_wedges_f;
  }

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num = freq_arr[i+1] - freq_arr[i] - 1;
    granular_for (j, freq_arr[i], freq_arr[i+1], (freq_arr[i+1]-freq_arr[i] > 1000), {
	cs.butterflies_seq_intt.A[(long)2*j] = make_tuple(eti[offsetsU[wedges[j].v1] + wedges[j].j], num);
	cs.butterflies_seq_intt.A[(long)2*j+1] = make_tuple(offsetsV[wedges[j].u] + wedges[j].k, num);
    });
  }
  //rehashWedgesTimer.stop();
  //retrieveCountsTimer.start();

  free(freq_arr);

  // now, we need to collate by our indices
  auto b_freq_pair = getFreqs(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;
  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    long num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]), num_freq, tupleAdd<uintE,long>());
    // These are our butterfly counts
    butterflies[eltsPerCacheLine*get<0>(reduce)] += get<1>(reduce);
  }

  free(b_freq_arr);

  //retrieveCountsTimer.stop();
  return wedges_pair.second;
}

intT CountESort(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long* butterflies_u, long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;

  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, use_v, UWedgeCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;

  //getWedgesTimer.start();
  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;
  //getWedgesTimer.stop();

  //rehashWedgesTimer.start();
  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num = freq_arr[i+1] - freq_arr[i] - 1;
    if (num > 0) {
    granular_for (j, freq_arr[i], freq_arr[i+1], (freq_arr[i+1]-freq_arr[i] > 1000), {
	writeAdd(&butterflies_u[eltsPerCacheLine*(offsetsU[wedges[j].v1] + wedges[j].j)], num);
	writeAdd(&butterflies[eltsPerCacheLine*(offsetsV[wedges[j].u] + wedges[j].k)], num);
    });
    }
  }
  //rehashWedgesTimer.stop();
  free(freq_arr);
  return wedges_pair.second;
}

intT CountESortCE(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges,
		  long* wedge_idxs, intT curr_idx=0) {
  const intT eltsPerCacheLine = 64/sizeof(long);
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx, num_wedges,
						    wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  using X = tuple<uintE, long>;
  if (cs.butterflies_seq_intt.n < 2*num_wedges_f) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*num_wedges_f);
    cs.butterflies_seq_intt.n = 2*num_wedges_f;
  }

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];
    if (num_butterflies > 1){
      parallel_for(long j=freq_arr[i]; j<freq_arr[i+1]; ++j) { //JS: test granular_for
	cs.butterflies_seq_intt.A[2*j] = make_tuple((GA.offsets[wedges[j].u]+wedges[j].k), num_butterflies - 1);
	cs.butterflies_seq_intt.A[2*j+1] = make_tuple((GA.offsets[wedges[j].v1]+wedges[j].j), num_butterflies - 1);
      }
    }
    else {
      parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) { //JS: test granular_for
        cs.butterflies_seq_intt.A[2*j] = make_tuple(UINT_E_MAX, 0);
        cs.butterflies_seq_intt.A[2*j+1] = make_tuple(UINT_E_MAX, 0);
      }
    }
  }

  free(freq_arr);

  // now, we need to collate by our indices
  auto b_freq_pair = getFreqs<long>(cs.butterflies_seq_intt.A, 2*num_wedges_f, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;

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

intT CountESort(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx, num_wedges,
						    wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());

  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for (long i = 0; i < freq_pair.second-1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];

    if (num_butterflies > 1) {
      granular_for (j, freq_pair.first[i], freq_pair.first[i+1], (freq_pair.first[i+1]-freq_pair.first[i] > 1000), {
	  // GA.offsets[wedges[j].u]+wedges[j].k --> this is v_offset + k
        // GA.offsets[wedges[j].v1]+wedges[j].j --> this is u_offset + j
	  writeAdd(&butterflies[eltsPerCacheLine*(GA.offsets[wedges[j].u]+wedges[j].k)], num_butterflies - 1);
	  writeAdd(&butterflies[eltsPerCacheLine*(GA.offsets[wedges[j].v1]+wedges[j].j)], num_butterflies - 1);
      });
    }
  }

  free(freq_pair.first);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

intT CountEHash(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
		intT curr_idx=0) {
  const intT eltsPerCacheLine = 64/sizeof(long);

  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1] - v_offset;
	if (v > i) {
	  for (intT k=0; k < v_deg; ++k) { 
	    uintE u2 = GA.edges[v_offset+k] >> 1;
	    if (u2 > i) {
	      long to_find = ((((long) i) *GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	      long num_butterflies = cs.wedges_hash.find(to_find).second;
	      if (num_butterflies > 1) {
		writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], num_butterflies - 1);
		writeAdd(&butterflies[eltsPerCacheLine*(u_offset+j)], num_butterflies - 1);
	      } // store on edges -- note must sum eti and ite versions in the end
	    }
	    else break;
	  }
	}
    }
  }
  return next_idx;
}

intT CountEHashCE(CountESpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
		  intT curr_idx=0) {
  using T = pair<uintE, long>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) { //JS: test granular_for
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1] - v_offset;
	if (v > i) { 
	  for (intT k=0; k < v_deg; ++k) { 
	    uintE u2 = GA.edges[v_offset+k] >> 1;
	    if (u2 > i) {
	      long to_find = ((((long) i) *GA.n + (long) u2) << 1) + (GA.edges[v_offset+k] & 0b1);
	      long num_butterflies = cs.wedges_hash.find(to_find).second;
	      if (num_butterflies > 1) {
		cs.butterflies_hash.insert(T((v_offset+k), num_butterflies-1));
		cs.butterflies_hash.insert(T((u_offset+j), num_butterflies-1));
	      } // store on edges -- note must sum eti and ite versions in the end
	    }
	    else break;
	  }
	}
    }
  }
  size_t num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.wedges_seq_intp);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterfly_pair = cs.wedges_seq_intp.A[i];
    butterflies[eltsPerCacheLine*butterfly_pair.first] += butterfly_pair.second;
  }
  return next_idx;
}

intT CountEHash(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long* butterflies_u, long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  //getWedgesTimer.start();
  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
  //getWedgesTimer.stop();

  //rehashWedgesTimer.start();
  const intT eltsPerCacheLine = 64/sizeof(long);
  // TODO modularize w/hashce
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    //granular_for (j, 0, u_deg, (u_deg > 1000), {
    parallel_for(intT j=0;j<u_deg;j++) { //JS: test granular_for
	uintE v = edgesU[u_offset + j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1] - v_offset;
	// Find all seagulls with center v and endpoint u
	for (intT k=0; k < v_deg; ++k) {
	  uintE u2 = edgesV[v_offset + k];
	  if (u2 < i) {
	    long num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
	    writeAdd(&butterflies_u[eltsPerCacheLine*(u_offset + j)],num_butterflies); 
	    writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)],num_butterflies); 
	  }
	  else break;
	}
    }
  }
  //rehashWedgesTimer.stop();
  
  return next_idx;
}

intT CountEHashCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  //getWedgesTimer.start();
  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
  //getWedgesTimer.stop();
  //rehashWedgesTimer.start();
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    //granular_for (j, 0, u_deg, (u_deg > 1000), {
    parallel_for(intT j=0;j<u_deg;j++) { //JS: test granular_for
	uintE v = edgesU[u_offset + j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1] - v_offset;
	// Find all seagulls with center v and endpoint u
	for (intT k=0; k < v_deg; ++k) {
	  uintE u2 = edgesV[v_offset + k];
	  if (u2 < i) {
	    long num_butterflies = cs.wedges_hash.find(((long) i)*nu + (long)u2).second - 1;
	    if (num_butterflies > 0) {
	    cs.butterflies_hash.insert(make_pair(eti[u_offset + j], num_butterflies));
	    cs.butterflies_hash.insert(make_pair(v_offset + k, num_butterflies));
	    }
	  }
	  else break;
	}
    }
  }
  //rehashWedgesTimer.stop();
  //retrieveCountsTimer.start();
  size_t num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.wedges_seq_intp);

  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterflies_pair = cs.wedges_seq_intp.A[i];
    butterflies[eltsPerCacheLine*butterflies_pair.first] += butterflies_pair.second;
  }
  //retrieveCountsTimer.stop();
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

void CountEOrigCompactParallel(uintE* eti, long* butterflies, long* butterflies_u, bipartiteCSR& GA, bool use_v, long max_array_size) {
  timer t1,t2,t3;
  //cout << "Original Parallel for Edges" << endl;
  t1.start();
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = min<long>(getWorkers() * 15, max_array_size/nu); //15 tunable parameter
  //cout << stepSize << " " << getWorkers() * 60 << " " << max_array_size / nu << endl;
  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);
  //cout << nu*stepSize << endl;
  t1.reportTotal("preprocess (malloc)");
  t3.start();
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);
  
  t3.reportTotal("preprocess (initialize)");

  t2.start();

  for(long step = 0; step < (nu+stepSize-1)/stepSize; step++) {
    parallel_for_1(long i=step*stepSize; i < min((step+1)*stepSize,nu); ++i){
      intT used_idx = 0;
      long shift = nu*(i-step*stepSize);
      intT u_offset  = offsetsU[i];
      intT u_deg = offsetsU[i+1]-u_offset;
      for (long j=0; j < u_deg; ++j ) { //JS: test granular_for
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1]-offsetsV[v];
	for (long k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx < i) {
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	  }
	  else break;
	}
      }
      for(long j=0; j < u_deg; ++j) { //JS: test granular_for
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

      for(long j=0; j < used_idx; ++j) { wedges[used[shift+j]+shift] = 0; }
    }
  }
  t2.reportTotal("counting (main loop)");
  
  free(wedges);
  free(used);
}

void CountEWorkEfficientParallel(graphCSR& GA, long* butterflies, long max_array_size) {
  timer t1,t2,t3;
  //cout << "Original Work-efficient Parallel for Edges" << endl;
  t1.start();

  long stepSize = min<long>(getWorkers() * 7.5, max_array_size/GA.n); //15 tunable parameter
  uintE* wedges = newA(uintE, GA.n*stepSize);
  uintE* used = newA(uintE, GA.n*stepSize);
  t1.reportTotal("preprocess (malloc)");
  t3.start();
  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);

  t3.reportTotal("preprocess (initialize)");

  t2.start();

  for(long step = 0; step < (GA.n+stepSize-1)/stepSize; step++) {
    parallel_for_1(long i=step*stepSize; i < min((step+1)*stepSize,GA.n); ++i){
      intT used_idx = 0;
      long shift = GA.n*(i-step*stepSize);
      intT u_offset  = GA.offsets[i];
      intT u_deg = GA.offsets[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) { //JS: test granular_for
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
	if (v <= i) break; 
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
	  if (u2_idx > i) {
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = GA.edges[v_offset+k];
	  }
	  else break;
	}
      }

      for (long j=0; j < u_deg; ++j ) { //JS: test granular_for
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
	if (v <= i) break; 
	for (long k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
	  if (u2_idx > i) { //TODO combine into one graph
	    if (wedges[shift+u2_idx] > 1) {
	      writeAdd(&butterflies[eltsPerCacheLine*(v_offset+k)], (long)(wedges[shift+u2_idx]-1));
	      writeAdd(&butterflies[eltsPerCacheLine*(u_offset+j)], (long)(wedges[shift+u2_idx]-1));
	    }
	  }
	  else break;
	}
      }

      for(long j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j] >> 1;
        wedges[shift+u2_idx] = 0;
      }
    }
  }
  t2.reportTotal("counting (main loop)");
  
  free(wedges);
  free(used);
}

long* CountERank(uintE* eti, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, long type,
		 uintE* ranks, uintE* rankV, uintE* rankU) {
  timer t_rank;
  t_rank.start();
  auto g_pair = rankGraphEdges(GA, use_v, ranks, rankV, rankU);
  auto g = g_pair.first;
  auto edge_converter = g_pair.second;
  free(ranks); free(rankU); free(rankV); 
  t_rank.reportTotal("ranking");

  timer t_time, t_time2;
  t_time.start();

  const intT eltsPerCacheLine = 64/sizeof(long);
  long numEdges = g.offsets[g.n]; //TODO check is this g.numEdges*2
  long* rank_butterflies = newA(long,eltsPerCacheLine*numEdges);
  t_time.reportTotal("preprocess (malloc)");

  t_time2.start();
  granular_for(i,0,numEdges,numEdges > 10000, { rank_butterflies[eltsPerCacheLine*i] = 0; });

  long* wedge_idxs = (type == 11) ? nullptr : countWedgesScan(g);
  CountESpace cs = CountESpace(type, GA.numEdges, true);

  if (type == 11) CountEWorkEfficientParallel(g, rank_butterflies, max_array_size);
  else {
    if (max_wedges >= num_wedges) {
      if (type == 0) CountESort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == 1) CountESortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == 2) CountEHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == 3) CountEHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == 4) CountEHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
      else if (type == 6) CountEHistCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
    }
    else {
      intT curr_idx = 0;
      while(curr_idx < g.n) {
	if (type ==0) curr_idx = CountESort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type ==1) curr_idx = CountESortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type ==2) curr_idx = CountEHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == 3) curr_idx = CountEHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	else if (type == 4) curr_idx = CountEHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
  else if (type == 6) curr_idx = CountEHistCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
	cs.clear();
      }
    }
  }

  cs.del();
  if (type != 11) free(wedge_idxs);

  //uintE* rank_butterflies2 = newA(uintE,eltsPerCacheLine*g.n);
  //granular_for(i,0,g.n,g.n > 10000, { rank_butterflies2[eltsPerCacheLine*i] = 0; });
  //CountWorkEfficientParallel(g, rank_butterflies2);
  //for(long i=0; i < g.n; ++i) {assert(rank_butterflies2[eltsPerCacheLine*i] == rank_butterflies[eltsPerCacheLine*i]);}

  if (type != 11) t_time2.reportTotal("counting");

  timer t_convert;
  t_convert.start();

  long* butterflies = newA(long, eltsPerCacheLine*GA.numEdges);
  granular_for(i,0,GA.numEdges,GA.numEdges > 10000, { butterflies[eltsPerCacheLine * i] = 0; });

  parallel_for(intT i=0; i < g.n; ++i) {
    intT v_offset = g.offsets[i];
    intT v_deg = g.offsets[i+1] - v_offset;
    granular_for(j, 0, v_deg, v_deg > 1000, {
	uintE u = g.edges[v_offset + j] >> 1;
	// can get rid of writeadd by using 2 butterfly arrays and then summing them
	if (g.edges[v_offset + j] & 0b1) writeAdd(&butterflies[eltsPerCacheLine*(get<1>(edge_converter[v_offset+j]))], rank_butterflies[eltsPerCacheLine*(v_offset+j)]);
	else writeAdd(&butterflies[eltsPerCacheLine*eti[(get<1>(edge_converter[v_offset+j]))]], rank_butterflies[eltsPerCacheLine*(v_offset+j)]);
      });
  }
  free(rank_butterflies);
  free(edge_converter);
  g.del();

  t_convert.reportTotal("convert");

  return butterflies;
}

long* CountE(uintE* eti, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, long type=0, long tw=0) {
  const long nu = use_v ? GA.nu : GA.nv;

  const intT eltsPerCacheLine = 64/sizeof(long);
  if (tw !=0) {
    timer t_rank;
    t_rank.start();
    auto rank_tup = tw == 1 ? getCoCoreRanks(GA) : getApproxCoCoreRanks(GA);

    long num_ccwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
					       rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup), get<1>(rank_tup)));
    num_ccwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
					   rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup)));
  
    t_rank.reportTotal("ranking");

    //cout << "Side wedges: " << num_wedges << "\n";
    //cout << "Co Core wedges: " << num_ccwedges << "\n";

    if (num_ccwedges < num_wedges + 1000 || tw == 1 || tw == 2) return CountERank(eti, GA, use_v, num_ccwedges, max_wedges, max_array_size, type, get<0>(rank_tup), get<1>(rank_tup), get<2>(rank_tup));
    free(get<0>(rank_tup)); free(get<1>(rank_tup)); free(get<2>(rank_tup));
  }

  timer t_malloc; t_malloc.start();

  long* butterflies = newA(long, eltsPerCacheLine*GA.numEdges);
  parallel_for(long i=0; i < GA.numEdges; ++i){
    butterflies[eltsPerCacheLine*i] = 0;
  }

  long* butterflies_u = (type == 3 || type == 1 || type == 6) ? nullptr : newA(long, eltsPerCacheLine*GA.numEdges);
  if (type != 3 && type != 1 && type != 6) {
  parallel_for(long i=0; i < GA.numEdges; ++i){
    butterflies_u[eltsPerCacheLine*i] = 0;
  }
  }

  t_malloc.reportTotal("preprocess (malloc)"); timer t_time;
  t_time.start();

  long* wedge_idxs = type == 5 ? nullptr : countWedgesScan(GA, use_v, true);
  CountESpace cs = CountESpace(type, GA.numEdges, false);

  if (type == 5) CountEOrigCompactParallel(eti, butterflies, butterflies_u, GA, use_v, max_array_size);
  else {
    // TODO clean up utils, do overflow stuff (double check youtube needs it)
    // TODO check correctness against each other
    // TODO why is hash so slow
    if (max_wedges >= num_wedges) {
      if (type == 0) CountESort(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti);
      else if (type == 1) CountESortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
      else if (type == 2) CountEHash(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti);
      else if (type == 3) CountEHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
      else if (type == 4) CountEHist(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti);
      else if (type == 6) CountEHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    }
    else {
      intT curr_idx = 0;
      while(curr_idx < nu) {
	if (type == 0) curr_idx = CountESort(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == 1) curr_idx = CountESortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == 2) curr_idx = CountEHash(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == 3) curr_idx = CountEHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
	else if (type == 4) curr_idx = CountEHist(cs, GA, use_v, num_wedges, butterflies, butterflies_u, max_wedges, wedge_idxs, eti, curr_idx);
  else if (type == 6) curr_idx = CountEHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
	cs.clear();
      }
    }
  }
  cs.del();

  if (type != 5) free(wedge_idxs);
  if (type != 5) t_time.reportTotal("counting");

  if (type !=3 && type != 1 && type != 6) {
  timer t_convert;
  t_convert.start();
  parallel_for(long i=0; i < GA.numEdges; ++i) {
    butterflies[eltsPerCacheLine*eti[i]] += butterflies_u[eltsPerCacheLine*i];
  }
  free(butterflies_u);
  t_convert.reportTotal("convert");
  }

  return butterflies;
}

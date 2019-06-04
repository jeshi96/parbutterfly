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

intT CountEHist(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  using X = tuple<long,uintE>;
  auto cons = UVertexPairIntCons(nu);
getWedgesTimer.start();
  pair<long, intT> wedges_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_list = wedges_pair.first;
  intT next_idx = wedges_pair.second;

  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_list);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);
getWedgesTimer.stop();
rehashWedgesTimer.start();
  //TODO save this space too
  cs.wedges_hash.resize(wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }
rehashWedgesTimer.stop();
retrieveCountsTimer.start();

const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    //granular_for (j, 0, u_deg, (u_deg > 1000), {
    parallel_for (intT j = 0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find all seagulls with center v and endpoint u
      for (intT k=0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset+k];
        if (u2 < i) {
          // TODO should take out -1, put in if to store only if > 0
          uintE num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
          writeAdd(&butterflies[eltsPerCacheLine*eti[u_offset + j]],num_butterflies); 
          writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)],num_butterflies); 
        }
        else break;
      }
    }
  }
retrieveCountsTimer.stop();
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

intT CountESortCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE* eti, intT curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
getWedgesTimer.start();
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, use_v, UWedgeCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;
getWedgesTimer.stop();
rehashWedgesTimer.start();
  using X = tuple<uintE,uintE>;
  if (cs.butterflies_seq_intt.n < 2*num_wedges_f) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*num_wedges_f);
    cs.butterflies_seq_intt.n = 2*num_wedges_f;
  }

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq());
  uintE* freq_arr = freq_pair.first;

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i] - 1;
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      cs.butterflies_seq_intt.A[2*j] = make_tuple(eti[offsetsU[wedges[j].v1] + wedges[j].j], num);
      cs.butterflies_seq_intt.A[2*j+1] = make_tuple(offsetsV[wedges[j].u] + wedges[j].k, num);
    }
  }
rehashWedgesTimer.stop();
retrieveCountsTimer.start();

  free(freq_arr);

  // now, we need to collate by our indices
  pair<uintE*, long> b_freq_pair = getFreqs(cs.butterflies_seq_intt.A, 2*num_wedges_f, uintETupleLt(), uintETupleEq());
  uintE* b_freq_arr = b_freq_pair.first;
const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    uintE num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]), num_freq, uintETupleAdd());
    // These are our butterfly counts
    butterflies[eltsPerCacheLine*get<0>(reduce)] += get<1>(reduce);
  }

  free(b_freq_arr);
retrieveCountsTimer.stop();
  return wedges_pair.second;
}

intT CountESort(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE* eti, intT curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;

  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, use_v, UWedgeCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;

getWedgesTimer.start();
  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq());
  uintE* freq_arr = freq_pair.first;
getWedgesTimer.stop();

rehashWedgesTimer.start();
const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i] - 1;
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      writeAdd(&butterflies[eltsPerCacheLine*eti[offsetsU[wedges[j].v1] + wedges[j].j]], num);
      writeAdd(&butterflies[eltsPerCacheLine*(offsetsV[wedges[j].u] + wedges[j].k)], num);
    }
  }
rehashWedgesTimer.stop();
  free(freq_arr);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

intT CountEHash(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

getWedgesTimer.start();
  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
getWedgesTimer.stop();

rehashWedgesTimer.start();
const intT eltsPerCacheLine = 64/sizeof(long);
  // TODO modularize w/hashce
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    //granular_for (j, 0, u_deg, (u_deg > 1000), {
    parallel_for (intT j = 0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find all seagulls with center v and endpoint u
      for (intT k=0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset + k];
        if (u2 < i) {
          uintE num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
          writeAdd(&butterflies[eltsPerCacheLine*eti[u_offset + j]],num_butterflies); 
          writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)],num_butterflies); 
        }
        else break;
      }
    }
  }
rehashWedgesTimer.stop();
  
  return next_idx;
}

intT CountEHashCE(CountESpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE* eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

getWedgesTimer.start();
  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
getWedgesTimer.stop();

rehashWedgesTimer.start();
  cs.butterflies_hash.resize(2*cs.wedges_hash.count());
  parallel_for (intT i = curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    parallel_for (intT j = 0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find all seagulls with center v and endpoint u
      for (long k=0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset + k];
        if (u2 < i) {
          uintE num_butterflies = cs.wedges_hash.find(i*nu + u2).second - 1;
          cs.butterflies_hash.insert(make_pair(eti[u_offset + j], num_butterflies));
          cs.butterflies_hash.insert(make_pair(v_offset + k, num_butterflies));
        }
        else break;
      }
    }
  }
rehashWedgesTimer.stop();

retrieveCountsTimer.start();
  intT num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.wedges_seq_intp);
const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for(long i=0; i < num_wedges_seq; ++i) {
      pair<uintE,uintE> butterflies_pair = cs.wedges_seq_intp.A[i];
      butterflies[eltsPerCacheLine*butterflies_pair.first] += butterflies_pair.second;
  }
retrieveCountsTimer.stop();
  
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

void CountEOrigCompactParallel(uintE* eti, uintE* butterflies, bipartiteCSR& GA, bool use_v) {
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
  
  //long* results = newA(long,eltsPerCacheLine*stepSize); //one entry per cache line
  //granular_for(i,0,stepSize,stepSize > 10000, {results[eltsPerCacheLine*i] = 0;});
  
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
	          //results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
	          wedges[shift+u2_idx]++;
	          if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = shift+u2_idx;
	        }
	        else break;
	      }
      }
      for(intT j=0; j < u_deg; ++j) {
        intT v = edgesU[u_offset+j];
        intT v_offset = offsetsV[v];
        intT v_deg = offsetsV[v+1] - v_offset;
        for(intT k=0; k < v_deg; ++k) {
          uintE u2_idx = edgesV[v_offset+k];
          if (u2_idx < i) {
          writeAdd(&butterflies[eltsPerCacheLine*eti[u_offset + j]], wedges[shift+u2_idx]-1);
          writeAdd(&butterflies[eltsPerCacheLine*(v_offset + k)], wedges[shift+u2_idx]-1);
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
}

uintE* CountE(uintE* eti, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long type=0) {
  const long nu = use_v ? GA.nu : GA.nv;

  const intT eltsPerCacheLine = 64/sizeof(long);
  uintE* butterflies = newA(uintE, eltsPerCacheLine*GA.numEdges);
  parallel_for(long i=0; i < GA.numEdges; ++i){
    butterflies[eltsPerCacheLine*i] = 0;
  }

  uintE* wedge_idxs = type == 5 ? nullptr : countWedgesScan(GA, use_v, true);
  CountESpace cs = CountESpace(type, nu);

  if (type == 5) CountEOrigCompactParallel(eti, butterflies, GA, use_v);
  else {
  // TODO clean up utils, do overflow stuff (double check youtube needs it)
  // TODO check correctness against each other
  // TODO why is hash so slow
  if (max_wedges >= num_wedges) {
    if (type == 0) CountESort(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    else if (type == 1) CountESortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    else if (type == 2) CountEHash(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    else if (type == 3) CountEHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    else if (type == 4) CountEHist(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
  }
  else {
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 0) curr_idx = CountESort(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    else if (type == 1) curr_idx = CountESortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    else if (type == 2) curr_idx = CountEHash(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    else if (type == 3) curr_idx = CountEHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    else if (type == 4) curr_idx = CountEHist(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    cs.clear();
  }
  }
  }

  if (type != 5) free(wedge_idxs);
  return butterflies;
}

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

intT CountEHist(_seq<uintE>& wedges, sparseAdditiveSet<uintE>& sizes, pbbsa::sequence<tuple<uintE, uintE>> tmp, pbbsa::sequence<tuple<uintE, uintE>> out,
  bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, edgeToIdx& eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  using X = tuple<uintE,uintE>;
  auto cons = UVertexPairIntCons(nu);

  pair<long, intT> wedges_pair = getWedges<uintE>(wedges, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
  uintE* wedges_list = wedges.A;
  long num_wedges_list = wedges_pair.first;
  intT next_idx = wedges_pair.second;

  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges_list);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nu*nu, tmp, out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  //TODO save this space too
  sizes.resize(wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    sizes.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

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
          uintE num_butterflies = sizes.find(cons(i, u2, v)).second - 1;
          writeAdd(&butterflies[eti(v, i)],num_butterflies); 
          writeAdd(&butterflies[eti(v, u2)],num_butterflies); 
        }
        else break;
      }
    }
  }

  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

intT CountESortCE(_seq<UWedge>& wedges_seq, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, edgeToIdx& eti, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UWedge>(wedges_seq, GA, use_v, UWedgeCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UWedge* wedges = wedges_seq.A;
  long num_wedges_f = wedges_pair.first;

  using X = tuple<uintE,uintE>;
  X* b_freqs = newA(X, 2*num_wedges_f);

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq());
  uintE* freq_arr = freq_pair.first;

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i] - 1;
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      b_freqs[2*j] = make_tuple(eti(wedges[j].u, wedges[j].v1), num);
      b_freqs[2*j+1] = make_tuple(eti(wedges[j].u, wedges[j].v2), num);
    }
  }

  free(freq_arr);

  // now, we need to collate by our indices
  pair<uintE*, long> b_freq_pair = getFreqs(b_freqs, 2*num_wedges_f, uintETupleLt(), uintETupleEq());
  uintE* b_freq_arr = b_freq_pair.first;

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    uintE num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(b_freqs[b_freq_arr[i-1]]), num_freq, uintETupleAdd());
    // These are our butterfly counts
    butterflies[get<0>(reduce)] += get<1>(reduce);
  }

  free(b_freq_arr);
  free(b_freqs);

  return wedges_pair.second;
}

intT CountESort(_seq<UWedge>& wedges_seq, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, edgeToIdx& eti, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UWedge>(wedges_seq, GA, use_v, UWedgeCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UWedge* wedges = wedges_seq.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq());
  uintE* freq_arr = freq_pair.first;

  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i] - 1;
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      writeAdd(&butterflies[eti(wedges[j].u, wedges[j].v1)], num);
      writeAdd(&butterflies[eti(wedges[j].u, wedges[j].v2)], num);
    }
  }

  free(freq_arr);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

intT CountEHash(sparseAdditiveSet<uintE>& wedges, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, edgeToIdx& eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(wedges, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);
  
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
          uintE num_butterflies = wedges.find(cons(i, u2, v)).second - 1;
          writeAdd(&butterflies[eti(v, i)],num_butterflies); 
          writeAdd(&butterflies[eti(v, u2)],num_butterflies); 
        }
        else break;
      }
    }
  }
  
  return next_idx;
}

intT CountEHashCE(sparseAdditiveSet<uintE>& wedges, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, edgeToIdx& eti, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  UVertexPairIntCons cons = UVertexPairIntCons(nu);
  intT next_idx = getWedgesHash(wedges, GA, use_v, cons, max_wedges, curr_idx, num_wedges, wedge_idxs);

  sparseAdditiveSet<uintE> butterflies_set = sparseAdditiveSet<uintE>(eti.numEdges,1,UINT_E_MAX);
  
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
          uintE num_butterflies = wedges.find(cons(i, u2, v)).second - 1;
          butterflies_set.insert(make_pair(eti(v, i), num_butterflies));
          butterflies_set.insert(make_pair(eti(v, u2), num_butterflies));
        }
        else break;
      }
    }
  }

  _seq<pair<uintE,uintE>> butterflies_seq = butterflies_set.entries();

  parallel_for(long i=0; i < butterflies_seq.n; ++i) {
      pair<uintE,uintE> butterflies_pair = butterflies_seq.A[i];
      butterflies[butterflies_pair.first] += butterflies_pair.second;
  }

  butterflies_seq.del();
  butterflies_set.del();
  
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

void CountEHash_helper(edgeToIdx& eti, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  initHashTimer.start();

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(eti.numEdges,1,UINT_E_MAX);

  if (max_wedges >= num_wedges) {
    if (type == 2) CountEHash(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    else CountEHashCE(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    wedges.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type ==2) curr_idx = CountEHash(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    else curr_idx = CountEHashCE(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    wedges.clear();
  }
  wedges.del();
}

void CountEHist_helper(edgeToIdx& eti, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  pbbsa::sequence<tuple<uintE, uintE>> tmp = pbbsa::sequence<tuple<uintE, uintE>>();
  pbbsa::sequence<tuple<uintE, uintE>> out = pbbsa::sequence<tuple<uintE, uintE>>();
  _seq<uintE> wedges_seq = _seq<uintE>(newA(uintE, nu), nu);
  sparseAdditiveSet<uintE> sizes = sparseAdditiveSet<uintE>(nu, 1, UINT_E_MAX);

  if (max_wedges >= num_wedges) {
    CountEHist(wedges_seq, sizes, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    wedges_seq.del();
    sizes.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    curr_idx = CountEHist(wedges_seq, sizes, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    sizes.clear();
  }
  wedges_seq.del();
  sizes.del();
}

void CountESort_helper(edgeToIdx& eti, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  _seq<UWedge> wedges_seq = _seq<UWedge>(newA(UWedge, nu), nu);

  if (max_wedges >= num_wedges) {
    if (type == 0) CountESort(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    else CountESortCE(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti);
    wedges_seq.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 0) curr_idx = CountESort(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
    else curr_idx = CountESortCE(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, eti, curr_idx);
  }
  wedges_seq.del();
}

uintE* CountE(edgeToIdx& eti, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long type=0) {
  const long nu = use_v ? GA.nu : GA.nv;

  uintE* butterflies = newA(uintE, GA.numEdges);
  parallel_for(long i=0; i < GA.numEdges; ++i){
    butterflies[i] = 0;
  }

  uintE* wedge_idxs = countWedgesScan(GA, use_v, true);

  // TODO clean up utils, do overflow stuff (double check youtube needs it)
  // TODO check correctness against each other
  // TODO why is hash so slow
  if (type == 2 || type == 3) CountEHash_helper(eti, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
  else if (type == 4) CountEHist_helper(eti, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
  else if (type == 0 || type == 1) CountESort_helper(eti, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);

  free(wedge_idxs);
  return butterflies;
}

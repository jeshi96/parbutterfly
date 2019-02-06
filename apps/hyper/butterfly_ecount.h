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

template <class vertex>
long CountEHistCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, edgeToIdx<vertex> eti, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  using X = tuple<uintE,uintE>;
  pair<pair<X*,long>, long> wedges_pair = getWedges2<X>(nu, V, U, UWedgeIntCons(nu), max_wedges, curr_idx, num_wedges,
    make_tuple(UINT_E_MAX, UINT_E_MAX), nonMaxTupleF(), true);
  X* wedges_list = wedges_pair.first.first;
  long num_wedges_list = wedges_pair.first.second;
  pbbsa::sequence<X> wedges_seq = pbbsa::sequence<X>(wedges_list, num_wedges_list);

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(num_wedges_list, 1, UINT_E_MAX);
  parallel_for(long i=0; i < num_wedges_list; ++i) {
    wedges.insert(make_pair(get<0>(wedges_list[i]), 1));
  }

  auto wedges_tuple = pbbsa::sparse_histogram_list<uintE, uintE>(wedges_seq, nu*nu + nu, wedges); //TODO this last one could be better
  tuple<uintE, pbbsa::sequentialHT<uintE,uintE>*>* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  X* b_freqs = newA(X, 2*num_wedges_list);
  long idx = 0;
  for (long i=0; i < wedge_freqs_n; ++i) {
    auto t = wedge_freqs[i];
    uintE v1 = get<0>(t) % nu;
    uintE v2 = get<0>(t) / nu;
    auto lst = get<1>(t)->compact();
    uintE num = get<0>(lst) - 1;
    for(long j=0; j < get<0>(lst); ++j) { //TODO can parallelize some prob
      uintE u = get<0>(get<1>(lst)[j]);
      b_freqs[idx++] = make_tuple(eti(u, v1), num);
      b_freqs[idx++] = make_tuple(eti(u, v2), num);
    }
    free(get<1>(lst));
    get<1>(t)->del();
  }
  free(wedge_freqs);
  free(wedges_list);

  pbbsa::sequence<tuple<uintE, uintE>> b_freqs_seq = pbbsa::sequence<tuple<uintE,uintE>>(b_freqs, 2*num_wedges_list);
  tuple<size_t,tuple<uintE, uintE>*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(b_freqs_seq, eti.num_edges, getAdd<uintE,uintE>, getAddReduce<uintE,uintE>);
  tuple<uintE,uintE>* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  /*uintE* butterflies = newA(uintE, nu*(nv-1)+nu-1);
  parallel_for(long i=0;i<nu*(nv-1)+nu-1;++i){
    butterflies[i] = 0;
  }*/
  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  free(b_freqs);

  return wedges_pair.second;
}

template <class vertex>
long CountEHist(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, edgeToIdx<vertex> eti, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  using X = tuple<uintE,uintE>;
  pair<pair<X*,long>, long> wedges_pair = getWedges2<X>(nu, V, U, UWedgeIntCons(nu), max_wedges, curr_idx, num_wedges,
    make_tuple(UINT_E_MAX, UINT_E_MAX), nonMaxTupleF(), true);
  X* wedges_list = wedges_pair.first.first;
  long num_wedges_list = wedges_pair.first.second;
  pbbsa::sequence<X> wedges_seq = pbbsa::sequence<X>(wedges_list,num_wedges_list);

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(num_wedges_list, 1, UINT_E_MAX);
  parallel_for(long i=0; i < num_wedges_list; ++i) {
    wedges.insert(make_pair(get<0>(wedges_list[i]), 1));
  }

  auto wedges_tuple = pbbsa::sparse_histogram_list<uintE, uintE>(wedges_seq, nu*nu + nu, wedges); //TODO this last one could be better
  tuple<uintE, pbbsa::sequentialHT<uintE,uintE>*>* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);


  /*uintE* butterflies = newA(uintE, nu*(nv-1)+nu-1);
  parallel_for(long i=0;i<nu*(nv-1)+nu-1;++i){
    butterflies[i] = 0;
  }*/

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    auto t = wedge_freqs[i];
    uintE v1 = get<0>(t) % nu;
    uintE v2 = get<0>(t) / nu;
    auto lst = get<1>(t)->compact();
    uintE num = get<0>(lst) - 1;
    parallel_for(long j=0; j < get<0>(lst); ++j) {
      uintE u = get<0>(get<1>(lst)[j]);
      writeAdd(&butterflies[eti(u, v1)], num);
      writeAdd(&butterflies[eti(u, v2)], num);
    }
    free(get<1>(lst));
    get<1>(t)->del();
  }

  free(wedge_freqs);
  free(wedges_list);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

template <class vertex>
long CountESortCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, edgeToIdx<vertex> eti, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  pair<pair<UWedge*,long>, long> wedges_pair = getWedges2<UWedge>(nu, V, U, UWedgeCons(), max_wedges, curr_idx, num_wedges,
    UWedge(UINT_E_MAX, UINT_E_MAX, UINT_E_MAX), nonEmptyUWF(), true);
  UWedge* wedges = wedges_pair.first.first;
  long num_wedges_f = wedges_pair.first.second;
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
  free(wedges);

  /*uintE* butterflies = newA(uintE, nu*(nv-1)+nu-1);
  parallel_for(long i=0;i<nu*(nv-1)+nu-1;++i){
    butterflies[i] = 0; 
  }*/

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

template <class vertex>
long CountESort(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, edgeToIdx<vertex> eti, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  pair<pair<UWedge*,long>, long> wedges_pair = getWedges2<UWedge>(nu, V, U, UWedgeCons(), max_wedges, curr_idx, num_wedges,
    UWedge(UINT_E_MAX, UINT_E_MAX, UINT_E_MAX), nonEmptyUWF(), true);
  UWedge* wedges = wedges_pair.first.first;
  long num_wedges_f = wedges_pair.first.second;

  /*uintE* butterflies = newA(uintE, nu*(nv-1)+nu-1);
  parallel_for(long i=0;i<nu*(nv-1)+nu-1;++i){
    butterflies[i] = 0; 
  }*/

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq());
  uintE* freq_arr = freq_pair.first;

  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i] - 1;
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      //fflush(stdout);
      writeAdd(&butterflies[eti(wedges[j].u, wedges[j].v1)], num);
      writeAdd(&butterflies[eti(wedges[j].u, wedges[j].v2)], num);
    }
  }

  free(freq_arr);
  free(wedges);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

template<class vertex, class wedgeCons>
sparseSet<sparseSet<uintE>*> allocateWedgesHash(sparseAdditiveSet<uintE>& wedges, const long nv, vertex* V, vertex* U, wedgeCons cons,
  long curr_idx, long next_idx) {
  using T=sparseSet<uintE>*;
  _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
  sparseSet<T> cwedges = sparseSet<T>(wedge_freqs.n, 1, NULL);
  // Allocate the space for all of our wedge lists
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
    pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
    T cset = new sparseSet<uintE>(wedge_freq_pair.second, 1, UINT_E_MAX);
    cwedges.insert(make_pair(wedge_freq_pair.first, cset));
  }

  // Count number of wedges by their key
  parallel_for (long i = curr_idx; i < next_idx; ++i) {
    const vertex u = U[i];
    parallel_for (long j = 0; j < u.getOutDegree(); ++j) {
      const uintE v_idx = u.getOutNeighbor(j);
      const vertex v = V[v_idx];
      // Find all seagulls with center v and endpoint u
      for (long k=0; k < v.getOutDegree(); ++k) {
        const uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx < i) {
          sparseSet<uintE>* wedge_pair = cwedges.find(cons(i, u2_idx, v_idx)).second;
          wedge_pair->insert(pair<uintE,uintE>(v_idx,1));
        }
      }
    }
  }

  return cwedges;
}

template <class writeAddOp, class vertex>
void countButterfliesEHash(sparseSet<sparseSet<uintE>*>& cwedges, const long nu, const long nv, writeAddOp op, edgeToIdx<vertex> eti) {
  auto cwedge_freqs = cwedges.entries();
  // Retrieve count on each key; that number minus 1 is the number of butterflies
  parallel_for (long i=0; i < cwedge_freqs.n; ++i) {
    uintE key = cwedge_freqs.A[i].first;
    sparseSet<uintE>* centers = cwedge_freqs.A[i].second;
    _seq<pair<uintE,uintE>> centers_seq = centers->entries();

    uintE u2 = key % nu;
    uintE u1 = key / nu;
    uintE num_butterflies = centers_seq.n - 1;
    parallel_for(long j=0; j < centers_seq.n; ++j) {
      uintE v = centers_seq.A[j].first;
      op(eti(v,u2), num_butterflies);
      op(eti(v,u1), num_butterflies);
    }

    centers_seq.del();
    centers->del();
  }
  cwedge_freqs.del();
}

template <class vertex>
long CountEHash(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, edgeToIdx<vertex> eti, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(min(num_wedges,max_wedges), 1, UINT_E_MAX);
  long next_idx = getWedgesHash2(wedges, nu, V, U, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, true);
  // TODO idk if this nested stuff is best -- I needed a linkedlist append kind of thing but parallel
  sparseSet<sparseSet<uintE>*> cwedges = allocateWedgesHash(wedges,nv,V,U,UVertexPairIntCons(nu), curr_idx, next_idx);
  
  //TODO this should be init to # of edges --> need a good way to index edges???? unless we do it on all pairs
  /*uintE* butterflies = newA(uintE, nu*(nv-1)+nu-1);
  parallel_for(long i=0;i<nu*(nv-1)+nu-1;++i){
    butterflies[i] = 0;
  }*/

  countButterfliesEHash(cwedges,nu,nv, writeAddArr<uintE>(butterflies), eti);
  
  cwedges.del();
  wedges.del();
  return next_idx;
}

template <class vertex>
long CountEHashCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, edgeToIdx<vertex> eti, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(min(num_wedges,max_wedges), 1, UINT_E_MAX);
  long next_idx = getWedgesHash2(wedges, nu, V, U, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, true);
  // TODO idk if this nested stuff is best -- I needed a linkedlist append kind of thing but parallel
  sparseSet<sparseSet<uintE>*> cwedges = allocateWedgesHash(wedges,nv,V,U,UVertexPairIntCons(nu), curr_idx, next_idx);
  
  //TODO this should be init to # of edges --> need a good way to index edges???? unless we do it on all pairs
  //sparseAdditiveSet<uintE> butterflies_set = sparseAdditiveSet<uintE>(nu*(nv-1)+nu-1,1,UINT_E_MAX);
  sparseAdditiveSet<uintE> butterflies_set = sparseAdditiveSet<uintE>(eti.num_edges, 1, UINT_E_MAX);

  countButterfliesEHash(cwedges,nu,nv, writeAddSet<uintE>(butterflies_set), eti);

  cwedges.del();
  wedges.del();

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

template <class vertex>
uintE* CountE(edgeToIdx<vertex> eti, bipartiteGraph<vertex> GA, bool use_v, long num_wedges, long max_wedges, long type=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  //auto eti = edgeToIdx<vertex>(GA, use_v, max_wedges);

  uintE* butterflies = newA(uintE, eti.num_edges);
  parallel_for(long i=0; i < eti.num_edges; ++i){
    butterflies[i] = 0;
  }
  // TODO clean up utils, do overflow stuff (double check youtube needs it)
  // TODO check correctness against each other
  // TODO why is hash so slow
  if (max_wedges >= num_wedges) {
    if (type == 2) CountEHash(GA, use_v, num_wedges, butterflies, max_wedges, eti);
    else if (type == 3) CountEHashCE(GA, use_v, num_wedges, butterflies, max_wedges, eti);
    else if (type == 0) CountESort(GA, use_v, num_wedges, butterflies, max_wedges, eti);
    else if (type==1) CountESortCE(GA, use_v, num_wedges, butterflies, max_wedges, eti);
    else if (type==4) CountEHist(GA, use_v, num_wedges, butterflies, max_wedges, eti);
    else CountEHistCE(GA, use_v, num_wedges, butterflies, max_wedges, eti);
    return butterflies;
  }
  long curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 2) curr_idx = CountEHash(GA, use_v, num_wedges, butterflies, max_wedges, eti, curr_idx);
    else if (type == 3) curr_idx = CountEHashCE(GA, use_v, num_wedges, butterflies, max_wedges, eti, curr_idx);
    else if (type == 0) curr_idx = CountESort(GA, use_v, num_wedges, butterflies, max_wedges, eti, curr_idx);
    else if (type==1) curr_idx = CountESortCE(GA, use_v, num_wedges, butterflies, max_wedges, eti, curr_idx);
    else if (type==4) curr_idx = CountEHist(GA, use_v, num_wedges, butterflies, max_wedges, eti, curr_idx);
    else curr_idx = CountEHistCE(GA, use_v, num_wedges, butterflies, max_wedges, eti, curr_idx);
  }
  return butterflies;
}

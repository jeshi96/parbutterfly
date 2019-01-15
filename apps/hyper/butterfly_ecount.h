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
uintE* CountEHistCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  using X = tuple<uintE,uintE>;
  X* wedges_list = getWedges<X>(nv, V, WedgeIntCons(nu));
  pbbsa::sequence<X> wedges_seq = pbbsa::sequence<X>(wedges_list,num_wedges);

  auto wedges_tuple = pbbsa::sparse_histogram_list<uintE, uintE>(wedges_seq, nu*nu + nu, nv); //TODO this last one could be better
  tuple<uintE, pbbsa::sequentialHT<uintE,uintE>*>* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  X* b_freqs = newA(X, 2*num_wedges);
  long idx = 0;
  for (long i=0; i < wedge_freqs_n; ++i) {
    auto t = wedge_freqs[i];
    uintE v1 = get<0>(t) % nu;
    uintE v2 = get<0>(t) / nu;
    auto lst = get<1>(t)->compact();
    uintE num = get<0>(lst) - 1;
    for(long j=0; j < get<0>(lst); ++j) {
      uintE u = get<0>(get<1>(lst)[j]);
      b_freqs[idx++] = make_tuple(nv*u + v1, num);
      b_freqs[idx++] = make_tuple(nv*u + v2, num);
    }
    free(get<1>(lst));
    get<1>(t)->del();
  }
  free(wedge_freqs);
  free(wedges_list);

  pbbsa::sequence<tuple<uintE, uintE>> b_freqs_seq = pbbsa::sequence<tuple<uintE,uintE>>(b_freqs, 2*num_wedges);
  tuple<size_t,tuple<uintE, uintE>*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(b_freqs_seq, nu*nv, getAdd<uintE,uintE>, getAddReduce<uintE,uintE>);
  tuple<uintE,uintE>* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  uintE* butterflies = newA(uintE, nu*nv);
  parallel_for(long i=0;i<nu*nv;++i){
    butterflies[i] = 0;
  }
  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[get<0>(butterflies_l[i])] = get<1>(butterflies_l[i]);
  }

  free(b_freqs);

  return butterflies;
}

template <class vertex>
uintE* CountEHist(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  using X = tuple<uintE,uintE>;
  X* wedges_list = getWedges<X>(nv, V, WedgeIntCons(nu));
  pbbsa::sequence<X> wedges_seq = pbbsa::sequence<X>(wedges_list,num_wedges);

  auto wedges_tuple = pbbsa::sparse_histogram_list<uintE, uintE>(wedges_seq, nu*nu + nu, nv); //TODO this last one could be better
  tuple<uintE, pbbsa::sequentialHT<uintE,uintE>*>* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);


  uintE* butterflies = newA(uintE, nu*nv);
  parallel_for(long i=0;i<nu*nv;++i){
    butterflies[i] = 0;
  }

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    auto t = wedge_freqs[i];
    uintE v1 = get<0>(t) % nu;
    uintE v2 = get<0>(t) / nu;
    auto lst = get<1>(t)->compact();
    uintE num = get<0>(lst) - 1;
    parallel_for(long j=0; j < get<0>(lst); ++j) {
      uintE u = get<0>(get<1>(lst)[j]);
      writeAdd(&butterflies[nv * u + v1], num);
      writeAdd(&butterflies[nv * u + v2], num);
    }
    free(get<1>(lst));
    get<1>(t)->del();
  }

  free(wedge_freqs);
  free(wedges_list);
  return butterflies;
}

//********************************************************************************************
//********************************************************************************************

template <class vertex>
uintE* CountESortCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  Wedge* wedges = getWedges<Wedge>(nv, V, WedgeCons());
  using X = tuple<uintE,uintE>;
  X* b_freqs = newA(X, 2*num_wedges);

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges, WedgeCmp(nu), WedgeEq());
  uintE* freq_arr = freq_pair.first;

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i] - 1;
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      b_freqs[j] = make_tuple(nv * wedges[j].u + wedges[j].v1, num);
      b_freqs[num_wedges + j] = make_tuple(nv * wedges[j].u + wedges[j].v2, num);
    }
  }

  free(freq_arr);
  free(wedges);

  uintE* butterflies = newA(uintE, nu*nv);
  parallel_for(long i=0;i<nu*nv;++i){
    butterflies[i] = 0; 
  }

  // now, we need to collate by our indices
  pair<uintE*, long> b_freq_pair = getFreqs(b_freqs, 2*num_wedges, uintETupleLt(), uintETupleEq());
  uintE* b_freq_arr = b_freq_pair.first;

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    uintE num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(b_freqs[b_freq_arr[i-1]]), num_freq, uintETupleAdd());
    // These are our butterfly counts
    butterflies[get<0>(reduce)] = get<1>(reduce);
  }

  free(b_freq_arr);
  free(b_freqs);

  return butterflies;
}

template <class vertex>
uintE* CountESort(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  Wedge* wedges = getWedges<Wedge>(nv, V, WedgeCons());

  uintE* butterflies = newA(uintE, nu*nv);
  parallel_for(long i=0;i<nu*nv;++i){
    butterflies[i] = 0; 
  }

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges, WedgeCmp(nu), WedgeEq());
  uintE* freq_arr = freq_pair.first;

  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i] - 1;
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      writeAdd(&butterflies[nv * wedges[j].u + wedges[j].v1], num);
      writeAdd(&butterflies[nv * wedges[j].u + wedges[j].v2], num);
    }
  }

  free(freq_arr);
  free(wedges);
  return butterflies;
}

//********************************************************************************************
//********************************************************************************************

template<class vertex, class wedgeCons>
sparseListSet<uintE> allocateWedgesHash(sparseAdditiveSet<uintE>& wedges, const long nv, vertex* V, wedgeCons cons) {
  _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
  sparseListSet<uintE> cwedges = sparseListSet<uintE>(wedge_freqs.n, 1);
  // Allocate the space for all of our wedge lists
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
    pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
    sparseKeySet<uintE>* cset = new sparseKeySet<uintE>(wedge_freq_pair.second,1,UINT_E_MAX,wedge_freq_pair.first);
    cwedges.insert(cset);
  }

  // Count number of wedges by their key
  for (long i = 0; i < nv; ++i) {
    const vertex v = V[i];
    const uintE v_deg = v.getOutDegree();
    for (long j = 0; j < v_deg; ++j) {
      for (long k = j+1; k < v_deg; ++k) {
        sparseKeySet<uintE>* wedge_pair = cwedges.find(cons(v.getOutNeighbor(j), v.getOutNeighbor(k), i));
        wedge_pair->insert(pair<uintE,uintE>(i,1));
      }
    }
  }

  return cwedges;
}

template <class writeAddOp>
void countButterfliesEHash(sparseListSet<uintE>& cwedges, const long nu, const long nv, writeAddOp op) {
  _seq<sparseKeySet<uintE>*> cwedge_freqs = cwedges.entries();
  // Retrieve count on each key; that number minus 1 is the number of butterflies
  parallel_for (long i=0; i < cwedge_freqs.n; ++i) {
    sparseKeySet<uintE>* centers = cwedge_freqs.A[i];
    _seq<pair<uintE,uintE>> centers_seq = centers->entries();

    uintE u2 = centers->mas_key % nu;
    uintE u1 = centers->mas_key / nu;
    uintE num_butterflies = centers_seq.n - 1;
    parallel_for(long j=0; j < centers_seq.n; ++j) {
      uintE v = centers_seq.A[j].first;
      op(v*nv + u2, num_butterflies);
      op(v*nv + u1, num_butterflies);
    }

    centers_seq.del();
    centers->del();
  }
  cwedge_freqs.del();
}

template <class vertex>
uintE* CountEHash(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(num_wedges, 1, UINT_E_MAX);
  getWedgesHash(wedges, nv, V, VertexPairIntCons(nu));
  // TODO idk if this nested stuff is best -- I needed a linkedlist append kind of thing but parallel
  sparseListSet<uintE> cwedges = allocateWedgesHash(wedges,nv,V,VertexPairIntCons(nu));
  
  //TODO this should be init to # of edges --> need a good way to index edges???? unless we do it on all pairs
  uintE* butterflies = newA(uintE, nu*nv);
  parallel_for(long i=0;i<nu*nv;++i){
    butterflies[i] = 0;
  }

  countButterfliesEHash(cwedges,nu,nv, writeAddArr<uintE>(butterflies));
  
  cwedges.del();
  wedges.del();
  return butterflies;
}

template <class vertex>
uintE* CountEHashCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(num_wedges, 1, UINT_E_MAX);
  getWedgesHash(wedges, nv, V, VertexPairIntCons(nu));
  // TODO idk if this nested stuff is best -- I needed a linkedlist append kind of thing but parallel
  sparseListSet<uintE> cwedges = allocateWedgesHash(wedges,nv,V,VertexPairIntCons(nu));
  
  //TODO this should be init to # of edges --> need a good way to index edges???? unless we do it on all pairs
  sparseAdditiveSet<uintE> butterflies_set = sparseAdditiveSet<uintE>(nu*nv,1,UINT_E_MAX);

  countButterfliesEHash(cwedges,nu,nv, writeAddSet<uintE>(butterflies_set));

  cwedges.del();
  wedges.del();

  _seq<pair<uintE,uintE>> butterflies_seq = butterflies_set.entries();

  uintE* butterflies = newA(uintE, nu*nv);
  parallel_for(long i=0;i<nu*nv;++i){
    butterflies[i] = 0;
  }

  parallel_for(long i=0; i < butterflies_seq.n; ++i) {
      pair<uintE,uintE> butterflies_pair = butterflies_seq.A[i];
      butterflies[butterflies_pair.first] = butterflies_pair.second;
  }

  butterflies_seq.del();
  butterflies_set.del();

  return butterflies;
}
//********************************************************************************************
//********************************************************************************************

template <class vertex>
uintE* CountE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, long type=0) {
  if (type == 0) return CountEHash(GA, use_v, num_wedges);
  else if (type == 1) return CountEHashCE(GA, use_v, num_wedges);
  else if (type == 2) return CountESort(GA, use_v, num_wedges);
  else if (type==3) return CountESortCE(GA, use_v, num_wedges);
  else if (type==4) return CountEHist(GA, use_v, num_wedges);
  else return CountEHistCE(GA, use_v, num_wedges);
}
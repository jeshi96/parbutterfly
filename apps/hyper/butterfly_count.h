#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "vertex.h"
#include "sequence.h"
#include "sparseSet.h"
#include "sampleSort.h"
#include "../../lib/histogram.h"
#include "../../lib/gbbs-histogram.h"

#include "butterfly_utils.h"

using namespace std;

// Not parallel
void storeButterfliesSortCE_seq(uintE* butterflies, VertexPair* wedges, uintE* wedge_freqs_f, 
                            uintE* par_idxs_f, long i, bool use_v1) {
  // Retrieve the frequency counts for each distinct key by taking the difference between
  // the wedge_freqs_f indices
  // Because of how we set up the par_idxs_f, we know that every element in this range has the same
  // first vertex
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // Store butterfly counts with each vertex in U
  for (long idx = par_idxs_f[i]; idx < par_idxs_f[i+1]; ++idx) { //TODO this is without reduce -- do thresholding
      uintE num_butterflies = wedge_freqs_f[idx+1] - wedge_freqs_f[idx];
      long wedge_idx = wedge_freqs_f[idx];
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      
      VertexPair vs = wedges[wedge_idx];

      if (use_v1) butterflies[vs.v1] += num_butterflies;
      else butterflies[vs.v2] += num_butterflies;
    }
}

// Parallel
void storeButterfliesSortCE(uintE* butterflies, long nu, VertexPair* wedges, uintE* wedge_freqs_f, 
                            uintE* par_idxs_f, long i, bool use_v1) {
  // Retrieve the frequency counts for each distinct key by taking the difference between
  // the wedge_freqs_f indices
  // Because of how we set up the par_idxs_f, we know that every element in this range has the same
  // first vertex, and a different second vertex
  // Store the number of butterflies on each key by the second vertex (no parallel issues b/c distinct)
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // The sum gives the number of butterflies on the first vertex
  uintE* butterflies_v = newA(uintE, nu);
  parallel_for(long i=0; i < nu; ++i) {butterflies_v[i] = 0; }
  parallel_for (long idx = par_idxs_f[i]; idx < par_idxs_f[i+1]; ++idx) {
      uintE num_butterflies = wedge_freqs_f[idx+1] - wedge_freqs_f[idx];
      long wedge_idx = wedge_freqs_f[idx];
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      
      VertexPair vs = wedges[wedge_idx];

      if (use_v1) butterflies_v[vs.v2] = num_butterflies;
      else butterflies_v[vs.v1] = num_butterflies;
  }
  long sum = sequence::plusReduce(butterflies_v, nu);
  VertexPair start_vs = wedges[wedge_freqs_f[par_idxs_f[i]]];
  if (use_v1) butterflies[start_vs.v1] += sum;
  else butterflies[start_vs.v2] += sum;
}

// Retrieve frequency counts for all wedges with the same key, and for all wedges with the same first vertex
// (if use_v1; otherwise, second vertex -- in comments, we assume first vertex)
// First, retrieve a list of indices where consecutive wedges have different keys
pair<uintE*, long> getWedgeFreqs(long nv, long nu, VertexPair* wedges, long num_wedges, bool use_v1) {
    if (use_v1) return getFreqs(wedges, num_wedges, VertexPairCmp(nu), VertexPairEq());
    return getFreqs(wedges, num_wedges, VertexPairCmp2(nu), VertexPairEq());
    //TODO check if num_wedges > 2^32 to see if overflow; if it is, compile w/env var LONG/EDGELONG
}

void countButterfliesSortCE(uintE* butterflies, long nv, long nu, VertexPair* wedges, uintE* wedge_freqs_f, 
            long num_wedge_freqs_f, bool use_v1) {
//timer t;
//t.start();
  // Given this list wedge_freqs_f, retrieve a list of indices for wedge_freqs_f such that
  // consecutive wedges have different first vertices
  // That is to say, the list we construct par_idxs_f is a list of indices on wedge_freqs_f,
  // such that wedges[wedge_freqs_f[par_idxs_f[i]]] and wedges[wedge_freqs_f[par_idxs_f[i+1]]]
  // have different first vertices
  pair<uintE*, long> par_idxs_pair = getFreqs(wedge_freqs_f,num_wedge_freqs_f-1, 
     NestedVPCmp(wedges,use_v1), NestedVPEq(wedges,use_v1), false);
  uintE* par_idxs_f = par_idxs_pair.first;
  long num_par_idxs_f = par_idxs_pair.second;

  // Use these two lists to retrieve the number of butterflies in a cache-efficient manner
  // Start by iterating through our par_idxs_f list
  parallel_for (long i = 0; i < num_par_idxs_f-1; ++i) {
    // The difference between consecutive elements tells us how many distinct values in wedge_freqs_f
    // have the same first vertex; use this to threshold whether to parallelize or not
    long range = par_idxs_f[i+1] - par_idxs_f[i];
    if (range > 10000) storeButterfliesSortCE(butterflies,nu,wedges,wedge_freqs_f,par_idxs_f,i,use_v1);
    else storeButterfliesSortCE_seq(butterflies,wedges,wedge_freqs_f,par_idxs_f,i,use_v1); //TODO set threshold var
  }

//t.stop();
//t.reportTotal("\tcountButterfliesSortCE:");
}

void countButterfliesSort(uintE* butterflies, VertexPair* wedges, uintE* wedge_freqs_f, long num_wedge_freqs_f) {
//timer t;
//t.start();
  // Then, retrieve the frequency counts for each distinct key by taking the difference between
  // these indices
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // Store butterfly counts with each vertex in U
  parallel_for (long i = 1; i < num_wedge_freqs_f; ++i) {
    uintE num_butterflies = wedge_freqs_f[i] - wedge_freqs_f[i-1];
    long wedge_idx = wedge_freqs_f[i-1];

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    VertexPair vs = wedges[wedge_idx];

    writeAdd(&butterflies[vs.v1],num_butterflies); 
    writeAdd(&butterflies[vs.v2],num_butterflies);
  }
//t.stop();
//t.reportTotal("\tcountButterfliesSort:");
}

// This is the original compute function, without the more cache-efficient sorting method
template <class vertex>
uintE* CountSort(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  VertexPair* wedges = getWedges<VertexPair>(nv, V, VertexPairCons());

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0; 
  }

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
//timer t2;
//t2.start();
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges, true);
//t2.stop();
//t2.reportTotal("\tgetWedgeFreqs:");

  countButterfliesSort(butterflies, wedges, freq_pair.first, freq_pair.second);

  free(freq_pair.first);

  free(wedges);
  return butterflies;
}

// This is the new compute function, with cache efficient sorting
template <class vertex>
uintE* CountSortCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  VertexPair* wedges = getWedges<VertexPair>(nv, V, VertexPairCons());

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0; 
  }

  // Retrieve butterflies + store with first vertex
//timer t2;
//t2.start();
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges, true);
//t2.stop();
//t2.reportTotal("\tgetWedgeFreqs:");
  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair.first, freq_pair.second , true);
  free(freq_pair.first);

  // Retrieve butterflies + store with second vertex
//timer t3;
//t3.start();
  pair<uintE*, long> freq_pair2 = getWedgeFreqs(nv, nu, wedges, num_wedges, false);
//t3.stop();
//t3.reportTotal("\tgetWedgeFreqs:");
  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair2.first, freq_pair2.second , false);
  free(freq_pair2.first);

  free(wedges);
  return butterflies;

}

//********************************************************************************************
//********************************************************************************************


template <class vertex>
uintE* CountHash(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  //float f = ((float)num_wedges)/((float) (nu*nu+nu));
  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(num_wedges, 1, UINT_E_MAX);
  getWedgesHash(wedges,nv, V, VertexPairIntCons(nu));
  
  _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }
//timer t;
//t.start();
  // Retrieve count on each key; that number choose 2 is the number of butterflies
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
    pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
    uintE num_butterflies = wedge_freq_pair.second;
    uintE v2 = wedge_freq_pair.first % nu;
    uintE v1 = wedge_freq_pair.first / nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    writeAdd(&butterflies[v1],num_butterflies);
    writeAdd(&butterflies[v2],num_butterflies);
  }
//t.stop();
//t.reportTotal("\tcountButterfliesHash:");

  wedge_freqs.del();
  wedges.del();
  return butterflies;
}

template <class vertex>
uintE* CountHashCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  //float f = ((float)num_wedges)/((float) (nu*nu+nu));
  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(num_wedges,1,UINT_E_MAX);
  getWedgesHash(wedges, nv, V, VertexPairIntCons(nu));
  
  _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
  sparseAdditiveSet<uintE> butterflies_set = sparseAdditiveSet<uintE>(nu,1,UINT_E_MAX);

//timer t;
//t.start();
  // Retrieve count on each key; that number choose 2 is the number of butterflies
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
    pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
    uintE num_butterflies = wedge_freq_pair.second;
    uintE v2 = wedge_freq_pair.first % nu;
    uintE v1 = wedge_freq_pair.first / nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    butterflies_set.insert(pair<uintE,uintE>(v1, num_butterflies));
    butterflies_set.insert(pair<uintE,uintE>(v2, num_butterflies));
  }

  _seq<pair<uintE,uintE>> butterflies_seq = butterflies_set.entries();
  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }
  parallel_for(long i=0; i < butterflies_seq.n; ++i) {
    pair<uintE, uintE> butterfly_pair = butterflies_seq.A[i];
    butterflies[butterfly_pair.first] = butterfly_pair.second;
  }
//t.stop();
//t.reportTotal("\tcountButterfliesHashCE:");
  
  butterflies_seq.del();
  butterflies_set.del();
  wedge_freqs.del();
  wedges.del();
  return butterflies;
}

//********************************************************************************************
//********************************************************************************************

// TODO use more efficient hist from laxman
template <class vertex>
uintE* CountHist(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, bool gbbs) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  uintE* wedges_list = getWedges<uintE>(nv, V, VertexPairIntCons(nu));
  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.
  tuple<uintE, uintE>* wedge_freqs;
  size_t wedge_freqs_n;

//timer t2;
//t2.start();
  if (gbbs) {
    tuple<size_t,tuple<uintE, uintE>*> wedges_tuple = 
      pbbsa::gbbs::sparse_histogram<uintE, uintE>(wedges_seq, nu*nu + nu);
    wedge_freqs = get<1>(wedges_tuple);
     wedge_freqs_n = get<0>(wedges_tuple);
  }
  else {
    tuple<size_t,tuple<uintE, uintE>*> wedges_tuple = 
      pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nu*nu + nu);
    wedge_freqs = get<1>(wedges_tuple);
    wedge_freqs_n = get<0>(wedges_tuple);
  }
//t2.stop();
//t2.reportTotal("\thistWedges:");

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }

//timer t;
//t.start();
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    tuple<uintE,uintE> wedge_freq_pair = wedge_freqs[i];
    uintE num_butterflies = get<1>(wedge_freq_pair);
    uintE wedge_freq_pair_first = get<0>(wedge_freq_pair);

    uintE v2 = wedge_freq_pair_first % nu;
    uintE v1 = wedge_freq_pair_first / nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    writeAdd(&butterflies[v1],num_butterflies);
    writeAdd(&butterflies[v2],num_butterflies);
  }
//t.stop();
//t.reportTotal("\tcountButterfliesHist:");

  free(wedge_freqs);
  free(wedges_list);
  return butterflies;
}

// TODO use more efficient hist from laxman
template <class vertex>
uintE* CountHistCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  uintE* wedges_list = getWedges<uintE>(nv, V, VertexPairIntCons(nu));
  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.
//timer t2;
//t2.start();
  tuple<size_t,tuple<uintE, uintE>*> wedges_tuple = 
    pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nv*nu + nu);
  tuple<uintE, uintE>* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);
using X = tuple<uintE,uintE>;
  X* wedge_freqs_i = newA(X, 2*wedge_freqs_n);
  parallel_for(long i = 0; i < wedge_freqs_n; i++) {
    tuple<uintE,uintE> wedge_freq_pair = wedge_freqs[i];
    uintE num = get<1>(wedge_freq_pair);
    num = (num * (num-1)) / 2;
    uintE wedge_num = get<0>(wedge_freq_pair);
    wedge_freqs_i[2*i] = make_tuple(wedge_num % nu, num);
    wedge_freqs_i[2*i + 1] = make_tuple(wedge_num / nu, num);
  }
  pbbsa::sequence<tuple<uintE, uintE>> wedge_freqs_i_seq = pbbsa::sequence<tuple<uintE,uintE>>(wedge_freqs_i,2*wedge_freqs_n);
  tuple<size_t,tuple<uintE, uintE>*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(wedge_freqs_i_seq,nu,getAdd<uintE,uintE>, getAddReduce<uintE,uintE>);
  tuple<uintE,uintE>* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);
//t2.stop();
//t2.reportTotal("\thistWedgesCE:");

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }
  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[get<0>(butterflies_l[i])] = get<1>(butterflies_l[i]);
  }

  free(wedge_freqs);
  free(wedges_list);

  return butterflies;
}


//********************************************************************************************
//********************************************************************************************

template <class vertex>
uintE* Count(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, long type=0) {
  if (type == 0) return CountSort(GA, use_v, num_wedges);
  else if (type == 1) return CountSortCE(GA, use_v, num_wedges);
  else if (type == 2) return CountHash(GA, use_v, num_wedges);
  else if (type == 3) return CountHashCE(GA, use_v, num_wedges);
  else if(type == 4) return CountHist(GA, use_v, num_wedges, false);
  else if(type == 5) return CountHist(GA, use_v, num_wedges, true);
  return CountHistCE(GA, use_v, num_wedges);

}
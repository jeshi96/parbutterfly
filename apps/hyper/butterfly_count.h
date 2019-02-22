#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <time.h>

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

timer nextWedgeTimer, hashInsertTimer, numButterfliesHashInsertTimer, getWedgesFromHashTimer, initHashTimer;
timer seqLoopTimer, seqWriteTimer;

#include "butterfly_utils.h"

using namespace std;



// Not parallel
void storeButterfliesSortCE_seq(uintE* butterflies, UVertexPair* wedges, uintE* wedge_freqs_f, 
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
      
      UVertexPair vs = wedges[wedge_idx];

      if (use_v1) butterflies[vs.v1] += num_butterflies;
      else butterflies[vs.v2] += num_butterflies;
    }
}

// Parallel
void storeButterfliesSortCE(uintE* butterflies, long nu, UVertexPair* wedges, uintE* wedge_freqs_f, 
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
      
      UVertexPair vs = wedges[wedge_idx];

      if (use_v1) butterflies_v[vs.v2] += num_butterflies;
      else butterflies_v[vs.v1] += num_butterflies;
  }
  long sum = sequence::plusReduce(butterflies_v, nu);
  UVertexPair start_vs = wedges[wedge_freqs_f[par_idxs_f[i]]];
  if (use_v1) butterflies[start_vs.v1] += sum;
  else butterflies[start_vs.v2] += sum;

  free(butterflies_v);
}

// Retrieve frequency counts for all wedges with the same key, and for all wedges with the same first vertex
// (if use_v1; otherwise, second vertex -- in comments, we assume first vertex)
// First, retrieve a list of indices where consecutive wedges have different keys
pair<uintE*, long> getWedgeFreqs(long nv, long nu, UVertexPair* wedges, long num_wedges, bool use_v1) {
    if (use_v1) return getFreqs(wedges, num_wedges, UVertexPairCmp(nu), UVertexPairEq());
    return getFreqs(wedges, num_wedges, UVertexPairCmp2(nu), UVertexPairEq());
}

void countButterfliesSortCE(uintE* butterflies, long nv, long nu, UVertexPair* wedges, uintE* wedge_freqs_f, 
            long num_wedge_freqs_f, bool use_v1) {
  // Given this list wedge_freqs_f, retrieve a list of indices for wedge_freqs_f such that
  // consecutive wedges have different first vertices
  // That is to say, the list we construct par_idxs_f is a list of indices on wedge_freqs_f,
  // such that wedges[wedge_freqs_f[par_idxs_f[i]]] and wedges[wedge_freqs_f[par_idxs_f[i+1]]]
  // have different first vertices
  pair<uintE*, long> par_idxs_pair = getFreqs(wedge_freqs_f,num_wedge_freqs_f-1, 
     NestedUVPCmp(wedges,use_v1), NestedUVPEq(wedges,use_v1), false);
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
  free(par_idxs_f);
}

void countButterfliesSort(uintE* butterflies, UVertexPair* wedges, uintE* wedge_freqs_f, long num_wedge_freqs_f) {
  // Then, retrieve the frequency counts for each distinct key by taking the difference between
  // these indices
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // Store butterfly counts with each vertex in U
  parallel_for (long i = 0; i < num_wedge_freqs_f-1; ++i) {
    uintE num_butterflies = wedge_freqs_f[i+1] - wedge_freqs_f[i];
    long wedge_idx = wedge_freqs_f[i];

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    writeAdd(&butterflies[wedges[wedge_idx].v1],num_butterflies); 
    writeAdd(&butterflies[wedges[wedge_idx].v2],num_butterflies); 
  }
}

// This is the original compute function, without the more cache-efficient sorting method
template <class vertex>
long CountSort(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE** nbhd_idxs, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  pair<pair<UVertexPair*,long>, long> wedges_pair  = getWedges<UVertexPair>(nu, V, U, UVertexPairCons(), max_wedges, curr_idx, num_wedges,
    wedge_idxs, nbhd_idxs);

  UVertexPair* wedges = wedges_pair.first.first;
  long num_wedges_curr = wedges_pair.first.second;


  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges_curr, true);

  countButterfliesSort(butterflies, wedges, freq_pair.first, freq_pair.second);

  free(freq_pair.first);

  free(wedges);
  return wedges_pair.second;
}

// This is the new compute function, with cache efficient sorting
template <class vertex>
long CountSortCE(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE** nbhd_idxs, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  pair<pair<UVertexPair*,long>, long> wedges_pair = getWedges<UVertexPair>(nu, V, U, UVertexPairCons(), max_wedges,
    curr_idx, num_wedges, wedge_idxs, nbhd_idxs);
  UVertexPair* wedges = wedges_pair.first.first;
  long num_wedges_curr = wedges_pair.first.second;

  // Retrieve butterflies + store with first vertex
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges_curr, true);

  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair.first, freq_pair.second , true);

  pair<uintE*, long> freq_pair2 = getWedgeFreqs(nv, nu, wedges, num_wedges_curr, false);

  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair2.first, freq_pair2.second , false);

  free(freq_pair.first);
  free(freq_pair2.first);
  free(wedges);

  return wedges_pair.second;

}

//********************************************************************************************
//********************************************************************************************

// TODO modularize this stuff better
template <class vertex>
long CountHash(sparseAdditiveSet<uintE> wedges, bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  using T = pair<uintE,uintE>;

  long next_idx = getWedgesHash(wedges, nu, V, U, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  
  _seq<T> wedge_freqs = wedges.entries();

  // Retrieve count on each key; that number choose 2 is the number of butterflies  
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
    T wedge_freq_pair = wedge_freqs.A[i];
    uintE num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;
    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    if (num_butterflies > 0){
      writeAdd(&butterflies[v1],num_butterflies);
      writeAdd(&butterflies[v2],num_butterflies);
    }
  }

  wedge_freqs.del();
  return next_idx;
}

template <class vertex>
long CountHashCE(sparseAdditiveSet<uintE>& wedges, bipartiteGraph<vertex>& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  using T = pair<uintE,uintE>;

  long next_idx = getWedgesHash(wedges, nu, V, U, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  getWedgesFromHashTimer.start();
  _seq<T> wedge_freqs = wedges.entries();
  getWedgesFromHashTimer.stop();
  numButterfliesHashInsertTimer.start();

  //sparseAdditiveSet<uintE> butterflies_set = sparseAdditiveSet<uintE>(nu,1,UINT_E_MAX);
  wedges.clear();


  // Retrieve count on each key; that number choose 2 is the number of butterflies
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
    T wedge_freq_pair = wedge_freqs.A[i];
    uintE num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    if (num_butterflies > 0) {
      wedges.insert(T(v1, num_butterflies));
      wedges.insert(T(v2, num_butterflies));
    }
  }

  _seq<T> butterflies_seq = wedges.entries();

  parallel_for(long i=0; i < butterflies_seq.n; ++i) {
    T butterfly_pair = butterflies_seq.A[i];
    butterflies[butterfly_pair.first] += butterfly_pair.second;
  }
  numButterfliesHashInsertTimer.stop();  
  butterflies_seq.del();
  //butterflies_set.del();
  wedge_freqs.del();
  
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

template <class vertex>
long CountHist(pbbsa::sequence<tuple<uintE, uintE>> tmp, pbbsa::sequence<tuple<uintE, uintE>> out, bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE** nbhd_idxs, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;
  using T = tuple<uintE,uintE>;

  pair<pair<uintE*,long>, long> wedges_list_pair = getWedges<uintE>(nu, V, U, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs, nbhd_idxs);
  uintE* wedges_list = wedges_list_pair.first.first;
  long num_wedges_curr = wedges_list_pair.first.second;
  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges_curr);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.

  tuple<size_t,T*> wedges_tuple = pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nu*nu, tmp, out);
  T* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    T wedge_freq_pair = wedge_freqs[i];
    uintE num_butterflies = get<1>(wedge_freq_pair);
    uintE wedge_freq_pair_first = get<0>(wedge_freq_pair);

    uintE v1 = wedge_freq_pair_first / nu;
    uintE v2 = wedge_freq_pair_first % nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    if (num_butterflies > 0) {
      writeAdd(&butterflies[v1],num_butterflies);
      writeAdd(&butterflies[v2],num_butterflies);
    }
  }

  //free(wedge_freqs);
  free(wedges_list);
  return wedges_list_pair.second;
}

template <class vertex>
long CountHistCE(pbbsa::sequence<tuple<uintE, uintE>>& tmp, pbbsa::sequence<tuple<uintE, uintE>>& out, bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE** nbhd_idxs, long curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  using X = tuple<uintE,uintE>;

  pair<pair<uintE*,long>, long> wedges_list_pair = getWedges<uintE>(nu, V, U, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs, nbhd_idxs);
  uintE* wedges_list = wedges_list_pair.first.first;
  long num_wedges_curr = wedges_list_pair.first.second;
  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges_curr);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nv*nu + nu, tmp, out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  X* wedge_freqs_i = newA(X, 2*wedge_freqs_n);
  parallel_for(long i = 0; i < wedge_freqs_n; i++) {
    X wedge_freq_pair = wedge_freqs[i];
    uintE num = get<1>(wedge_freq_pair);
    num = (num * (num-1)) / 2;
    uintE wedge_num = get<0>(wedge_freq_pair);
    wedge_freqs_i[2*i] = make_tuple(wedge_num % nu, num);
    wedge_freqs_i[2*i + 1] = make_tuple(wedge_num / nu, num);
  }
  //free(wedge_freqs);
  free(wedges_list);

  pbbsa::sequence<X> wedge_freqs_i_seq = pbbsa::sequence<X>(wedge_freqs_i,2*wedge_freqs_n);
  tuple<size_t, X*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(wedge_freqs_i_seq,nu,getAdd<uintE,uintE>, getAddReduce<uintE,uintE>, tmp, out);
  X* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }
  
  free(wedge_freqs_i);
  //free(butterflies_l);

  return wedges_list_pair.second;
}

//********************************************************************************************
//********************************************************************************************

template <class E, class I>
struct nestA { 
  E* arr1;
  E* arr2;
  nestA(E* _arr1, E* _arr2) : arr1(_arr1), arr2(_arr2) {}
  inline E operator() (const I& i) const {
    return arr1[arr2[i]] * (arr1[arr2[i]]-1) / 2;
  }
};

template <class vertex>
void CountSeq_helper(long curr_idx, long next_idx, uintE* wedges, uintE* used, uintE* butterflies,
  const long nv, const long nu, const vertex* V, const vertex* U) {
  parallel_for_1(long i=curr_idx; i < next_idx; ++i){
    //seqLoopTimer.start();
    //parallel_for(long i=0; i < nu; ++i) { wedges[i] = 0; }
    //seqLoopTimer.stop();
    long used_idx=0;
    long shift = (i-curr_idx) * nu;
    vertex u = U[i];
    for (long j=0; j < u.getOutDegree(); ++j ) {
      vertex v = V[u.getOutNeighbor(j)];
      // Find all seagulls with center v and endpoint u
      for (long k=0; k < v.getOutDegree(); ++k) { 
        uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx < i) {
	  
          writeAdd(&butterflies[i*16], wedges[shift+u2_idx]);
          writeAdd(&butterflies[u2_idx*16], wedges[shift+u2_idx]);
          //writeAdd(&butterflies[8*(i-curr_idx)],wedges[shift + u2_idx]);
	  wedges[shift+u2_idx]++;
          if (wedges[shift + u2_idx] == 1) {used[shift + used_idx] = shift + u2_idx; used_idx++;}
        }
        else break;
      }
    }
    //seqWriteTimer.start();
    //uintE total = used_idx == 0 ? 0 : sequence::reduce<long>((long) shift, (long) shift + used_idx, addF<uintE>(), nestA<uintE, long>(wedges, used));
    //cout << (double)used_idx/nu << endl;
    //writeAdd(&butterflies[0], 2*total);
    /* parallel_for(long j=shift; j < shift + used_idx; ++j) { */
    /*   uintE num_butterflies = wedges[used[j]] * (wedges[used[j]] - 1) / 2; */
    /*   writeAdd(&butterflies[used[j] - shift], num_butterflies); */
    /* } */
    //seqWriteTimer.stop();
    parallel_for(long j=shift; j < shift + used_idx; ++j) { wedges[used[j]] = 0; }
  }
}

template <class vertex>
void CountSeq(bipartiteGraph<vertex> GA, bool use_v, uintE* butterflies, long max_wedges) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  long step = 1000;//max_wedges / nu;
  uintE* wedges = newA(uintE, step * nu);
  uintE* used = newA(uintE, step * nu);
  parallel_for(long i=0; i < step * nu; ++i) { wedges[i] = 0; }
  for(long i=0; i < nu; i += step) {
    CountSeq_helper(i, i+step < nu ? i+step : nu, wedges, used, butterflies, nv, nu, V, U);
  }

  free(wedges);
  free(used);
}

template <class vertex>
void CountOrig(bipartiteGraph<vertex> GA, bool use_v, uintE* butterflies) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  long results = 0;

  uintE* wedges = newA(uintE, nu);
  uintE* used = newA(uintE, nu);
  for(long i=0; i < nu; ++i) { wedges[i] = 0; }

  for(long i=0; i < nu; ++i){
    long used_idx = 0;
    vertex u = U[i];
    for (long j=0; j < U[i].getOutDegree(); ++j ) {
      vertex v = V[U[i].getOutNeighbor(j)];
      for (long k=0; k < v.getOutDegree(); ++k) { 
        uintE u2_idx = v.getOutNeighbor(k);
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
    for(long j=0; j < used_idx; ++j) { wedges[used[j]] = 0; }
  }
  
  free(wedges);
  free(used);
  cout << "num: " << results << "\n";
}

//********************************************************************************************
//********************************************************************************************

template <class vertex>
void CountHash_helper(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  initHashTimer.start();

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(nu,1,UINT_E_MAX);
  // TODO fix hist so that it saves storage
  // TODO also can config float factor for hash table
  initHashTimer.stop();

  if (max_wedges >= num_wedges) {
    if (type == 2) CountHash(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else CountHashCE(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    wedges.del();
    return;
  }
  long curr_idx = 0;
  while(curr_idx < nu) {
    if (type ==2) curr_idx = CountHash(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else curr_idx = CountHashCE(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    wedges.clear();
  }
  wedges.del();

  nextWedgeTimer.reportTotal("getting next wedge");
  initHashTimer.reportTotal("init first hash table");
  hashInsertTimer.reportTotal("inserting into first hash table");
  getWedgesFromHashTimer.reportTotal("extract from first hash table");
  numButterfliesHashInsertTimer.reportTotal("init and inserting into second hash table");
}

template <class vertex>
void CountHist_helper(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, uintE** nbhd_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  pbbsa::sequence<tuple<uintE, uintE>> tmp = pbbsa::sequence<tuple<uintE, uintE>>();
  pbbsa::sequence<tuple<uintE, uintE>> out = pbbsa::sequence<tuple<uintE, uintE>>();

  if (max_wedges >= num_wedges) {
    if (type == 4) CountHist(tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs);
    else CountHistCE(tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs);
    return;
  }
  long curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 4) curr_idx = CountHist(tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs, curr_idx);
    else curr_idx = CountHistCE(tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs, curr_idx);
  }
}

template <class vertex>
uintE* Count(bipartiteGraph<vertex> GA, bool use_v, long num_wedges, long max_wedges, long type=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  uintE* butterflies = newA(uintE, 16*nu);
  parallel_for(long i=0;i<nu*16;++i){
    butterflies[i] = 0;
  }

  uintE* wedge_idxs = nullptr;
  uintE** nbhd_idxs = nullptr;

  pair<uintE*, uintE**> idxs_pair = countWedgesScan(nu, V, U, true, type != 2 && type != 3);
  wedge_idxs = idxs_pair.first;
  nbhd_idxs = idxs_pair.second;

  if (type > 8 || type == 5) return butterflies;

  if (type == 2 || type == 3) {
    CountHash_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
    return butterflies;
  }
  else if (type == 4 || type == 6) {
    CountHist_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs, type);
    return butterflies;
  }
  else if (type == 7) {
    CountSeq(GA, use_v, butterflies, max_wedges);
    seqWriteTimer.reportTotal("seq write time");
    seqLoopTimer.reportTotal("seq allocate space time");
    return butterflies;
  }
  else if (type ==8) {
    CountOrig(GA, use_v, butterflies);
    return butterflies;
  }

  if (max_wedges >= num_wedges) {
    if (type == 0) CountSort(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs);
    else if (type == 1) CountSortCE(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs);
    return butterflies;
  }
  long curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 0) curr_idx = CountSort(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs, curr_idx);
    else if (type == 1) curr_idx = CountSortCE(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, nbhd_idxs, curr_idx);
  }

  if (wedge_idxs) free(wedge_idxs);
  if (nbhd_idxs) free(nbhd_idxs);
  return butterflies;

}

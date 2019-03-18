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
  const intT eltsPerCacheLine = 64/sizeof(long);
  for (long idx = par_idxs_f[i]; idx < par_idxs_f[i+1]; ++idx) { //TODO this is without reduce -- do thresholding
      uintE num_butterflies = wedge_freqs_f[idx+1] - wedge_freqs_f[idx];
      long wedge_idx = wedge_freqs_f[idx];
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      
      UVertexPair vs = wedges[wedge_idx];

      if (use_v1) butterflies[eltsPerCacheLine*vs.v1] += num_butterflies;
      else butterflies[eltsPerCacheLine*vs.v2] += num_butterflies;
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
  const intT eltsPerCacheLine = 64/sizeof(long);
  uintE* butterflies_v = newA(uintE, nu);
  parallel_for(intT i=0; i < nu; ++i) {butterflies_v[i] = 0; }
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
  if (use_v1) butterflies[eltsPerCacheLine*start_vs.v1] += sum;
  else butterflies[eltsPerCacheLine*start_vs.v2] += sum;

  free(butterflies_v);
}

// Retrieve frequency counts for all wedges with the same key, and for all wedges with the same first vertex
// (if use_v1; otherwise, second vertex -- in comments, we assume first vertex)
// First, retrieve a list of indices where consecutive wedges have different keys
pair<uintE*, long> getWedgeFreqs(long nv, long nu, UVertexPair* wedges, long num_wedges, bool use_v1) {
    if (use_v1) return getFreqs(wedges, num_wedges, UVertexPairCmp(), UVertexPairEq());
    return getFreqs(wedges, num_wedges, UVertexPairCmp2(), UVertexPairEq());
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
  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for (long i = 0; i < num_wedge_freqs_f-1; ++i) {
    uintE num_butterflies = wedge_freqs_f[i+1] - wedge_freqs_f[i];
    long wedge_idx = wedge_freqs_f[i];

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v1],num_butterflies); 
    writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v2],num_butterflies); 
  }
}

// This is the original compute function, without the more cache-efficient sorting method
intT CountSort(_seq<UVertexPair>& wedges_seq, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  pair<long, intT> wedges_pair  = getWedges<UVertexPair>(wedges_seq, GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges,
    wedge_idxs);

  UVertexPair* wedges = wedges_seq.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges_curr, true);

  countButterfliesSort(butterflies, wedges, freq_pair.first, freq_pair.second);

  free(freq_pair.first);
  //free(wedges);
  return wedges_pair.second;
}

// This is the new compute function, with cache efficient sorting
/*intT CountSortCE(_seq<UVertexPair>& wedges_seq, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  pair<long, intT> wedges_pair = getWedges<UVertexPair>(wedges_seq, GA, use_v, UVertexPairCons(), max_wedges,
    curr_idx, num_wedges, wedge_idxs);
  UVertexPair* wedges = wedges_seq.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve butterflies + store with first vertex
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges_curr, true);

  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair.first, freq_pair.second , true);

  pair<uintE*, long> freq_pair2 = getWedgeFreqs(nv, nu, wedges, num_wedges_curr, false);

  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair2.first, freq_pair2.second , false);

  free(freq_pair.first);
  free(freq_pair2.first);
  //free(wedges);

  return wedges_pair.second;

}*/

intT CountSortCE(_seq<UVertexPair>& wedges_seq, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UVertexPair>(wedges_seq, GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UVertexPair* wedges = wedges_seq.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges_f, UVertexPairCmp(), UVertexPairEq());
  uintE* freq_arr = freq_pair.first;

  using X = tuple<uintE,uintE>;
  X* b_freqs = newA(X, 2*freq_pair.second);

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i];
    num = num * (num-1) / 2;
    //parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
    long j = freq_arr[i];
    b_freqs[2*i] = make_tuple(wedges[j].v1, num);
    b_freqs[2*i+1] = make_tuple(wedges[j].v2, num);
    //}
  }

  free(freq_arr);

  // now, we need to collate by our indices
  pair<uintE*, long> b_freq_pair = getFreqs(b_freqs, 2*freq_pair.second, uintETupleLt(), uintETupleEq());
  uintE* b_freq_arr = b_freq_pair.first;
  const intT eltsPerCacheLine = 64/sizeof(long);

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    uintE num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(b_freqs[b_freq_arr[i-1]]), num_freq, uintETupleAdd());
    // These are our butterfly counts
    butterflies[get<0>(reduce)*eltsPerCacheLine] += get<1>(reduce);
  }

  free(b_freq_arr);
  free(b_freqs);

  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

// TODO modularize this stuff better
intT CountHash(sparseAdditiveSet<uintE>& wedges, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  using T = pair<uintE,uintE>;

  intT next_idx = getWedgesHash(wedges, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  
  _seq<T> wedge_freqs = wedges.entries();

  // Retrieve count on each key; that number choose 2 is the number of butterflies  
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
    T wedge_freq_pair = wedge_freqs.A[i];
    uintE num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;
    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    if (num_butterflies > 0){
      writeAdd(&butterflies[eltsPerCacheLine*v1],num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*v2],num_butterflies);
    }
  }

  wedge_freqs.del();
  return next_idx;
}

intT CountHashCE(sparseAdditiveSet<uintE>& wedges, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  using T = pair<uintE,uintE>;

  intT next_idx = getWedgesHash(wedges, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  //getWedgesFromHashTimer.start();
  _seq<T> wedge_freqs = wedges.entries();
  //getWedgesFromHashTimer.stop();
  //numButterfliesHashInsertTimer.start();

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

   wedge_freqs.del();
  _seq<T> butterflies_seq = wedges.entries();

  parallel_for(long i=0; i < butterflies_seq.n; ++i) {
    T butterfly_pair = butterflies_seq.A[i];
    butterflies[eltsPerCacheLine*butterfly_pair.first] += butterfly_pair.second;
  }
  //numButterfliesHashInsertTimer.stop();  
  butterflies_seq.del();
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

intT CountHist(_seq<uintE>& wedges, pbbsa::sequence<tuple<uintE, uintE>> tmp, pbbsa::sequence<tuple<uintE, uintE>> out,
  bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  using T = tuple<uintE,uintE>;

  pair<long, intT> wedges_list_pair = getWedges<uintE>(wedges, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  uintE* wedges_list = wedges.A;
  long num_wedges_curr = wedges_list_pair.first;
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
      writeAdd(&butterflies[eltsPerCacheLine*v1],num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*v2],num_butterflies);
    }
  }

  return wedges_list_pair.second;
}

//TODO can reuse more space
intT CountHistCE(_seq<uintE>& wedges, pbbsa::sequence<tuple<uintE, uintE>>& tmp, pbbsa::sequence<tuple<uintE, uintE>>& out, 
  bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  const long nv = use_v ? GA.nv : GA.nu;
  const intT eltsPerCacheLine = 64/sizeof(long);

  using X = tuple<uintE,uintE>;

  pair<long, intT> wedges_list_pair = getWedges<uintE>(wedges, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  uintE* wedges_list = wedges.A;
  long num_wedges_curr = wedges_list_pair.first;
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

  pbbsa::sequence<X> wedge_freqs_i_seq = pbbsa::sequence<X>(wedge_freqs_i,2*wedge_freqs_n);
  tuple<size_t, X*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(wedge_freqs_i_seq,nu,getAdd<uintE,uintE>, getAddReduce<uintE,uintE>, tmp, out);
  X* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }
  
  free(wedge_freqs_i);

  return wedges_list_pair.second;
}

//********************************************************************************************
//********************************************************************************************

void CountHash_helper(bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  //initHashTimer.start();

  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(nu,1,UINT_E_MAX);
  // TODO fix hist so that it saves storage
  // TODO also can config float factor for hash table
  //initHashTimer.stop();

  if (max_wedges >= num_wedges) {
    if (type == 2) CountHash(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else CountHashCE(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    wedges.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type ==2) curr_idx = CountHash(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else curr_idx = CountHashCE(wedges, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    wedges.clear();
  }
  wedges.del();

  /*nextWedgeTimer.reportTotal("getting next wedge");
  initHashTimer.reportTotal("init first hash table");
  hashInsertTimer.reportTotal("inserting into first hash table");
  getWedgesFromHashTimer.reportTotal("extract from first hash table");
  numButterfliesHashInsertTimer.reportTotal("init and inserting into second hash table");*/
}

void CountHist_helper(bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  pbbsa::sequence<tuple<uintE, uintE>> tmp = pbbsa::sequence<tuple<uintE, uintE>>();
  pbbsa::sequence<tuple<uintE, uintE>> out = pbbsa::sequence<tuple<uintE, uintE>>();
  _seq<uintE> wedges_seq = _seq<uintE>(newA(uintE, nu), nu);

  if (max_wedges >= num_wedges) {
    if (type == 4) CountHist(wedges_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else CountHistCE(wedges_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    wedges_seq.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 4) curr_idx = CountHist(wedges_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else curr_idx = CountHistCE(wedges_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
  }
  wedges_seq.del();
}

void CountSort_helper(bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  _seq<UVertexPair> wedges_seq = _seq<UVertexPair>(newA(UVertexPair, nu), nu);

  if (max_wedges >= num_wedges) {
    if (type == 0) CountSort(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else CountSortCE(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    wedges_seq.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 0) curr_idx = CountSort(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else curr_idx = CountSortCE(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
  }
  wedges_seq.del();
}

uintE* Count(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long type=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  //TODO put in thing where we don't have contention b./c of cache line
  const intT eltsPerCacheLine = 64/sizeof(long);
  uintE* butterflies = newA(uintE, eltsPerCacheLine*nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i*eltsPerCacheLine] = 0;
  }

  uintE* wedge_idxs = countWedgesScan(GA, use_v, true);

  // TODO make enums?
  // TODO make CountSorthelper so that this is more streamlined
  if (type == 2 || type == 3) CountHash_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
  else if (type == 4 || type == 6) CountHist_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
  else if (type == 0 || type == 1) CountSort_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);

  free(wedge_idxs);
  return butterflies;

}

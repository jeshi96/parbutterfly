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

intT CountSortCE(_seq<UVertexPair>& wedges_seq, _seq<tuple<uintE, uintE>>& butterflies_seq, bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UVertexPair>(wedges_seq, GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UVertexPair* wedges = wedges_seq.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  pair<uintE*, long> freq_pair = getFreqs(wedges, num_wedges_f, UVertexPairCmp(), UVertexPairEq());
  uintE* freq_arr = freq_pair.first;

  using X = tuple<uintE,uintE>;
  if (butterflies_seq.n < 2*freq_pair.second) {
    free(butterflies_seq.A);
    butterflies_seq.A = newA(X, 2*freq_pair.second);
    butterflies_seq.n = 2*freq_pair.second;
  }

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    uintE num = freq_arr[i+1] - freq_arr[i];
    num = num * (num-1) / 2;
    //parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
    long j = freq_arr[i];
    butterflies_seq.A[2*i] = make_tuple(wedges[j].v1, num);
    butterflies_seq.A[2*i+1] = make_tuple(wedges[j].v2, num);
    //}
  }

  free(freq_arr);

  // now, we need to collate by our indices
  pair<uintE*, long> b_freq_pair = getFreqs(butterflies_seq.A, 2*freq_pair.second, uintETupleLt(), uintETupleEq());
  uintE* b_freq_arr = b_freq_pair.first;
  const intT eltsPerCacheLine = 64/sizeof(long);

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    uintE num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(butterflies_seq.A[b_freq_arr[i-1]]), num_freq, uintETupleAdd());
    // These are our butterfly counts
    butterflies[get<0>(reduce)*eltsPerCacheLine] += get<1>(reduce);
  }

  free(b_freq_arr);

  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

// TODO modularize this stuff better
intT CountHash(sparseAdditiveSet<uintE>& wedges, _seq<pair<uintE,uintE>>& wedges_seq, bipartiteCSR& GA, 
  bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  using T = pair<uintE,uintE>;

  intT next_idx = getWedgesHash(wedges, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  
  intT num_wedges_seq = wedges.entries_no_init(wedges_seq);

  // Retrieve count on each key; that number choose 2 is the number of butterflies  
  parallel_for (long i=0; i < num_wedges_seq; ++i) {
    T wedge_freq_pair = wedges_seq.A[i];
    uintE num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;
    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    if (num_butterflies > 0){
      writeAdd(&butterflies[eltsPerCacheLine*v1],num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*v2],num_butterflies);
    }
  }

  return next_idx;
}

intT CountHashCE(sparseAdditiveSet<uintE>& wedges, _seq<pair<uintE,uintE>>& wedges_seq, bipartiteCSR& GA, 
  bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  using T = pair<uintE,uintE>;

  intT next_idx = getWedgesHash(wedges, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  //getWedgesFromHashTimer.start();
  intT num_wedges_seq = wedges.entries_no_init(wedges_seq);
  //getWedgesFromHashTimer.stop();
  //numButterfliesHashInsertTimer.start();

  wedges.clear();
  wedges.resize(nu);

  // Retrieve count on each key; that number choose 2 is the number of butterflies
  parallel_for (long i=0; i < num_wedges_seq; ++i) {
    T wedge_freq_pair = wedges_seq.A[i];
    uintE num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    if (num_butterflies > 0) {
      wedges.insert(T(v1, num_butterflies));
      wedges.insert(T(v2, num_butterflies));
    }
  }
  num_wedges_seq = wedges.entries_no_init(wedges_seq);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    T butterfly_pair = wedges_seq.A[i];
    butterflies[eltsPerCacheLine*butterfly_pair.first] += butterfly_pair.second;
  }
  //numButterfliesHashInsertTimer.stop();  
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
intT CountHistCE(_seq<uintE>& wedges, _seq<tuple<uintE, uintE>>& butterflies_seq, pbbsa::sequence<tuple<uintE, uintE>>& tmp, pbbsa::sequence<tuple<uintE, uintE>>& out, 
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

  if (butterflies_seq.n < 2*wedge_freqs_n) {
    free(butterflies_seq.A);
    butterflies_seq.A = newA(X, 2*wedge_freqs_n);
    butterflies_seq.n = 2*wedge_freqs_n;
  }

  parallel_for(long i = 0; i < wedge_freqs_n; i++) {
    X wedge_freq_pair = wedge_freqs[i];
    uintE num = get<1>(wedge_freq_pair);
    num = (num * (num-1)) / 2;
    uintE wedge_num = get<0>(wedge_freq_pair);
    butterflies_seq.A[2*i] = make_tuple(wedge_num % nu, num);
    butterflies_seq.A[2*i + 1] = make_tuple(wedge_num / nu, num);
  }

  pbbsa::sequence<X> wedge_freqs_i_seq = pbbsa::sequence<X>(butterflies_seq.A,2*wedge_freqs_n);
  tuple<size_t, X*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(wedge_freqs_i_seq,nu,getAdd<uintE,uintE>, getAddReduce<uintE,uintE>, tmp, out);
  X* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  return wedges_list_pair.second;
}

//********************************************************************************************
//********************************************************************************************

void CountHash_helper(bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  //initHashTimer.start();

  using T = pair<uintE,uintE>;
  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(nu,1,UINT_E_MAX);
  _seq<T> wedges_seq = _seq<T>(newA(T, nu), nu);
  // TODO fix hist so that it saves storage
  // TODO also can config float factor for hash table
  //initHashTimer.stop();

  if (max_wedges >= num_wedges) {
    if (type == 2) CountHash(wedges, wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else CountHashCE(wedges, wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    wedges.del();
    wedges_seq.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type ==2) curr_idx = CountHash(wedges, wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else curr_idx = CountHashCE(wedges, wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    wedges.clear();
  }
  wedges.del();
  wedges_seq.del();

  /*nextWedgeTimer.reportTotal("getting next wedge");
  initHashTimer.reportTotal("init first hash table");
  hashInsertTimer.reportTotal("inserting into first hash table");
  getWedgesFromHashTimer.reportTotal("extract from first hash table");
  numButterfliesHashInsertTimer.reportTotal("init and inserting into second hash table");*/
}

void CountHist_helper(bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  using T = tuple<uintE,uintE>;
  pbbsa::sequence<tuple<uintE, uintE>> tmp = pbbsa::sequence<tuple<uintE, uintE>>();
  pbbsa::sequence<tuple<uintE, uintE>> out = pbbsa::sequence<tuple<uintE, uintE>>();
  _seq<uintE> wedges_seq = _seq<uintE>(newA(uintE, nu), nu);
  _seq<T> butterflies_seq = _seq<T>(newA(T, 1), 1);

  if (max_wedges >= num_wedges) {
    if (type == 4) CountHist(wedges_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else CountHistCE(wedges_seq, butterflies_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    wedges_seq.del();
    butterflies_seq.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 4) curr_idx = CountHist(wedges_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else curr_idx = CountHistCE(wedges_seq, butterflies_seq, tmp, out, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
  }
  wedges_seq.del();
  butterflies_seq.del();
}

void CountSort_helper(bipartiteCSR& GA, bool use_v, long num_wedges, uintE* butterflies, long max_wedges, uintE* wedge_idxs, long type) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  using T = tuple<uintE, uintE>;
  _seq<UVertexPair> wedges_seq = _seq<UVertexPair>(newA(UVertexPair, nu), nu);
  _seq<T> butterflies_seq = _seq<T>(newA(T, 1), 1);

  if (max_wedges >= num_wedges) {
    if (type == 0) CountSort(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else CountSortCE(wedges_seq, butterflies_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    wedges_seq.del();
    butterflies_seq.del();
    return;
  }
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 0) curr_idx = CountSort(wedges_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else curr_idx = CountSortCE(wedges_seq, butterflies_seq, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
  }
  wedges_seq.del();
  butterflies_seq.del();
}

//********************************************************************************************
//********************************************************************************************


void CountOrigCompactParallel(bipartiteCSR& GA, bool use_v, uintE* butterflies) {
  timer t1,t2;
  t1.start();
  cout << "Original Parallel" << endl;
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = 1000; //tunable parameter

  long* wedges = newA(long, nu*stepSize);
  long* used = newA(long, nu*stepSize);

  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);
  
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
	    //butterflies[i*eltsPerCacheLine] += wedges[shift+u2_idx];
	    //butterflies[u2_idx*eltsPerCacheLine] += wedges[shift+u2_idx];
	    //results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	  }
	  else break;
	}
      }
      for(intT j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j];
        writeAdd(&butterflies[i*eltsPerCacheLine],  (uintE)((uintE) wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
        writeAdd(&butterflies[u2_idx*eltsPerCacheLine], (uintE)((uintE) wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
        wedges[u2_idx+shift] = 0;
      }
    }
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);
}

uintE* CountWorkEfficientParallel(graphCSR& GA) {
  timer t1,t2;
  t1.start();

  long stepSize = getWorkers() * 15; //tunable parameter
  uintE* wedges = newA(uintE, GA.n*stepSize);
  uintE* used = newA(uintE, GA.n*stepSize);

  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);
  
  uintE* butterflies = newA(uintE,eltsPerCacheLine*GA.n);
  granular_for(i,0,GA.n,GA.n > 10000, { butterflies[eltsPerCacheLine*i] = 0; });
  t1.reportTotal("preprocess");

  t2.start();

  for(long step = 0; step < (GA.n+stepSize-1)/stepSize; step++) {
    parallel_for_1(long i=step*stepSize; i < min((step+1)*stepSize,GA.n); ++i){
      intT used_idx = 0;
      long shift = GA.n*(i-step*stepSize);
      intT u_offset  = GA.offsets[i];
      intT u_deg = GA.offsets[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	intT v = GA.edges[u_offset+j] >> 1;
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

      for (long j=0; j < u_deg; ++j ) {
	intT v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
	if (v <= i) break;
  if (!(GA.edges[u_offset+j] & 0b1)) continue;
	for (long k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
	  if (u2_idx > i) { //TODO combine into one graph
	    if (wedges[shift+u2_idx] > 1) writeAdd(&butterflies[eltsPerCacheLine*v], (uintE)(wedges[shift+u2_idx]-1));
	  }
	  else break;
	}
      }

      for(long j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j] >> 1;
        if(used[shift+j] & 0b1) {
        writeAdd(&butterflies[i*eltsPerCacheLine],  (uintE)((uintE)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
        writeAdd(&butterflies[u2_idx*eltsPerCacheLine], (uintE)((uintE)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
        }
        wedges[shift+u2_idx] = 0;
      }
    }
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);

  return butterflies;
}

uintE* Count(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long type=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  //TODO put in thing where we don't have contention b./c of cache line
  const intT eltsPerCacheLine = 64/sizeof(long);
  uintE* butterflies = newA(uintE, eltsPerCacheLine*nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[eltsPerCacheLine * i] = 0; });

  uintE* wedge_idxs = countWedgesScan(GA, use_v, true);

  // TODO make enums?
  // TODO make CountSorthelper so that this is more streamlined
  if (type == 2 || type == 3) CountHash_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
  else if (type == 4 || type == 6) CountHist_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
  else if (type == 0 || type == 1) CountSort_helper(GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, type);
  else if (type == 7) CountOrigCompactParallel(GA, use_v, butterflies);
  else if (type == 11) {
    timer t_rank;
    t_rank.start();
    auto rank_tup = getDegRanks(GA);
    auto rank_tup2 = getCoreRanks(GA);
    auto rank_tup3 = getCoCoreRanks(GA);
    t_rank.reportTotal("ranking");

    long num_rwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
      rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup), get<1>(rank_tup)));
    num_rwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
      rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup)));

    long num_cwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
      rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup2), get<1>(rank_tup2)));
    num_cwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
      rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup2), get<2>(rank_tup2)));

    long num_ccwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
      rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup3), get<1>(rank_tup3)));
    num_ccwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
      rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup3), get<2>(rank_tup3)));

    cout << "Rank wedges: " << num_rwedges << "\n";
    cout << "Side wedges: " << num_wedges << "\n";
    cout << "Core wedges: " << num_cwedges << "\n";
    cout << "Co Core wedges: " << num_ccwedges << "\n";
    
    if (num_rwedges < num_wedges + 1000) {
      auto g = rankGraph(GA, use_v, get<0>(rank_tup), get<1>(rank_tup), get<2>(rank_tup));
      free(get<0>(rank_tup));

      uintE* rank_butterflies = CountWorkEfficientParallel(g);
      auto rank_converter = use_v ? get<2>(rank_tup) : get<1>(rank_tup);
      granular_for(i,0,nu,nu > 10000, { 
        butterflies[eltsPerCacheLine*i] = rank_butterflies[eltsPerCacheLine*rank_converter[i]];
      });
      free(get<1>(rank_tup)); free(get<2>(rank_tup)); free(rank_butterflies);
    }
    else  {
      free(get<0>(rank_tup)); free(get<1>(rank_tup)); free(get<2>(rank_tup));
      CountOrigCompactParallel(GA, use_v, butterflies);
    }
    
  }

  free(wedge_idxs);
  return butterflies;

}

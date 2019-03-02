#ifndef _BUTILS_
#define _BUTILS_

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "hypergraphIO.h"
#include "parseCommandLine.h"
#include "vertex.h"
#include "sequence.h"
#include "binary_search.h"

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }


using namespace std;

template <class E>
struct seagullSumHelper { 
  uintE u;
  uintT* offsetsU;
  uintT* offsetsV;
  uintE* edgesU;
  seagullSumHelper(uintE _u, uintT* _offsetsU, uintT* _offsetsV, uintE* _edgesU) : u(_u), offsetsU(_offsetsU), offsetsV(_offsetsV), edgesU(_edgesU) {}
  inline E operator() (const E& i) const {
    intT u_offset = offsetsU[u];
    uintE v = edgesU[u_offset + i];
	  return (E) (offsetsV[v+1] - offsetsV[v] - 1);
  }
};

template <class E>
struct seagullSum { 
  uintT* offsetsU;
  uintT* offsetsV;
  uintE* edgesU;
  uintE* active;
  seagullSum(uintT* _offsetsU, uintT* _offsetsV, uintE* _edgesU, uintE* _active) : offsetsU(_offsetsU), offsetsV(_offsetsV), edgesU(_edgesU), active(_active) {}
  inline E operator() (const E& i) const {
    /*uintE u = active[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[active[i]+1] - offsetsU[active[i]];
    E ret=0;
    for (long k=0; k < u_deg; ++k) {
      uintE v = edgesU[u_offset + k];
      ret += (offsetsV[v+1] - offsetsV[v] - 1);
    }
  return ret;*/
    intT u_deg = offsetsU[active[i]+1] - offsetsU[active[i]];
	return sequence::reduce<E>((E) 0, u_deg, addF<E>(), seagullSumHelper<E>(active[i], offsetsU, offsetsV, edgesU));
  }
};

struct UWedgeIntCons {
  long nu;
  UWedgeIntCons(long _nu) : nu(_nu) {}
  inline tuple<uintE,uintE> operator() (uintE v1, uintE v2, uintE c) {
    return make_tuple(v1 * nu + v2, c);
  }
};

struct UWedge {
  uintE v1;
  uintE v2;
  uintE u;
  UWedge(uintE _v1, uintE _v2, uintE _u) : v1(_v1), v2(_v2), u(_u) {}
};

struct UWedgeCons { inline UWedge operator() (uintE v1, uintE v2, uintE c) { return UWedge(v1, v2, c); }};

struct UWedgeCmp {
  inline bool operator() (UWedge vs1, UWedge vs2) {
  	if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
  	return vs1.v1 < vs2.v1;
  }
};

struct UWedgeEq { inline bool operator() (UWedge vs1, UWedge vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Represents a pair of vertices on one side of a bipartite graph (ordered, with least vertex first)

// Represents a pair of vertices on one side of a bipartite graph (unordered, stored based on constructor order)
struct UVertexPair {
  uintE v1;
  uintE v2;
  UVertexPair(uintE _v1, uintE _v2) : v1(_v1), v2(_v2) {}
};

//TODO get rid of legacy _nv
// Comparer for VertexPair based on least vertex in pair and then greatest vertex in pair

// Comparer for VertexPair based on greatest vertex in pair and then least vertex in pair

// Comparer for UVertexPair
struct UVertexPairCmp2{
  long nv;
  UVertexPairCmp2(long _nv) : nv(_nv) {}
  inline bool operator() (UVertexPair vs1, UVertexPair vs2) {
    if (vs1.v2 == vs2.v2) return vs1.v1 < vs2.v1;
    return vs1.v2 < vs2.v2;
  }
};

struct UVertexPairCmp {
  long nv;
  UVertexPairCmp(long _nv) : nv(_nv) {}
  inline bool operator() (UVertexPair vs1, UVertexPair vs2) {
    if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
    return vs1.v1 < vs2.v1;
  }
};

// Equality for VertexPair and UVertexPair
struct UVertexPairEq { inline bool operator() (UVertexPair vs1, UVertexPair vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Constructs a VertexPair and UVertexPair
struct UVertexPairCons { inline UVertexPair operator() (uintE v1, uintE v2, uintE c) { return UVertexPair(v1, v2); }};

// Constructs a uintE form of a VertexPair and UVertexPair
struct UVertexPairIntCons {
  long nu;
  UVertexPairIntCons(long _nu) : nu(_nu) {}
  inline uintE operator() (uintE v1, uintE v2, uintE c) {
    return v1 * nu + v2;
  }
};

// Comparer for indices for VertexPairs in nest, by v1 or v2 (based on use_v1)
struct NestedUVPCmp {
  UVertexPair* nest;
  bool use_v1;
  NestedUVPCmp(UVertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  inline bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 < nest[idx2].v1;
    return nest[idx1].v2 < nest[idx2].v2;
  }
};
struct NestedUVPEq {
  UVertexPair* nest;
  bool use_v1;
  NestedUVPEq(UVertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  inline bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 == nest[idx2].v1;
    return nest[idx1].v2 == nest[idx2].v2;
  }
};

struct nonZeroF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != 0);}};
struct nonZeroPairF{inline bool operator() (pair<uintE,uintE> &a) {return (a.second != 0);}};
struct greaterOneF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) > 1);}};
template<class T> struct cmpF{inline bool operator() (T a, T b) {return a < b;}};

struct nonEmptyUVPF{inline bool operator() (UVertexPair &a) {return (a.v1 != UINT_E_MAX || a.v2 != UINT_E_MAX);}};
struct nonEmptyUWF{inline bool operator() (UWedge &a) {return (a.v1 != UINT_E_MAX || a.v2 != UINT_E_MAX || a.u != UINT_E_MAX);}};

struct nonMaxTupleF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != UINT_E_MAX || get<0>(a) != UINT_E_MAX);}};

struct uintELt {inline bool operator () (uintE a, uintE b) {return a < b;};};
struct uintEEq {inline bool operator() (uintE a, uintE b) {return a == b;};};
struct uintETupleLt {inline bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {return get<0>(a) < get<0>(b);} };
struct uintETupleEq {inline bool operator() (tuple<uintE,uintE> a,tuple<uintE,uintE> b) {return get<0>(a) == get<0>(b); }};
struct uintETupleAdd {
  inline tuple<uintE,uintE> operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) const {
    return make_tuple(get<0>(a), get<1>(a) + get<1>(b));
  };
};
struct uintETupleLtBoth {
  inline bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {
    if (get<0>(a) == get<0>(b)) return get<1>(a) < get<1>(b);
      return get<0>(a) < get<0>(b);
  }
};

template<class T> struct refl{inline T operator() (T obj) {return obj;}};
template<class T> struct reflCount{inline long operator() (T obj) {return 1;}};
struct UVertexPairV2{inline uintE operator() (UVertexPair obj) {return obj.v2;}};
struct choose2{inline uintE operator() (uintE obj) {return obj*(obj-1)/2;}};
struct uintETupleGet0{inline uintE operator() (tuple<uintE,uintE> obj) {return get<0>(obj);}};
struct uintECount{inline long operator() (tuple<uintE,uintE> obj) {return get<1>(obj);}};

template <class E>
struct writeAddArr {
  E* arr;
  writeAddArr(E* _arr) : arr(_arr) {}
  inline void operator() (long idx, E num) {
    writeAdd(&arr[idx], num);
  }
};

template <class E>
struct writeAddSet {
  sparseAdditiveSet<E> set;
  writeAddSet(sparseAdditiveSet<E> _set) : set(_set) {}
  inline void operator() (long idx, E num) {
    set.insert(pair<uintE, E>(idx, num));
  }
};

template<class E, class K>
E getAdd (E curr, tuple<K,E> v) {
  return curr + get<1>(v);
}

template<class E, class K>
tuple<K,E> getAddReduce (tuple<K,E> curr, tuple<K,E> v) {
  return make_tuple(get<0>(curr),get<1>(curr) + get<1>(v));
}

template<class E, class K>
pair<K,E> getAddReducePair (pair<K,E> curr, pair<K,E> v) {
  return make_pair(v.first, curr.second + v.second);
}

//***********************************************************************************************
//***********************************************************************************************


//symmetric compact bipartite
struct bipartiteCSR {
  uintT *offsetsV, *offsetsU;
  uintE *edgesV, *edgesU;
  long nv, nu, numEdges;

  bipartiteCSR(uintT* _offsetsV, uintT* _offsetsU, uintE* _edgesV, uintE* _edgesU, long _nv, long _nu, long _ne) :
    offsetsV(_offsetsV), offsetsU(_offsetsU), edgesV(_edgesV), edgesU(_edgesU), nv(_nv), nu(_nu), numEdges(_ne)
  {}

  void del() {
    free(offsetsV); free(offsetsU); free(edgesV); free(edgesU);
  }
};

bipartiteCSR readBipartite(char* fname) {
  words W;
  _seq<char> S = readStringFromFile(fname);
  W = stringToWords(S.A, S.n);

  if (W.Strings[0] != (string) "AdjacencyHypergraph") {
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m -1;
  long nv = atol(W.Strings[1]);
  long mv = atol(W.Strings[2]);
  long nu = atol(W.Strings[3]);
  long mu = atol(W.Strings[4]);

  if ((len != nv + mv + nu + mu + 4) | (mv != mu)) {
    cout << "Bad input file" << endl;
    abort();
  }

  uintT* offsetsV = newA(uintT,nv+1);
  uintT* offsetsU = newA(uintT,nu+1);
  uintE* edgesV = newA(uintE,mv);
  uintE* edgesU = newA(uintE,mu);

  {parallel_for(long i=0; i < nv; i++) offsetsV[i] = atol(W.Strings[i + 5]);}
  offsetsV[nv] = mv;
  
  {parallel_for(long i=0; i<mv; i++) {
      edgesV[i] = atol(W.Strings[i+nv+5]);
      if(edgesV[i] < 0 || edgesV[i] >= nu) { cout << "edgesV out of range: nu = " << nu << " edge = " << edgesV[i] << endl; exit(0); }
    }}

  {parallel_for(long i=0; i < nu; i++) offsetsU[i] = atol(W.Strings[i + nv + mv + 5]);}
  offsetsU[nu] = mu;
  
  {parallel_for(long i=0; i<mu; i++) {
      edgesU[i] = atol(W.Strings[i+nv+mv+nu+5]);
      if(edgesU[i] < 0 || edgesU[i] >= nv) { cout << "edgesU out of range: nv = " << nv << " edge = " << edgesU[i] << endl; exit(0); }
    }}

  {parallel_for(long i=0; i < nu; ++i) {
  	// Sort from edgesU[offsetsU[i]] to edgesU[offsetsU[i+1]-1]
  	//quickSort(&edgesU[offsetsU[i]], offsetsU[i+1] - offsetsU[i], less<uintE>());
  }}

  {parallel_for(long i=0; i < nv; ++i) {
  	// Sort from edgesU[offsetsU[i]] to edgesU[offsetsU[i+1]-1]
  	//quickSort(&edgesV[offsetsV[i]], offsetsV[i+1] - offsetsV[i], less<uintE>());
  }}
  //W.del(); // to deal with performance bug in malloc
  return bipartiteCSR(offsetsV,offsetsU,edgesV,edgesU,nv,nu,mv);  
}

// Takes the elements of a vertex array, and returns the out degree choose 2
template <class E>
struct wedgeF { 
  uintT* offsets;
  wedgeF(uintT* _offsets) : offsets(_offsets) {}
  inline E operator() (const uintT& i) const {
    uintE v_deg = offsets[i+1]-offsets[i];
    return (E) ((v_deg * (v_deg-1)) / 2); 
  }
};

tuple<bool,long,long*> cmpWedgeCounts(bipartiteCSR & GA) {
  const long nv = GA.nv, nu = GA.nu;
  long num_wedges_v = sequence::reduce<long>((long) 0, nv, addF<long>(), wedgeF<long>(GA.offsetsV));
  long num_wedges_u = sequence::reduce<long>((long) 0, nu, addF<long>(), wedgeF<long>(GA.offsetsU));
  long* workPrefixSum;
  if(num_wedges_v <= num_wedges_u) {
    workPrefixSum = newA(long,nu);
    parallel_for(intT i=0;i<nu;i++) {
      long u_work = 0;
      intT u_offset = GA.offsetsU[i];
      intT u_deg = GA.offsetsU[i+1]-GA.offsetsU[i];
      for(intT j=0; j<u_deg;j++) {
  intT v = GA.edgesU[u_offset+j];
  intT v_offset = GA.offsetsV[v];
  intT v_deg = GA.offsetsV[v+1]-GA.offsetsV[v];
  // for(intT k=0; k<v_deg;k++) {
  //   uintE u2_idx = GA.edgesV[v_offset+k];
  //   if(u2_idx < i) u_work++;
  //   else break;
  
  // }
  u_work += GA.offsetsV[v+1]-GA.offsetsV[v];
      }
      workPrefixSum[i] = u_work;
    }
    sequence::plusScan(workPrefixSum,workPrefixSum,nu);
  } else {
    workPrefixSum = newA(long,nv);
    parallel_for(intT i=0;i<nv;i++) {
      long v_work = 0;
      intT v_offset = GA.offsetsV[i];
      intT v_deg = GA.offsetsV[i+1]-GA.offsetsV[i];
      for(intT j=0; j<v_deg;j++) {
  intT u = GA.edgesV[v_offset+j];
  intT u_offset = GA.offsetsU[u];
  intT u_deg = GA.offsetsU[u+1]-GA.offsetsU[u];
  // for(intT k=0; k<u_deg; k++) {
  //   uintE v2_idx = GA.edgesU[u_offset+k];
  //   if(v2_idx < i) v_work++;
  //   else break;
  // }
  v_work += GA.offsetsU[u+1]-GA.offsetsU[u];
      }
      workPrefixSum[i] = v_work;
    }
    sequence::plusScan(workPrefixSum,workPrefixSum,nv);
    //sequence::scan<long>(tuplePrefixSum,(long) 0, nu, addF<long>(),wedgeF<long>(GA.offsetsU), 0, false, false);
  }
  return make_tuple((num_wedges_v <= num_wedges_u), num_wedges_v <= num_wedges_u ? num_wedges_v : num_wedges_u, workPrefixSum);
}

pair<bool,long> cmpWedgeCountsSeq(bipartiteCSR & GA) {
  const long nv = GA.nv, nu = GA.nu;

  long num_wedges_v = 0;
  for(long i=0; i < nv; ++i) {
    uintE deg_v = GA.offsetsV[i+1]-GA.offsetsV[i];
    num_wedges_v += deg_v * (deg_v - 1) / 2;
  }

  long num_wedges_u = 0;
  for(long i=0; i < nu; ++i) {
    uintE deg_u = GA.offsetsU[i+1]-GA.offsetsU[i];
    num_wedges_u += deg_u * (deg_u - 1) / 2;
  }
  return make_pair((num_wedges_v <= num_wedges_u), num_wedges_v <= num_wedges_u ? num_wedges_v : num_wedges_u);
}
//***********************************************************************************************
//***********************************************************************************************

/*
 *  Sort objects using cmp, and then retrieve indices of where consecutive objects are
 *  not equal (as given by eq) to give frequency counts of each object in objs.
 *  In other words, if we let arr be the returned array, arr[i+1] - arr[i] is the 
 *  frequency of objs[arr[i]] (objs[arr[i]] = ... = objs[arr[i+1]-1]).
 *  Iterating from 0 (inclusive) to the returned length - 1 (exclusive) will give
 *  all frequency counts.
 * 
 *  objs: Objects to count frequencies of
 *  num : Length of objs
 *  cmp : Comparator for T objects
 *  eq  : Equality comparator for T objects
 *  sort: If objs is already sorted, then this should be set to false so we don't resort;
 *        by default, we assume objs is not sorted.
 * 
 *  Returns: Array and length of array with frequency counts (as described above)
 */
template <class T, class Cmp, class Eq>
pair<uintE*, long> getFreqs(T* objs, long num, Cmp cmp, Eq eq, bool sort=true) {
  // Sort objects
  if (sort) sampleSort(objs, num, cmp);

  uintE* freqs = newA(uintE, num + 1);
  freqs[0] = 0;
  freqs[num] = num;

  // Retrieve indices where objects differ
  parallel_for(long i=1; i < num; ++i) {
    if (!eq(objs[i-1],objs[i])) freqs[i] = i;
    else freqs[i] = UINT_E_MAX;
  }
  uintE* freqs_f = newA(uintE, num+1);
  long num_freqs_f = sequence::filter(freqs, freqs_f, num+1, nonMaxF());
  free(freqs);
  return make_pair(freqs_f, num_freqs_f);
}

template <class S, class T, class Cmp, class Eq, class OpT, class OpuintE, class OpCount>
pair<tuple<S,uintE>*, long> getFreqs_seq(T* objs, long num, Cmp cmp, Eq eq, bool sort=true, OpT opt=refl<T>(),
  OpuintE opuinte=refl<uintE>(), OpCount opcount=reflCount<T>()) {
  if(sort) sampleSort(objs, num, cmp);

  using X = tuple<S,uintE>;
  X* freqs = newA(X, num);
  T prev = objs[0];
  T curr = objs[0];
  long idx = 0;
  long count = opcount(prev);
  for(long i=1; i < num; ++i) {
    curr = objs[i];
    if (!eq(prev, curr)) {
      freqs[idx] = make_tuple(opt(prev), opuinte(count));
      idx++;
      count = opcount(curr);
      prev = curr;
    }
    else {
      count += opcount(curr);
    }
  }
  freqs[idx] = make_tuple(opt(curr), opuinte(count));
  return make_pair(freqs, idx + 1);
}

//********************************************************************************************
//********************************************************************************************

template <class T>
tuple<long, T*> intersect_seq(T* a, T* b, size_t num_a, size_t num_b) {
  if (num_b < num_a) return intersect_seq(b, a, num_b, num_a);
  //sampleSort(a,num_a,cmpF<T>());
  uintE* ret = newA(uintE,num_a);
  long idx = 0;

  auto a_map = make_in_imap<uintT>(num_a, [&] (size_t i) { return a[i]; });
  auto lt = [] (const T& l, const T& r) { return l < r; };

  for(size_t i=0;i<num_b;++i) {
    size_t find_idx = pbbs::linear_search(a_map, b[i], lt);
    if(a[find_idx] == b[i]) ret[idx++] = b[i];
  }
  return make_tuple(idx, ret);
}

template <class T>
tuple<long, T*> intersect(T* a, T* b, size_t num_a, size_t num_b) {
  if (num_a < 1000 && num_b < 1000) return intersect_seq(a,b,num_a,num_b);
  if (num_b < num_a) return intersect(b, a, num_b, num_a);

  uintE* a_cpy = newA(uintE, num_a);
  parallel_for(size_t i=0; i < num_a; ++i) {a_cpy[i] = a[i];};

  sampleSort(a_cpy,num_a,cmpF<T>());
  uintE* ret = newA(uintE,num_a);
  long idx = 0;

  auto a_map = make_in_imap<uintT>(num_a, [&] (size_t i) { return a_cpy[i]; });
  auto lt = [] (const T& l, const T& r) { return l < r; };

  for(size_t i=0;i<num_b;++i) {
    size_t find_idx = pbbs::binary_search(a_map, b[i], lt);
    if(a[find_idx] == b[i]) ret[idx++] = b[i];
  }
  return make_tuple(idx, ret);
}

tuple<long, uintE*> intersect_hist(uintE* a, uintE* b, size_t num_a, size_t num_b, uintE max) {
  if (num_a < 1000 && num_b < 1000) return intersect_seq(a,b,num_a,num_b);
  using X = tuple<uintE,uintE>;
  size_t num_total = num_a + num_b;
  uintE* total = newA(uintE, num_total);

  parallel_for(long i=0; i < num_a; ++i) {total[i] = a[i]; }
  parallel_for(long i=0; i < num_b; ++i) {total[i+num_a] = b[i]; }

  pbbsa::sequence<uintE> seq_total = pbbsa::sequence<uintE>(total, num_total);
  tuple<size_t,X*> hist_total = pbbsa::sparse_histogram<uintE,uintE>(seq_total, max);

  X* arr = newA(X, get<0>(hist_total));
  long len = sequence::filter(get<1>(hist_total), arr, get<0>(hist_total), greaterOneF());

  uintE* ret = newA(uintE,len);
  parallel_for(long i=0; i < len; ++i) { ret[i] = get<0>(arr[i]); }

  free(total);
  free(get<1>(hist_total));
  free(arr);

  return make_tuple(len, ret);
}

tuple<long, uintE*> intersect_hash_seq(uintE* a, uintE* b, size_t num_a, size_t num_b) {
  if (num_b < num_a) return intersect_hash_seq(b, a, num_b, num_a);

  sparseAdditiveSet<uintE> set = sparseAdditiveSet<uintE>(num_a, 1, UINT_E_MAX);
  for(size_t i=0; i<num_a; ++i) { set.insert(make_pair(a[i],0)); }

  uintE* ret = newA(uintE,num_a);
  long idx = 0;
  for(size_t i=0; i<num_b; ++i) {
    if(set.find(b[i]).first != UINT_E_MAX) ret[idx++] = b[i];
  }

  set.del();
  return make_tuple(idx, ret);
}

tuple<long, uintE*> intersect_hash(uintE* a, uintE* b, size_t num_a, size_t num_b) {
  //if (num_b < num_a) return intersect_hash(b, a, num_b, num_a);
  if (num_a < 1000 && num_b < 1000) return intersect_hash_seq(a,b,num_a,num_b);
  using X=pair<uintE,uintE>;

  sparseAdditiveSet<uintE> set = sparseAdditiveSet<uintE>(num_a, 1, UINT_E_MAX);
  parallel_for(size_t i=0; i<num_a; ++i) { set.insert(make_pair(a[i],0)); }
  parallel_for(size_t i=0; i<num_b; ++i) { set.insert(make_pair(b[i],1)); }

  _seq<pair<uintE,uintE>> seq = set.entries();
  X* seq_f = newA(X, seq.n);
  long len = sequence::filter(seq.A, seq_f, seq.n, nonZeroPairF());

  uintE* ret = newA(uintE,len);
  parallel_for(long i=0; i < len; ++i) { ret[i] = seq_f[i].first; }

  free(seq_f);
  seq.del();
  set.del();
  return make_tuple(len, ret);
}

//********************************************************************************************
//********************************************************************************************

//TODO we don't use nbhrs or save nbhrs -- delete from everything
uintE* countWedgesScan(bipartiteCSR& GA, bool use_v, bool half=false) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  uintE* idxs = newA(uintE, nu + 1);
  idxs[nu] = 0;

  using T = uintT*;
  T* nbhd_idxs = newA(T, nu);

  parallel_for(intT i=0; i < nu; ++i) {
    idxs[i] = 0;
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;

    nbhd_idxs[i] = newA(uintE, u_deg + 1);
    parallel_for(intT j=0; j < u_deg+1; ++j) {(nbhd_idxs[i])[j] = 0;}

    parallel_for(intT j=0; j < u_deg; ++j) { //TODO can parallelize this too technically
      if (!half) {
        uintE v = edgesU[u_offset + j];
        (nbhd_idxs[i])[j] = offsetsV[v+1] - offsetsV[v] - 1;//V[U[i].getOutNeighbor(j)].getOutDegree() - 1;
      }
      else {
        (nbhd_idxs[i])[j] = 0;
        uintE v = edgesU[u_offset + j];
        intT v_offset = offsetsV[v];
        intT v_deg = offsetsV[v+1] - v_offset;
        for (intT k = 0; k < v_deg; ++k) { //TODO can parallelize this too technically
          if (edgesV[v_offset + k] < i) nbhd_idxs[i][j] ++;
          else break;
        }
      }
    }

    idxs[i] = sequence::plusReduce(nbhd_idxs[i], u_deg + 1);
    free(nbhd_idxs[i]);
  }
  free(nbhd_idxs);
  sequence::plusScan(idxs, idxs, nu + 1);

  return idxs;
}


//***************************************************************************************************
//***************************************************************************************************

template<class wedge, class wedgeCons, class Sequence>
void _getActiveWedges_seq(_seq<wedge>& wedges_seq, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges,
  intT curr_idx, intT next_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long idx = 0;
  for(intT i=curr_idx; i < next_idx; ++i) {
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset + k];
        if (u2 != u) {
          wedges_seq.A[idx] = cons(u, u2, v);
          ++idx;
        }
      }
    }
  }
}

template<class wedge, class wedgeCons, class Sequence>
void _getActiveWedges(_seq<wedge>& wedges_seq, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges,
  intT curr_idx=0, intT next_idx=INT_T_MAX) {

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (wedges_seq.n < num_wedges) {
    free(wedges_seq.A);
    wedges_seq.A = newA(wedge, num_wedges);
    wedges_seq.n = num_wedges;
  }

  if (next_idx == INT_T_MAX) next_idx = num_I;
  if (num_wedges < 10000) return _getActiveWedges_seq<wedge>(wedges_seq, I, num_I, GA, use_v, cons, num_wedges, curr_idx, next_idx);

  // First, we must retrive the indices at which to store each seagull (so that storage can be parallelized)
  // Array of indices associated with seagulls for each active vertex
  long* sg_idxs = newA(long, next_idx - curr_idx + 1);
  sg_idxs[next_idx-curr_idx] = 0;

  // 2D array of indices associated with seagulls for each neighbor of each active vertex
  using T = long*;
  T* nbhd_idxs = newA(T, next_idx-curr_idx); 

  parallel_for(intT i=curr_idx; i < next_idx; ++i) {
    // Set up for each active vertex
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    long u_deg = offsetsU[u+1] - u_offset;
    // Allocate space for indices associated with each neighbor of u
    nbhd_idxs[i-curr_idx] = newA(long, u_deg + 1);
    (nbhd_idxs[i-curr_idx])[u_deg] = 0;
    // Assign indices associated with each neighbor of u
    parallel_for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      (nbhd_idxs[i-curr_idx])[j] = offsetsV[v+1] - offsetsV[v] - 1;
    }
    sequence::plusScan(nbhd_idxs[i-curr_idx], nbhd_idxs[i-curr_idx], u_deg+1);
    // Set up indices associated with u
    sg_idxs[i-curr_idx] = (nbhd_idxs[i-curr_idx])[u_deg];
  }
  // Assign indices associated with each active vertex
  sequence::plusScan(sg_idxs, sg_idxs, next_idx-curr_idx + 1);

  // Store seagulls in parallel
  parallel_for(intT i=curr_idx; i < next_idx; ++i) {
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    long sg_idx = sg_idxs[i-curr_idx];
    // Consider each neighbor v of active vertex u
    parallel_for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      long nbhd_idx = (nbhd_idxs[i-curr_idx])[j];
      long idx = 0;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset + k];
        if (u2 != u) {
          wedges_seq.A[sg_idx+nbhd_idx+idx] = cons(u, u2, v);
          ++idx;
        }
      }
    }
  }

  // Cleanup
  parallel_for(intT i=curr_idx; i<next_idx;++i) { 
    free(nbhd_idxs[i-curr_idx]); 
  }
  free(nbhd_idxs);
  free(sg_idxs);
}

template<class wedgeCons, class T, class Sequence>
void _getActiveWedgesHash(T& wedges, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, 
  intT curr_idx=0, intT next_idx=INT_T_MAX) {

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (next_idx == INT_T_MAX) next_idx = num_I;
  wedges.resize(num_wedges);

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    parallel_for (intT j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find all seagulls with center v and endpoint u
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2 = edgesV[v_offset + k];
        if (u2 != u) wedges.insert(make_pair(cons(u,u2,v),1));
	    }
    }
  }
}

template<class Sequence>
pair<long,intT> getNextActiveWedgeIdx_seq(Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long num_wedges = 0;
  for(intT i=curr_idx; i < num_I; ++i) {
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    long num = 0;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_deg = offsetsV[v+1] - offsetsV[v];
      num += (v_deg - 1);
      if (num_wedges + num > max_wedges) {
        if (i == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
        return make_pair(num_wedges, i);
      }
    }
    num_wedges += num;
  }
  return make_pair(num_wedges, num_I);
}

// TODO based on half can also optimize this stuff? but will be hard to parallelize later so maybe not
//JS: there looks to be some inefficiency here. we should have gotten the prefix sum of wedges for all vertices at the beginning, so we don't need to recompute it here. binary search can use a doubling search instead, starting at the last index until getting the right range, then binary search inside.
template<class Sequence>
pair<long, intT> getNextActiveWedgeIdx(Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (num_I - curr_idx < 10000) return getNextActiveWedgeIdx_seq(I, num_I, GA, use_v, max_wedges, curr_idx);
  long* idxs = newA(long, num_I - curr_idx + 1);
  idxs[num_I-curr_idx] = 0;
  parallel_for(intT i=curr_idx; i < num_I; ++i) {
    idxs[i-curr_idx] = 0;
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_deg = offsetsV[v+1] - offsetsV[v];
      idxs[i-curr_idx] += (v_deg - 1);
    }
  }
  sequence::plusScan(idxs, idxs, num_I - curr_idx + 1);

  auto idx_map = make_in_imap<long>(num_I - curr_idx, [&] (size_t i) { return idxs[i+1]; });
  auto lte = [] (const long& l, const long& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx; //this rets first # > searched num
  long num = idxs[find_idx - curr_idx];
  free(idxs);
  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }

  return make_pair(num, find_idx); //TODO make sure right
}

template<class wedge, class wedgeCons, class Sequence>
pair<long, intT> getActiveWedges(_seq<wedge>& wedges_seq, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons,
  long max_wedges, intT curr_idx, long num_wedges) {
  if (max_wedges >= num_wedges) {
    _getActiveWedges<wedge>(wedges_seq, I, num_I, GA, use_v, cons, num_wedges);
    return make_pair(num_wedges, num_I);
  }

  pair<long, intT> p = getNextActiveWedgeIdx(I, num_I, GA, use_v, max_wedges, curr_idx);
  _getActiveWedges<wedge>(wedges_seq, I, num_I, GA, use_v, cons, p.first, curr_idx, p.second);
  return p;
}

template<class wedgeCons, class T, class Sequence>
long getActiveWedgesHash(T& wedges, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons, long max_wedges, intT curr_idx, long num_wedges) {
  if (max_wedges >= num_wedges) {
    _getActiveWedgesHash(wedges, I, num_I, GA, use_v, cons, num_wedges);
    return num_I;
  }
  pair<long, intT> p = getNextActiveWedgeIdx(I, num_I, GA, use_v, max_wedges, curr_idx);
  _getActiveWedgesHash(wedges, I, num_I, GA, use_v, cons, p.first, curr_idx, p.second);
  return p.second;
}

//***************************************************************************************************
//***************************************************************************************************

template<class wedge, class wedgeCons>
void _getWedges_seq(wedge* wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, intT curr_idx, intT next_idx) {
  const long nu = use_v ? GA.nu : GA.nv;
  const long nv = use_v ? GA.nv : GA.nu;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long idx = 0;
  for(intT i=curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset+k];
        if (u2 < i) {
          wedges[idx] = cons(i, u2, v);
          ++idx;
        }
        else break; 
      }
    }
  }
}

template<class wedge, class wedgeCons>
void _getWedges(_seq<wedge>& wedges_seq, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, uintE* wedge_idxs, intT curr_idx=0, intT next_idx=INT_T_MAX) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Allocate space for seagull storage
  if (wedges_seq.n < num_wedges) {
    free(wedges_seq.A);
    wedges_seq.A = newA(wedge, num_wedges);
    wedges_seq.n = num_wedges;
  }

  if (next_idx == INT_T_MAX) next_idx = nu;
  if (num_wedges < 10000) return _getWedges_seq<wedge>(wedges_seq.A, GA, use_v, cons, num_wedges, curr_idx, next_idx);
 
  // Store seagulls in parallel
  parallel_for(intT i=curr_idx; i < next_idx; ++i) {
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    // Consider each neighbor v of active vertex u
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    long idx = 0;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        intT u2 = edgesV[v_offset+k];
        if (u2 < i) {
          wedges_seq.A[wedge_idx+idx] = cons(i, u2, v);
          ++idx;
        }
        else break;
      }
    }
  }
}

template<class wedgeCons, class T>
void _getWedgesHash(T& wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, intT curr_idx=0, intT next_idx=INT_T_MAX) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (next_idx == INT_T_MAX) next_idx = nu;
  wedges.resize(num_wedges);
  //hashInsertTimer.start();
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    parallel_for (long j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find all seagulls with center v and endpoint u
      for (long k=0; k < v_deg; ++k) { 
        uintE u2 = edgesV[v_offset+k];
        if (u2 < i) wedges.insert(make_pair(cons(i,u2,v),1));
        else break;
      }
    }
  }
  //hashInsertTimer.stop();
}

intT getNextWedgeIdx_seq(bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long orig = max_wedges;
  for(intT i=curr_idx; i < nu; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      uintE num = 0;
      for (intT k=0; k < v_deg; ++k) {
        if (edgesV[v_offset+k] < i) num ++;
        else break;
      }
      if (num > max_wedges) {
        if (i == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
        return i;
      }
      else { max_wedges -= num; }
    }
  }
  return nu;
}

// TODO doubling search
intT getNextWedgeIdx(bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx, uintE* wedge_idxs) {
  const long nu = use_v ? GA.nu : GA.nv;
  nextWedgeTimer.start();
  if (nu - curr_idx < 2000) return getNextWedgeIdx_seq(GA, use_v, max_wedges, curr_idx);

  auto idx_map = make_in_imap<uintE>(nu - curr_idx, [&] (size_t i) { return wedge_idxs[curr_idx+i+1] - wedge_idxs[curr_idx]; });
  auto lte = [] (const uintE& l, const uintE& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx; //this rets first # > searched num
  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
  nextWedgeTimer.stop();
  return find_idx; //TODO make sure right
}

//TODO 3 tuple instead of nested pairs
template<class wedge, class wedgeCons>
pair<long, intT> getWedges(_seq<wedge>& wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long max_wedges, intT curr_idx, long num_wedges, uintE* wedge_idxs) {
  const long nu = use_v ? GA.nu : GA.nv;
  if (max_wedges >= num_wedges) {
    _getWedges<wedge>(wedges, GA, use_v, cons, num_wedges, wedge_idxs);
    return make_pair(num_wedges, nu);
  }
  long next_idx = getNextWedgeIdx(GA, use_v, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  _getWedges<wedge>(wedges, GA, use_v, cons, num_wedges, wedge_idxs, curr_idx, next_idx);
  return make_pair(num_wedges, next_idx);
}

template<class wedgeCons, class T>
intT getWedgesHash(T& wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long max_wedges, intT curr_idx, long num_wedges, uintE* wedge_idxs) {
  const long nu = use_v ? GA.nu : GA.nv;
  if (max_wedges >= num_wedges) {
    _getWedgesHash(wedges, GA, use_v, cons, num_wedges);
    return nu;
  }
  intT next_idx = getNextWedgeIdx(GA, use_v, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  _getWedgesHash(wedges, GA, use_v, cons, num_wedges, curr_idx, next_idx);
  return next_idx;
}

//***************************************************************************************************
//***************************************************************************************************

pair<tuple<uintE,uintE>*,long> updateBuckets(uintE* update_idxs, long num_updates, uintE* butterflies, 
  array_imap<uintE> D, buckets<array_imap<uintE>> b, uintE k) {
  using X = tuple<uintE,uintE>;
  X* update = newA(X,num_updates);

    // Filter for bucket updates
  parallel_for(long i=0; i < num_updates; ++i) {
    const uintE u_idx = update_idxs[i];
    uintE old_b = D.s[u_idx];

    if (old_b > k) {
        uintE new_b = max(butterflies[u_idx],k);
        D.s[u_idx] = new_b;
        uintE new_bkt = b.get_bucket(old_b, new_b);
        update[i] = make_tuple(u_idx, new_bkt);
    }
    else {update[i] = make_tuple(UINT_E_MAX,UINT_E_MAX);}
  }

  X* update_filter = newA(X, num_updates);
  long num_updates_filter = sequence::filter(update,update_filter,num_updates, nonMaxTupleF());
  free(update);
  return make_pair(update_filter,  num_updates_filter);
}

pair<tuple<uintE,uintE>*,long> updateBuckets_seq(uintE* update_idxs, long num_updates, uintE* butterflies, 
  array_imap<uintE> D, buckets<array_imap<uintE>> b, uintE k) {
  using X = tuple<uintE, uintE>;
  X* update = newA(X, num_updates);
  long idx = 0;
  for(long i=0; i < num_updates; ++i) {
    uintE u_idx = update_idxs[i];
    uintE old_b = D.s[u_idx];

    if(old_b > k) {
      uintE new_b = max(butterflies[u_idx], k);
      D.s[u_idx] = new_b;
      uintE new_bkt = b.get_bucket(old_b, new_b);
      update[idx] = make_tuple(u_idx, new_bkt);
      ++idx;
    }
  }
  return make_pair(update, idx);
}

struct edgeToIdx { 
  long nv;
  long nu;
  uintT* offsetsV;
  uintT* offsetsU;
  uintE* edgesV;
  uintE* edgesU;

  long numEdges;
  long max_wedges;
  sparseAdditiveSet<uintE> edges;

  edgeToIdx(bipartiteCSR& GA, bool use_v, long _max_wedges) : max_wedges(_max_wedges) {
    nu = use_v ? GA.nu : GA.nv;
    nv = use_v ? GA.nv : GA.nu;
    offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
    offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
    edgesV = use_v ? GA.edgesV : GA.edgesU;
    edgesU = use_v ? GA.edgesU : GA.edgesV;
    numEdges = GA.numEdges;

    edges = sparseAdditiveSet<uintE>(numEdges, 1, UINT_E_MAX);
    
    if (nu*nv < max_wedges) {
      bool* edges_bool = newA(bool, nu*nv);
      parallel_for(long i=0; i < nu*nv; ++i) {edges_bool[i] = false;}
      parallel_for(intT i=0; i < nv; ++i) {
      	intT v_offset = offsetsV[i];
      	intT v_deg = offsetsV[i+1] - v_offset;
        for (intT j = 0; j < v_deg; ++j) {
          intT u = edgesV[v_offset + j];
          edges_bool[nu*i + u] = true;
        }
      }
      auto f = [&] (size_t i) { return edges_bool[i]; };
      auto f_in = make_in_imap<bool>(nu*nv, f);
      auto out = pbbs::pack_index<uintE>(f_in);
      out.alloc = false;
      free(edges_bool);

      parallel_for(long i=0; i < out.size(); ++i) {
        edges.insert(make_pair(out.s[i], i));
      }

      free(out.s);
    }
    else {
      parallel_for (intT i=0; i < nv; ++i) {
        intT v_offset = offsetsV[i];
        intT v_deg = offsetsV[i+1] - v_offset;
        for (long j=0; j < v_deg; ++j) {
          intT u = edgesV[v_offset + j];
          edges.insert(make_pair(nu*i + u, 0));
        }
      }
      auto edges_seq = edges.entries();
      sequence::scan<pair<uintE,uintE>>(edges_seq.A, (long) 0, edges_seq.n, getAddReducePair<uintE,uintE>,
        sequence::getA<pair<uintE,uintE>, long>(edges_seq.A), make_pair(0,0), false, false);
      parallel_for(long i=0; i < edges_seq.n; ++i) { edges.insert(edges_seq.A[i]);}
      edges_seq.del();
    }
  }

  void del() {
    edges.del();
  }

  // Note: always in the format (V, U)
  inline uintE operator() (const uintE& i, const uintE& j) {
    return (edges.find((uintE) (nu*i + j))).second;
  }
};


#endif

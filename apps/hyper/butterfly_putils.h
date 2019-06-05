#ifndef _BPUTILS_
#define _BPUTILS_

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

#include "butterfly_utils.h"

#define MAX_STEP_SIZE 1000

using namespace std;

struct PeelESpace {
  long type;
  long nu;
  long stepSize;
  long n_side;
  _seq<tuple<uintE,uintE>> wedges_seq_tup;
  _seq<tuple<uintE,uintE>> wedges_seq_tup_fil;
  _seq<uintE> update_seq_int;
  sparseAdditiveSet<uintE> update_hash;
  _seq<uintE> wedges_seq_int;
  _seq<uintE> used_seq_int;
PeelESpace(long _type, long _nu, long _stepSize, long _n_side) : type(_type), nu(_nu), stepSize(_stepSize), n_side(_n_side) {
  using X = tuple<uintE,uintE>;
  update_seq_int = _seq<uintE>(newA(uintE, nu), nu);
  if (type == 0) update_hash = sparseAdditiveSet<uintE>(nu, (float) 1, UINT_E_MAX);
  else if (type == 1 || type == 2) {
    wedges_seq_tup = _seq<X>(newA(X, nu), nu);
    wedges_seq_tup_fil = _seq<X>(newA(X, nu), nu);
  }
  else {
    wedges_seq_int = _seq<uintE>(newA(uintE, n_side*stepSize), n_side*stepSize);
    granular_for(i,0,n_side*stepSize,n_side*stepSize > 10000, { wedges_seq_int.A[i] = 0; });
    used_seq_int = _seq<uintE>(newA(uintE, n_side*stepSize), n_side*stepSize);
  }
}
  void resize_update(size_t size) {
  	if (update_seq_int.n < size) {
      free(update_seq_int.A);
      update_seq_int.A = newA(uintE, size);
      update_seq_int.n = size;
    }
  }

  void clear() {
    if (type == 0) update_hash.clear();
  }

  void del() {
  	update_seq_int.del();
    if (type == 0) update_hash.del();
    else if (type == 1 || type == 2) { wedges_seq_tup.del(); wedges_seq_tup_fil.del(); }
    else {wedges_seq_int.del(); used_seq_int.del();}
  }
};

struct PeelSpace {
  long type;
  long nu;
  long stepSize;
  _seq<UVertexPair> wedges_seq_uvp;
  _seq<uintE> wedges_seq_int;
  _seq<long> wedges_seq_long;
  _seq<uintE> used_seq_int;
  _seq<uintE> update_seq_int;
  sparseAdditiveSet<uintE> update_hash;
  sparseAdditiveSet<long, long>* wedges_hash;
  sparseAdditiveSet<long, long>** wedges_hash_list;
  _seq<pair<long,long>> wedges_seq_intp;
  _seq<pair<uintE,uintE>> butterflies_seq_intp;
  intT num_wedges_hash;
PeelSpace(long _type, long _nu, long _stepSize) : type(_type), nu(_nu), stepSize(_stepSize) {
  using E = pair<long, long>;
  using X = pair<uintE,uintE>;
  update_seq_int = _seq<uintE>(newA(uintE, nu), nu);
  if (type == 0) {
    using T = sparseAdditiveSet<long, long>*;
    wedges_hash = new sparseAdditiveSet<long, long>(1,1, LONG_MAX, LONG_MAX);
    wedges_hash_list = newA(T, 1);
    wedges_hash_list[0] = wedges_hash;
    num_wedges_hash = 1;
    update_hash = sparseAdditiveSet<uintE>(nu, (float) 1, UINT_E_MAX);
    wedges_seq_intp = _seq<E>(newA(E, nu), nu);
    butterflies_seq_intp = _seq<X>(newA(X, nu), nu);
  }
  else if (type == 1) wedges_seq_uvp = _seq<UVertexPair>(newA(UVertexPair, nu), nu);
  else if (type == 2) wedges_seq_long = _seq<long>(newA(long, nu), nu);
  else {
    timer t1;
    wedges_seq_int = _seq<uintE>(newA(uintE, nu*stepSize), nu*stepSize);
    t1.start();
    granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges_seq_int.A[i] = 0; });
    t1.reportTotal("time for init wedges");
    used_seq_int = _seq<uintE>(newA(uintE, nu*stepSize), nu*stepSize);
  }
}
  void resize_update(size_t size) {
  	if (update_seq_int.n < size) {
      free(update_seq_int.A);
      update_seq_int.A = newA(uintE, size);
      update_seq_int.n = size;
    }
  }
  
  sparseAdditiveSet<long, long>* resize(size_t size) {
    if (type != 0) return nullptr;
    
    size_t find_idx = log2RoundUp(size);
    if (find_idx < num_wedges_hash) {
      wedges_hash = wedges_hash_list[find_idx];
      wedges_hash->clear(); 
      return wedges_hash_list[find_idx];
    }
    using T = sparseAdditiveSet<long, long>*;
    T* new_wedges_hash_list = newA(T, find_idx+1);
    parallel_for(intT i=0; i < num_wedges_hash; ++i) {
      new_wedges_hash_list[i] = wedges_hash_list[i];
    }
    parallel_for(intT i=num_wedges_hash; i < find_idx+1; ++i) {
      new_wedges_hash_list[i] = new sparseAdditiveSet<long, long>(1u << i,1, LONG_MAX, LONG_MAX);
    }
    free(wedges_hash_list);
    wedges_hash_list = new_wedges_hash_list;
    num_wedges_hash = find_idx+1;
    wedges_hash = wedges_hash_list[find_idx];
    return wedges_hash_list[find_idx];
    // binary search wedges_hash_list for right size
    // if none found, square top until we get to size
  }

  void clear() {
    if (type == 0) update_hash.clear();
  }

  void del() {
  	update_seq_int.del();
    if (type == 0) {
      parallel_for(intT i=0; i < num_wedges_hash; ++i) { free(wedges_hash_list[i]); }
      free(wedges_hash_list);
      update_hash.del();
      wedges_seq_intp.del(); 
      butterflies_seq_intp.del();
    }
    else if (type == 1) wedges_seq_uvp.del();
    else if (type == 2) wedges_seq_long.del();
    else {wedges_seq_int.del(); used_seq_int.del();}
  }
};


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
          wedges_seq.A[idx] = cons(u, u2, v, j, k);
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
          wedges_seq.A[sg_idx+nbhd_idx+idx] = cons(u, u2, v, j, k);
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

template<class wedgeCons, class Sequence>
void _getActiveWedgesHash(PeelSpace& ps, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, 
  intT curr_idx=0, intT next_idx=INT_T_MAX) {

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (next_idx == INT_T_MAX) next_idx = num_I;
  auto wedges = ps.resize(num_wedges);
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
        if (u2 != u) wedges->insert(make_pair(cons(u, u2, v, j, k),1));
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

template<class wedgeCons, class Sequence>
long getActiveWedgesHash(PeelSpace& ps, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons, long max_wedges, intT curr_idx, long num_wedges) {
  if (max_wedges >= num_wedges) {
    _getActiveWedgesHash(ps, I, num_I, GA, use_v, cons, num_wedges);
    return num_I;
  }
  pair<long, intT> p = getNextActiveWedgeIdx(I, num_I, GA, use_v, max_wedges, curr_idx);
  _getActiveWedgesHash(ps, I, num_I, GA, use_v, cons, p.first, curr_idx, p.second);
  return p.second;
}


//***************************************************************************************************
//***************************************************************************************************

// store as tuple uintE, uintE -- so hist can use
/*template<class Sequence>
void _getIntersectWedges_seq(_seq<tuple<uintE,uintE>>& wedges_seq, uintE* eti, uintE* ite, bool* current, Sequence I, intT num_I, 
  bipartiteCSR& GA, bool use_v, long num_wedges, intT curr_idx, intT next_idx) {
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
          wedges_seq.A[idx] = cons(u, u2, v, j, k);
          ++idx;
        }
      }
    }
  }
}*/

void intersect(_seq<tuple<uintE,uintE>>& wedges_seq, uintE* eti, uintE* ite, bool* current, 
  bipartiteCSR& GA, bool use_v, uintE u, uintE v, uintE u2, uintE idx_vu, uintE idx_vu2, intT wedges_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  // have u, v, u2
  intT u2_offset = offsetsU[u2];
  intT u2_deg = offsetsU[u2+1] - u2_offset;

  intT u_offset = offsetsU[u];
  intT u_deg = offsetsU[u+1] - u_offset;
  bool* same = newA(bool, u2_deg);
  //auto u_map = make_in_imap<uintT>(u_deg, [&] (size_t i) { return edgesU[u_offset + i]; });
  //auto lt = [] (const uintE& l, const uintE& r) { return l < r; };
  //parallel_
  intT int_offset = 0;
  parallel_for(intT j = 0; j < u2_deg; ++j) { same[j] = 0; }
  for(intT j = 0; j < u2_deg; ++j) {
    uintE v2 = edgesU[u2_offset + j];
    uintE idx_v2u2 = eti[u2_offset + j];
    if(v2 != UINT_E_MAX && v2 != v && (!current[idx_v2u2] || idx_v2u2 >= idx_vu)) {
      while(int_offset < u_deg && (edgesU[u_offset + int_offset] == UINT_E_MAX || edgesU[u_offset + int_offset] < v2)) { int_offset++; }
      if (int_offset < u_deg && edgesU[u_offset+int_offset] == v2 && (!current[eti[u_offset + int_offset]] || idx_vu < eti[u_offset + int_offset])) {
        same[j] = 1;//edgesU[u2_offset + j];
        //assert(wedges_idx+1 < wedges_seq.n);
        wedges_seq.A[wedges_idx++] = make_tuple(eti[u2_offset + j],1);
        //assert(wedges_idx+1 < wedges_seq.n);
        wedges_seq.A[wedges_idx++] = make_tuple(eti[u_offset + int_offset],1);
      }
      else if(int_offset >= u_deg) break;
    }
  }
  long num_same = sequence::sum(same, u2_deg);
  //assert(wedges_idx+1 < wedges_seq.n);
  if (num_same > 0) { wedges_seq.A[wedges_idx++] = make_tuple(idx_vu2, num_same); }
  free(same);
}

template<class Sequence>
long _getIntersectWedges(PeelESpace& ps, uintE* eti, uintE* ite, bool* current, Sequence I, intT num_I, 
  bipartiteCSR& GA, bool use_v, long num_wedges, long* idxs, intT curr_idx=0, intT next_idx=INT_T_MAX) {
  using X = tuple<uintE,uintE>;

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (ps.wedges_seq_tup.n < num_wedges*2+next_idx-curr_idx) {
    long num = num_wedges*2+next_idx-curr_idx;
    free(ps.wedges_seq_tup.A);
    free(ps.wedges_seq_tup_fil.A);
    ps.wedges_seq_tup.A = newA(X, num);
    ps.wedges_seq_tup_fil.A = newA(X, num);
    ps.wedges_seq_tup.n = num;
    ps.wedges_seq_tup_fil.n = num;
  }
  parallel_for(long i=0; i < num_wedges*2+next_idx-curr_idx; ++i) { ps.wedges_seq_tup.A[i] = make_tuple(UINT_E_MAX, UINT_E_MAX); }

  if (next_idx == INT_T_MAX) next_idx = num_I;
  //if (num_wedges < 10000)
  //  return _getIntersectWedges_seq(wedges_seq, eti, ite, current, I, num_I, GA, use_v, num_wedges, curr_idx, next_idx);
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    uintE idx = I[i];
    uintE u = edgesV[idx];
    uintE v = edgesU[ite[idx]];
    intT v_offset = offsetsV[v];
    intT v_deg = offsetsV[v+1] - v_offset;
  // we want to go through all neighbors of v
  // intersect each N(u') with N(u)
  // get all u, v, u2
  // intersection of N(u) and N(u2) - 1 is the num butterflies to subtract from v, u2
    parallel_for (intT k=0; k < v_deg; ++k) { 
      uintE u2 = edgesV[v_offset + k];
      if (u2 != u && u2 != UINT_E_MAX && (!current[v_offset+k] || idx < v_offset + k)) {//u2 != UINT_E_MAX - 1
        intersect(ps.wedges_seq_tup, eti, ite, current, GA, use_v, u, v, u2, idx, v_offset+k, idxs[i-curr_idx]*2 + i - curr_idx); //TODO make sure this is ok
      }
	  }
  }
  long num_wedges_fil = sequence::filter(ps.wedges_seq_tup.A, ps.wedges_seq_tup_fil.A, num_wedges*2+next_idx-curr_idx, nonMaxTupleF());
  return num_wedges_fil;
}

template<class Sequence>
tuple<long, intT, long*> getNextIntersectWedgeIdx_seq(uintE* eti, uintE* ite, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long num_wedges = 0;
  long curr_num_wedges = 0;
  long* idxs = newA(long, num_I - curr_idx + 1);
  idxs[0] = 0;
  for(intT i=curr_idx; i < num_I; ++i){
    // Set up for each active vertex
    uintE idx = I[i];
    uintE u = edgesV[idx];
    uintE v = edgesU[ite[idx]];
    intT v_offset = offsetsV[v];
    intT v_deg = offsetsV[v+1] - v_offset;
    intT u_deg = offsetsU[u+1] - offsetsU[u];
    curr_num_wedges = 0;
    for (intT k=0; k < v_deg; ++k) {
      uintE u2 = edgesV[v_offset + k];
      if (u2 != u && u2 != UINT_E_MAX && (u2 != UINT_E_MAX - 1 || idx < v_offset + k)) {
        intT u2_deg = offsetsU[u2+1] - offsetsU[u2];
        curr_num_wedges += min(u_deg, u2_deg);
      }
    }
    idxs[i-curr_idx+1] = num_wedges;
    if (2*(num_wedges+curr_num_wedges) > max_wedges) {
      if (i == curr_idx) {cout << "Space must accomodate seagulls originating from one edge\n"; exit(0);}
      return make_tuple(num_wedges, i, idxs);
    }
    num_wedges += curr_num_wedges;
  }
  return make_tuple(num_wedges, num_I, idxs);
}

template<class Sequence>
tuple<long, intT, long*> getNextIntersectWedgeIdx(uintE* eti, uintE* ite, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, 
  long max_wedges, intT curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (num_I - curr_idx < 10000) return getNextIntersectWedgeIdx_seq(eti, ite, I, num_I, GA, use_v, max_wedges, curr_idx);
  long* idxs = newA(long, num_I - curr_idx + 1);
  idxs[num_I-curr_idx] = 0;
  parallel_for(intT i=curr_idx; i < num_I; ++i) {
    idxs[i-curr_idx] = 0;
    uintE idx = I[i];
    uintE u = edgesV[idx];
    uintE v = edgesU[ite[idx]];
    intT v_offset = offsetsV[v];
    intT v_deg = offsetsV[v+1] - v_offset;
    intT u_deg = offsetsU[u+1] - offsetsU[u];
    for (intT k=0; k < v_deg; ++k) { 
      uintE u2 = edgesV[v_offset + k];
      if (u2 != u && u2 != UINT_E_MAX && (u2 != UINT_E_MAX - 1 || idx < v_offset + k)) {
        intT u2_deg = offsetsU[u2+1] - offsetsU[u2];
        idxs[i-curr_idx] += min(u_deg, u2_deg);
      }
    }
  }
  sequence::plusScan(idxs, idxs, num_I - curr_idx + 1);

  auto idx_map = make_in_imap<long>(num_I - curr_idx, [&] (size_t i) { return idxs[i+1]; });
  auto lte = [] (const long& l, const long& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx; //this rets first # > searched num
  long num = idxs[find_idx - curr_idx];
  //free(idxs);
  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }

  return make_tuple(num, find_idx, idxs); //TODO make sure right
}

template<class Sequence>
pair<long, intT> getIntersectWedges(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, Sequence I, intT num_I, bipartiteCSR& GA,
  bool use_v, long max_wedges, intT curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (curr_idx == 0){parallel_for(intT i=0; i < num_I; ++i){ current[I[i]] = 1; } }

  auto p = getNextIntersectWedgeIdx(eti, ite, I, num_I, GA, use_v, max_wedges, curr_idx);
  long num_wedges = _getIntersectWedges(ps, eti, ite, current, I, num_I, GA, use_v, get<0>(p), get<2>(p), 
    curr_idx, get<1>(p)); //TODO filter wedges_seq_tup
  free(get<2>(p));

  if (get<1>(p) == num_I) {
  parallel_for(intT i=0; i < num_I; ++i){
    edgesV[I[i]] = UINT_E_MAX; edgesU[ite[I[i]]] = UINT_E_MAX;
    current[I[i]] = 0;
  } }

  return make_pair(num_wedges, get<1>(p));
}

template <class E>
struct oneMaxF { 
  uintE* same;
  oneMaxF(uintE* _same) : same(_same) {}
  inline E operator() (const uintT& i) const {
    return (same[i] != UINT_E_MAX);
  }
};

void intersect_hash(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, bipartiteCSR& GA, bool use_v,
  uintE u, uintE v, uintE u2, uintE idx_vu, uintE idx_vu2) {

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  // have u, v, u2
  intT u2_offset = offsetsU[u2];
  intT u2_deg = offsetsU[u2+1] - u2_offset;

  intT u_offset = offsetsU[u];
  intT u_deg = offsetsU[u+1] - u_offset;
  bool* same = newA(bool, u2_deg);
  //auto u_map = make_in_imap<uintT>(u_deg, [&] (size_t i) { return edgesU[u_offset + i]; });
  //auto lt = [] (const uintE& l, const uintE& r) { return l < r; };
  //parallel_
  intT int_offset = 0;
  parallel_for(intT j = 0; j < u2_deg; ++j) { same[j] = 0; }
  for(intT j = 0; j < u2_deg; ++j) {
    uintE v2 = edgesU[u2_offset + j];
    uintE idx_v2u2 = eti[u2_offset + j];
    //if (v2 == UINT_E_MAX || v2 == v) same[j]=0;//same[j] = UINT_E_MAX;
    //else if (current[idx_v2u2] && idx_v2u2 < idx_vu) same[j] = 0; //v2 == UINT_E_MAX - 1
    //else {
    if(v2 != UINT_E_MAX && v2 != v && (!current[idx_v2u2] || idx_v2u2 >= idx_vu)) {
      while(int_offset < u_deg && (edgesU[u_offset + int_offset] == UINT_E_MAX || edgesU[u_offset + int_offset] < v2)) { int_offset++; }
      if (int_offset < u_deg && edgesU[u_offset+int_offset] == v2 && (!current[eti[u_offset + int_offset]] || idx_vu < eti[u_offset + int_offset])) {
        same[j] = 1;//edgesU[u2_offset + j];
        ps.update_hash.insert(make_pair(eti[u2_offset + j],1));
        ps.update_hash.insert(make_pair(eti[u_offset + int_offset],1));
      }
      else if(int_offset >= u_deg) break;
    /*intT find_idx = binary_search_mod(edgesU, u_offset, u_deg, v2);
    if (find_idx < u_deg && edgesU[u_offset + find_idx] == v2 &&
      (!current[eti[u_offset + find_idx]] || idx_vu < eti[u_offset + find_idx])) {
      same[j] = 1;//edgesU[u2_offset + j];
      ps.update_hash.insert(make_pair(eti[u2_offset + j],1));
      ps.update_hash.insert(make_pair(eti[u_offset + find_idx],1));
    }
    else same[j] = 0;//UINT_E_MAX;*/
    }
  }
  long num_same = sequence::sum(same, u2_deg);//sequence::reduce<long>((long) 0, u2_deg, addF<long>(), oneMaxF<long>(same));
  if (num_same > 0) { ps.update_hash.insert(make_pair(idx_vu2, num_same)); }
  free(same);
}

template<class Sequence>
void getIntersectWedgesHash(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, Sequence I, intT num_I, bipartiteCSR& GA,
  bool use_v) { //, long max_wedges, intT curr_idx
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  
  //TODO fix
  //if (curr_idx == 0) 
  //{parallel_for(intT i=0; i < num_I; ++i){ edgesV[I[i]] = UINT_E_MAX-1; edgesU[ite[I[i]]] = UINT_E_MAX-1; } }
  {parallel_for(intT i=0; i < num_I; ++i){ current[I[i]] = 1; } }

  //pair<long, intT> ret = getNextIntersectWedgeIdx(eti, ite, I, num_I, GA, use_v, max_wedges, curr_idx);
  auto wedges = ps.update_hash;
  parallel_for(intT i=0; i < num_I; ++i){
    // Set up for each active vertex
    uintE idx = I[i];
    uintE u = edgesV[idx];
    uintE v = edgesU[ite[idx]];
    intT v_offset = offsetsV[v];
    intT v_deg = offsetsV[v+1] - v_offset;
  // we want to go through all neighbors of v
  // intersect each N(u') with N(u)
  // get all u, v, u2
  // intersection of N(u) and N(u2) - 1 is the num butterflies to subtract from v, u2
    parallel_for (intT k=0; k < v_deg; ++k) { 
      uintE u2 = edgesV[v_offset + k];
      if (u2 != u && u2 != UINT_E_MAX && (!current[v_offset+k] || idx < v_offset + k)) {//u2 != UINT_E_MAX - 1
        intersect_hash(eti, ite, current, ps, GA, use_v, u, v, u2, idx, v_offset+k);
        //uintE tmp = intersect(&edgesU[u_offset], &edgesU[u2_offset], u_deg, u2_deg) - 1;
        //if (tmp > 0) {wedges.insert(make_pair(v_offset+k, tmp)); cout << "(" << u << ", " << v  <<", " << u2<< "), "; fflush(stdout);}
      }
	  }
  }
  //if (ret.second == num_I) {
  parallel_for(intT i=0; i < num_I; ++i){
    edgesV[I[i]] = UINT_E_MAX; edgesU[ite[I[i]]] = UINT_E_MAX;
    current[I[i]] = 0;
  } //}
  //return ret.second;
}


//***************************************************************************************************
//***************************************************************************************************

pair<tuple<uintE,uintE>*,long> _updateBuckets(uintE* update_idxs, long num_updates, uintE* butterflies, 
  array_imap<uintE> D, buckets<array_imap<uintE>> b, uintE k) {
  using X = tuple<uintE,uintE>;
  X* update = newA(X,num_updates);
  const intT eltsPerCacheLine = 64/sizeof(long);

  // Filter for bucket updates
  parallel_for(long i=0; i < num_updates; ++i) {
    const uintE u_idx = update_idxs[i];
    uintE old_b = D.s[u_idx];

    if (old_b > k) {
        uintE new_b = max(butterflies[eltsPerCacheLine*u_idx],k);
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

pair<tuple<uintE,uintE>*,long> _updateBuckets_seq(uintE* update_idxs, long num_updates, uintE* butterflies, 
  array_imap<uintE> D, buckets<array_imap<uintE>> b, uintE k) {
  using X = tuple<uintE, uintE>;
  X* update = newA(X, num_updates);
  const intT eltsPerCacheLine = 64/sizeof(long);
  long idx = 0;
  for(long i=0; i < num_updates; ++i) {
    uintE u_idx = update_idxs[i];
    uintE old_b = D.s[u_idx];

    if(old_b > k) {
      uintE new_b = max(butterflies[eltsPerCacheLine*u_idx], k);
      D.s[u_idx] = new_b;
      uintE new_bkt = b.get_bucket(old_b, new_b);
      update[idx] = make_tuple(u_idx, new_bkt);
      ++idx;
    }
  }
  return make_pair(update, idx);
}

void updateBuckets(array_imap<uintE>& D, buckets<array_imap<uintE>>& b, uintE k, uintE* butterflies, bool is_seq, long nu,
  uintE* update_idxs, long num_updates) {
  pair<tuple<uintE,uintE>*,long> bucket_pair;
  if (is_seq) bucket_pair = _updateBuckets_seq(update_idxs, num_updates, butterflies, D, b, k);
  else bucket_pair = _updateBuckets(update_idxs, num_updates, butterflies, D, b, k);

  vertexSubsetData<uintE> moved = vertexSubsetData<uintE>(nu, bucket_pair.second, bucket_pair.first);
  b.update_buckets(moved.get_fn_repr(), moved.size());

  moved.del();
}

#endif
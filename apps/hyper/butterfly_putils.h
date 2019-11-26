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

// Allocates and keeps all space needed for edge peeling algorithms, so that
// space can be reused between batches
struct PeelESpace {
  PeelType type;
  long nu;
  long stepSize;
  long n_side;
  _seq<tuple<uintE,long>> wedges_seq_tup;
  _seq<tuple<uintE,long>> wedges_seq_tup_fil;
  _seq<uintE> update_seq_int;
  sparseAdditiveSet<long, long> update_hash;
  _seq<uintE> wedges_seq_int;
  _seq<uintE> used_seq_int;

PeelESpace(PeelType _type, long _nu, long _stepSize, long _n_side) :
  type(_type), nu(_nu), stepSize(_stepSize), n_side(_n_side) {
    using X = tuple<uintE,long>;
    if (type != PBATCHS && type != PBATCHWA) update_seq_int = _seq<uintE>(newA(uintE, nu), nu);
    if (type == PHASH) update_hash = sparseAdditiveSet<long, long>(nu, (float) 1, LONG_MAX, LONG_MAX);
    else if (type == PHIST || type == PSORT) {
      wedges_seq_tup = _seq<X>(newA(X, nu), nu);
      wedges_seq_tup_fil = _seq<X>(newA(X, nu), nu);
    }
    else if (type == PBATCHS || type == PBATCHWA) {
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
    if (type != PBATCHS && type != PBATCHWA) update_seq_int.del();
    if (type == PHASH) update_hash.del();
    else if (type == PHIST || type == PSORT) { wedges_seq_tup.del(); wedges_seq_tup_fil.del(); }
    else if (type == PBATCHS || type == PBATCHWA) {wedges_seq_int.del(); used_seq_int.del();}
  }
};

// Allocates and keeps all space needed for vertex peeling algorithms, so that
// space can be reused between batches
struct PeelSpace {
  PeelType type;
  long nu;
  long stepSize;
  _seq<UVertexPair> wedges_seq_uvp;
  _seq<long> wedges_seq_int;
  _seq<long> wedges_seq_long;
  _seq<uintE> used_seq_int;
  _seq<uintE> update_seq_int;
  _seq<long> update_idx_seq_int;
  sparseAdditiveSet<long, long> update_hash;
  sparseAdditiveSet<long, long>* wedges_hash;
  sparseAdditiveSet<long, long>** wedges_hash_list;
  _seq<pair<long,long>> wedges_seq_intp;
  _seq<pair<long,long>> butterflies_seq_intp;
  long num_wedges_hash;
#ifndef OPENMP
  pbbsa::sequence<tuple<long, uintE>> tmp;
  pbbsa::sequence<tuple<long, uintE>> out;
#endif
PeelSpace(PeelType _type, long _nu, long _stepSize) : type(_type), nu(_nu), stepSize(_stepSize) {
  using E = pair<long, long>;
  using X = pair<uintE,long>;
  update_seq_int = _seq<uintE>(newA(uintE, nu), nu);

  if (type == PHASH) {
    using T = sparseAdditiveSet<long, long>*;
    wedges_hash = new sparseAdditiveSet<long, long>(1,1, LONG_MAX, LONG_MAX);
    wedges_hash_list = newA(T, 1);
    wedges_hash_list[0] = wedges_hash;
    num_wedges_hash = 1;
    update_hash = sparseAdditiveSet<long, long>(nu, (float) 1, LONG_MAX, LONG_MAX);
    wedges_seq_intp = _seq<E>(newA(E, nu), nu);
    butterflies_seq_intp = _seq<E>(newA(E, nu), nu);
  }
  else if (type == PSORT) wedges_seq_uvp = _seq<UVertexPair>(newA(UVertexPair, nu), nu);
#ifndef OPENMP
  else if (type == PHIST) {
    wedges_seq_long = _seq<long>(newA(long, nu), nu);
    tmp = pbbsa::sequence<tuple<long, uintE>>();
    out = pbbsa::sequence<tuple<long, uintE>>();
  }
#endif
  else {
    wedges_seq_int = _seq<long>(newA(long, nu*stepSize), nu*stepSize);
    granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges_seq_int.A[i] = 0; });
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
    if (type != PHASH) return nullptr;
    
    size_t find_idx = log2RoundUp(size);
    if (find_idx < num_wedges_hash) {
      wedges_hash = wedges_hash_list[find_idx];
      wedges_hash->clear(); 

      return wedges_hash_list[find_idx];
    }
    using T = sparseAdditiveSet<long, long>*;
    T* new_wedges_hash_list = newA(T, find_idx+1);
    parallel_for(long i=0; i < num_wedges_hash; ++i) {
      new_wedges_hash_list[i] = wedges_hash_list[i];
    }
    parallel_for(long i=num_wedges_hash; i < find_idx+1; ++i) {
      new_wedges_hash_list[i] = new sparseAdditiveSet<long, long>(1u << i, 1, LONG_MAX, LONG_MAX);
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
    if (type == PHASH) update_hash.clear();
  }

  void del() {
    update_seq_int.del();
    if (type == PHASH) {
      parallel_for(long i=0; i < num_wedges_hash; ++i) { wedges_hash_list[i]->del(); free(wedges_hash_list[i]); }
      free(wedges_hash_list);
      update_hash.del();
      wedges_seq_intp.del(); 
      butterflies_seq_intp.del();
    }
    else if (type == PSORT) wedges_seq_uvp.del();
    else if (type == PHIST) wedges_seq_long.del();
    else { wedges_seq_int.del(); used_seq_int.del(); }
  }
};


//***************************************************************************************************
//***************************************************************************************************

/*
 *  Retrieve all wedges in this batch based on active vertices, sequentially.
 * 
 *  wedges_seq: Array to store wedges
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  cons      : Constructor for a wedge; takes as input i, u2, v, j, and k, where i and u2 are endpoints, v is the
 *              center, and j and k are indices indicating the edge indices for (i, v) and (v, u2) respectively.
 *  num_wedges: Number of wedges in this batch
 *  curr_idx  : Starting index in I of this batch
 *  next_idx  : Ending index in I of this batch
 */
template<class wedge, class wedgeCons, class Sequence>
  void _getActiveWedges_seq(_seq<wedge>& wedges_seq, Sequence I, long num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons,
			    long num_wedges, long curr_idx, long next_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long idx = 0;
  // Store wedges sequentially
  // Iterate through each index i in this batch
  for(long i=curr_idx; i < next_idx; ++i) {
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    for(long j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find two-hop neighbors of u
      for (long k = 0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset + k];
        if (u2 != u) {
          wedges_seq.A[idx] = cons(u, u2, v, j, k);
          ++idx;
        }
      }
    }
  }
}

/*
 *  Retrieve all wedges in this batch based on active vertices, in parallel.
 * 
 *  wedges_seq: Array to store wedges
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  cons      : Constructor for a wedge; takes as input i, u2, v, j, and k, where i and u2 are endpoints, v is the
 *              center, and j and k are indices indicating the edge indices for (i, v) and (v, u2) respectively.
 *  num_wedges: Number of wedges in this batch
 *  curr_idx  : Starting index in I of this batch
 *  next_idx  : Ending index in I of this batch
 */
template<class wedge, class wedgeCons, class Sequence>
  void _getActiveWedges(_seq<wedge>& wedges_seq, Sequence I,long num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons,
			long num_wedges, long curr_idx=0, long next_idx=LONG_MAX) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Allocate space for wedge storage
  if (wedges_seq.n < num_wedges) {
    free(wedges_seq.A);
    wedges_seq.A = newA(wedge, num_wedges);
    wedges_seq.n = num_wedges;
  }

  // Set next index in batch
  if (next_idx == LONG_MAX) next_idx = num_I;
  // Run sequentially if batch is small enough
  if (num_wedges < 10000) return _getActiveWedges_seq<wedge>(wedges_seq, I, num_I, GA, use_v, cons, num_wedges,
                                                             curr_idx, next_idx);

  // First, we must retrive the indices at which to store each wedge (so that storage can be parallelized)
  // Array of indices associated with wedges for each active vertex
  long* sg_idxs = newA(long, next_idx - curr_idx + 1);
  sg_idxs[next_idx-curr_idx] = 0;

  // 2D array of indices associated with wedges for each neighbor of each active vertex
  using T = long*;
  T* nbhd_idxs = newA(T, next_idx-curr_idx); 

  parallel_for(long i=curr_idx; i < next_idx; ++i) {
    // Set up for each active vertex
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    long u_deg = offsetsU[u+1] - u_offset;
    // Allocate space for indices associated with each neighbor of u
    nbhd_idxs[i-curr_idx] = newA(long, u_deg + 1);
    (nbhd_idxs[i-curr_idx])[u_deg] = 0;
    // Assign indices associated with each neighbor of u
    parallel_for(long j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      (nbhd_idxs[i-curr_idx])[j] = offsetsV[v+1] - offsetsV[v] - 1;
    }
    sequence::plusScan(nbhd_idxs[i-curr_idx], nbhd_idxs[i-curr_idx], u_deg+1);
    // Set up indices associated with u
    sg_idxs[i-curr_idx] = (nbhd_idxs[i-curr_idx])[u_deg];
  }
  // Assign indices associated with each active vertex
  sequence::plusScan(sg_idxs, sg_idxs, next_idx-curr_idx + 1);

  // Store wedges in parallel
  parallel_for(long i=curr_idx; i < next_idx; ++i) {
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
  parallel_for(long i=curr_idx; i<next_idx;++i) { 
    free(nbhd_idxs[i-curr_idx]); 
  }
  free(nbhd_idxs);
  free(sg_idxs);
}

/*
 *  Retrieve and hash all wedges in this batch based on active vertices.
 * 
 *  ps        : Holds all array space needed, to be reused between buckets
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  cons      : Constructor for a wedge; takes as input i, u2, v, j, and k, where i and u2 are endpoints, v is the
 *              center, and j and k are indices indicating the edge indices for (i, v) and (v, u2) respectively.
 *  num_wedges: Number of wedges in this batch
 *  curr_idx  : Starting index in I of this batch
 *  next_idx  : Ending index in I of this batch
 */
template<class wedgeCons, class Sequence>
  void _getActiveWedgesHash(PeelSpace& ps, Sequence I, long num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons,
			    long num_wedges, long curr_idx=0, long next_idx=INT_T_MAX) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Set next index in batch
  if (next_idx == INT_T_MAX) next_idx = num_I;

  // Allocate hash table space
  auto wedges = ps.resize(num_wedges);

  parallel_for(long i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    parallel_for (long j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find all seagulls with center v and endpoint u
      for (long k=0; k < v_deg; ++k) { 
        uintE u2 = edgesV[v_offset + k];
        if (u2 != u) wedges->insert(make_pair(cons(u, u2, v, j, k),1));
      }
    }
  }
}
/*
 *  Identify the last index of the batch of wedges to be processed,
 *  sequentially.
 * 
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges: Maximum number of wedges in a batch
 *  curr_idx  : Starting index of this batch
 * 
 *  Returns the index starting the next batch.
 */
template<class Sequence>
pair<long,intT> getNextActiveWedgeIdx_seq(Sequence I, long num_I, bipartiteCSR& GA, bool use_v, long max_wedges,
                                          long curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long num_wedges = 0;
  // Iterate through each index i starting with curr_idx
  for(long i=curr_idx; i < num_I; ++i) {
    uintE u = I[i];
    intT u_offset = offsetsU[u];
    intT u_deg = offsetsU[u+1] - u_offset;
    long num = 0;
    for(long j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      intT v_deg = offsetsV[v+1] - offsetsV[v];
      // Determine the number of wedges that u contributes to the batch
      num += (v_deg - 1);
      // If we have exceeded the max batch size, end the batch here
      if (num_wedges + num > max_wedges) {
        if (i == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
        return make_pair(num_wedges, i);
      }
    }
    num_wedges += num;
  }
  return make_pair(num_wedges, num_I);
}

/*
 *  Compute number of wedges that each vertex in an active set contributes,
 *  so that wedges can be properly batched and processed in parallel.
 * 
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  curr_idx  : Starting index of this batch
 * 
 *  Returns the index starting the next batch.
 */
template<class Sequence>
long* getActiveWedgeIdxs(Sequence I, long num_I, bipartiteCSR& GA, bool use_v, long curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long* idxs = newA(long, num_I - curr_idx + 1);
  idxs[num_I-curr_idx] = 0;

  // Determine the number of wedges each vertex u contributes
  parallel_for(long i=curr_idx; i < num_I; ++i) {
    idxs[i-curr_idx] = 0;
    uintE u = I[i];
    uintT u_offset = offsetsU[u];
    uintT u_deg = offsetsU[u+1] - u_offset;
    for(long j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset + j];
      uintT v_deg = offsetsV[v+1] - offsetsV[v];
      idxs[i-curr_idx] += (v_deg - 1);
    }
  }
  sequence::plusScan(idxs, idxs, num_I - curr_idx + 1);

  return idxs;
}

/*
 *  Identify the last index of the batch of wedges to be processed,
 *  in parallel.
 * 
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges: Maximum number of wedges in a batch
 *  curr_idx  : Starting index of this batch
 * 
 *  Returns the index starting the next batch.
 */
template<class Sequence>
pair<long, long> getNextActiveWedgeIdx(Sequence I, long num_I, bipartiteCSR& GA, bool use_v, long max_wedges,
                                       long curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Determine if the remaining number of vertices is small enough to process sequentially
  if (num_I - curr_idx < 10000) return getNextActiveWedgeIdx_seq(I, num_I, GA, use_v, max_wedges, curr_idx);
  // Retrieve number of wedges on each index
  long* idxs = getActiveWedgeIdxs(I, num_I, GA, use_v, curr_idx);

  // Binary search for the first index that exceeds the max number of wedges in a batch, using idxs
  auto idx_map = make_in_imap<long>(num_I - curr_idx, [&] (size_t i) { return idxs[i+1]; });
  auto lte = [] (const long& l, const long& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx;
  long num = idxs[find_idx - curr_idx];
  free(idxs);

  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }

  return make_pair(num, find_idx);
}

/*
 *  Retrieve all wedges in this batch corresponding with active vertices in I.
 * 
 *  wedges_seq: Array to store wedges
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  cons      : Constructor for a wedge; takes as input i, u2, v, j, and k, where i and u2 are endpoints, v is the
 *              center, and j and k are indices indicating the edge indices for (i, v) and (v, u2) respectively.
 *  max_wedges: Max number of wedges in a batch
 *  curr_idx  : Starting index of this batch
 *  num_wedges: Number of wedges left in the graph total
 * 
 *  Returns a pair, the first of which indicates the number of wedges in this
 *  batch and the second of which indicates the index starting the next
 *  batch.
 */
template<class wedge, class wedgeCons, class Sequence>
  pair<long, long> getActiveWedges(_seq<wedge>& wedges_seq, Sequence I, long num_I, bipartiteCSR& GA, bool use_v,
				   wedgeCons cons, long max_wedges, long curr_idx, long num_wedges) {
  if (max_wedges >= num_wedges) {
    _getActiveWedges<wedge>(wedges_seq, I, num_I, GA, use_v, cons, num_wedges);
    return make_pair(num_wedges, num_I);
  }

  auto p = getNextActiveWedgeIdx(I, num_I, GA, use_v, max_wedges, curr_idx);
  _getActiveWedges<wedge>(wedges_seq, I, num_I, GA, use_v, cons, p.first, curr_idx, p.second);
  return p;
}

/*
 *  Retrieve and hash all wedges in this batch corresponding with active
 *  vertices in I.
 * 
 *  ps        : Holds all array space needed, to be reused between buckets
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  cons      : Constructor for a wedge; takes as input i, u2, v, j, and k, where i and u2 are endpoints, v is the
 *              center, and j and k are indices indicating the edge indices for (i, v) and (v, u2) respectively.
 *  max_wedges: Max number of wedges in a batch
 *  curr_idx  : Starting index of this batch
 *  num_wedges: Number of wedges left in the graph total
 * 
 *  Returns a pair, the first of which indicates the number of wedges in this
 *  batch and the second of which indicates the index starting the next
 *  batch.
 */
template<class wedgeCons, class Sequence>
  long getActiveWedgesHash(PeelSpace& ps, Sequence I, long num_I, bipartiteCSR& GA, bool use_v, wedgeCons cons,
			   long max_wedges, long curr_idx, long num_wedges) {
  if (max_wedges >= num_wedges) {
    _getActiveWedgesHash(ps, I, num_I, GA, use_v, cons, num_wedges);
    return num_I;
  }
  auto p = getNextActiveWedgeIdx(I, num_I, GA, use_v, max_wedges, curr_idx);
  _getActiveWedgesHash(ps, I, num_I, GA, use_v, cons, p.first, curr_idx, p.second);
  return p.second;
}


//***************************************************************************************************
//***************************************************************************************************

/*
 *  Given a wedge, find the intersection of the neighborhoods of the two
 *  endpoints u and u2. This gives updated butterfly counts per edge
 *  in the discovered butterflies.
 * 
 *  wedges_seq: Array to store edges with updated butterfly counts
 *  eti       : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite       : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current   : Array to keep track of active edges
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  u         : Endpoint of wedge
 *  v         : Center of wedge
 *  u2        : Endpoint of wedge
 *  idx_vu    : Index of u in the adjacency list for v
 *  idx_vu2   : Index of u2 in the adjacency list for v
 *  wedges_idx: Index in the wedges_seq array to add butterfly counts to
 */
void intersect(_seq<tuple<uintE,long>>& wedges_seq, uintE* eti, uintE* ite, bool* current, 
	       bipartiteCSR& GA, bool use_v, uintE u, uintE v, uintE u2, uintE idx_vu, uintE idx_vu2, intT wedges_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  intT u2_offset = offsetsU[u2];
  intT u2_deg = offsetsU[u2+1] - u2_offset;

  intT u_offset = offsetsU[u];
  intT u_deg = offsetsU[u+1] - u_offset;

  // Initialize array to find intersection between N(u) and N(u2)
  bool* same = newA(bool, u2_deg);
  intT int_offset = 0;
  parallel_for(intT j = 0; j < u2_deg; ++j) { same[j] = 0; }

  // Iterate through all neighbors v2 of u2
  for(intT j = 0; j < u2_deg; ++j) {
    uintE v2 = edgesU[u2_offset + j];
    uintE idx_v2u2 = eti[u2_offset + j];
    if(v2 != UINT_E_MAX && v2 != v && (!current[idx_v2u2] || idx_v2u2 >= idx_vu)) {
      // Find if v2 is also a neighbor of u
      while(int_offset < u_deg && (edgesU[u_offset + int_offset] == UINT_E_MAX || edgesU[u_offset + int_offset] < v2)) {
        int_offset++;
      }
      if (int_offset < u_deg && edgesU[u_offset+int_offset] == v2 &&
	  (!current[eti[u_offset + int_offset]] || idx_vu < eti[u_offset + int_offset])) {
        // Note the butterfly on (v2, u2) and (v2, u)
        same[j] = 1;
        wedges_seq.A[wedges_idx++] = make_tuple(eti[u2_offset + j], (long) 1);
        wedges_seq.A[wedges_idx++] = make_tuple(eti[u_offset + int_offset], (long) 1);
      }
      else if(int_offset >= u_deg) break;
    }
  }
  // Note the butterflies on (v, u2)
  long num_same = sequence::sum(same, u2_deg);
  if (num_same > 0) { wedges_seq.A[wedges_idx++] = make_tuple(idx_vu2, num_same); }
  free(same);
}

/*
 *  Find all wedges in this batch and use intersect to obtain the requisite
 *  butterfly updates per edge.
 * 
 *  ps        : Holds all array space needed, to be reused between buckets
 *  eti       : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite       : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current   : Array to keep track of active edges
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges: Number of wedges in this batch
 *  idxs      : Indices of butterfly updates in this batch so that they can be accessed in parallel
 *  curr_idx  : Starting index in I of this batch
 *  next_idx  : Ending index in I of this batch
 */
template<class Sequence>
long _getIntersectWedges(PeelESpace& ps, uintE* eti, uintE* ite, bool* current, Sequence I, intT num_I, 
                         bipartiteCSR& GA, bool use_v, long num_wedges, long* idxs, intT curr_idx=0,
                         intT next_idx=INT_T_MAX) {
  using X = tuple<uintE,long>;

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Set up space to store counts
  if (ps.wedges_seq_tup.n < num_wedges*2+next_idx-curr_idx) {
    long num = num_wedges*2+next_idx-curr_idx;
    free(ps.wedges_seq_tup.A);
    free(ps.wedges_seq_tup_fil.A);
    ps.wedges_seq_tup.A = newA(X, num);
    ps.wedges_seq_tup_fil.A = newA(X, num);
    ps.wedges_seq_tup.n = num;
    ps.wedges_seq_tup_fil.n = num;
  }
  parallel_for(long i=0; i < num_wedges*2+next_idx-curr_idx; ++i) {
    ps.wedges_seq_tup.A[i] = make_tuple(UINT_E_MAX, LONG_MAX);
  }
  auto nonMaxLF = [] (const X& r) { return get<0>(r) != UINT_E_MAX || get<1>(r) != LONG_MAX; };

  // Set up next index in batch
  if (next_idx == INT_T_MAX) next_idx = num_I;

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    uintE idx = I[i];
    uintE u = edgesV[idx];
    uintE v = edgesU[ite[idx]];
    intT v_offset = offsetsV[v];
    intT v_deg = offsetsV[v+1] - v_offset;
    // Iterate through all neighbors u2 of v
    parallel_for (intT k=0; k < v_deg; ++k) { 
      uintE u2 = edgesV[v_offset + k];
      // Find the intersection of N(u2) and N(u), and store the updated butterfly counts
      if (u2 != u && u2 != UINT_E_MAX && (!current[v_offset+k] || idx < v_offset + k)) {
        intersect(ps.wedges_seq_tup, eti, ite, current, GA, use_v, u, v, u2, idx, v_offset + k,
                  idxs[i-curr_idx]*2 + i - curr_idx);
      }
    }
  }
  // Filter out unused entries
  long num_wedges_fil = sequence::filter(ps.wedges_seq_tup.A, ps.wedges_seq_tup_fil.A,
                                         num_wedges * 2 + next_idx - curr_idx, nonMaxLF);
  return num_wedges_fil;
}

/*
 *  Identify the last index of the batch of wedges to be processed,
 *  sequentially.
 * 
 *  eti       : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite       : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges: Maximum number of wedges in a batch
 *  curr_idx  : Starting index of this batch
 * 
 *  Returns the number of wedges in this batch, the index in I starting the
 *  next batch, and the indices of butterfly updates in this batch so that they
 *  can be accessed in parallel.
 */
template<class Sequence>
tuple<long, intT, long*> getNextIntersectWedgeIdx_seq(uintE* eti, uintE* ite, Sequence I, intT num_I, bipartiteCSR& GA,
                                                      bool use_v, long max_wedges, intT curr_idx) {
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
    // Iterate through all neighbors u2 of v
    for (intT k=0; k < v_deg; ++k) {
      uintE u2 = edgesV[v_offset + k];
      // Find the max amount of space that updates from processing wedge (u, v, u2) could take
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

/*
 *  Computes number of butterfly updates that each vertex in an active set
 *  contributes, so that updates can be properly batched and processed in
 *  parallel.
 * 
 *  ite       : Array that converts edge indices from the other orientation to the first
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  curr_idx  : Starting index of this batch
 * 
 *  Returns the indices of butterfly updates in this batch so that they can
 *  be accessed in parallel.
 */
template<class Sequence>
long* getIntersectWedgeIdxs(uintE* ite, Sequence I, intT num_I, bipartiteCSR& GA, bool use_v, intT curr_idx=0) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Initialize butterfly update indices
  long* idxs = newA(long, num_I - curr_idx + 1);
  idxs[num_I-curr_idx] = 0;

  // Iterate through all vertices in active set, starting from index curr_idx
  parallel_for(intT i=curr_idx; i < num_I; ++i) {
    idxs[i-curr_idx] = 0;
    uintE idx = I[i];
    uintE u = edgesV[idx];
    uintE v = edgesU[ite[idx]];
    intT v_offset = offsetsV[v];
    intT v_deg = offsetsV[v+1] - v_offset;
    intT u_deg = offsetsU[u+1] - offsetsU[u];
    // Iterate through all neighbors u2 of v
    for (intT k=0; k < v_deg; ++k) { 
      uintE u2 = edgesV[v_offset + k];
      // Find the max amount of space that updates from processing wedge (u, v, u2) could take
      if (u2 != u && u2 != UINT_E_MAX && (u2 != UINT_E_MAX - 1 || idx < v_offset + k)) {
        intT u2_deg = offsetsU[u2+1] - offsetsU[u2];
        idxs[i-curr_idx] += min(u_deg, u2_deg);
      }
    }
  }
  sequence::plusScan(idxs, idxs, num_I - curr_idx + 1);

  return idxs;
}

/*
 *  Identify the last index of the batch of wedges to be processed,
 *  in parallel.
 * 
 *  eti       : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite       : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges: Maximum number of wedges in a batch
 *  curr_idx  : Starting index of this batch
 * 
 *  Returns the number of wedges in this batch, the index in I starting the
 *  next batch, and the indices of butterfly updates in this batch so that they
 *  can be accessed in parallel.
 */
template<class Sequence>
tuple<long, intT, long*> getNextIntersectWedgeIdx(uintE* eti, uintE* ite, Sequence I, intT num_I, bipartiteCSR& GA,
                                                  bool use_v, long max_wedges, intT curr_idx) {
  // Determine if the remaining number of vertices is small enough to process sequentially
  if (num_I - curr_idx < 10000)
    return getNextIntersectWedgeIdx_seq(eti, ite, I, num_I, GA, use_v, max_wedges, curr_idx);

  // Find indices of wedges in this batch so that they can be accessed in parallel
  long* idxs = getIntersectWedgeIdxs(ite, I, num_I, GA, use_v, curr_idx);

  // Binary search for the first index that exceeds the max number of wedges in a batch, using idxs
  auto idx_map = make_in_imap<long>(num_I - curr_idx, [&] (size_t i) { return idxs[i+1]; });
  auto lte = [] (const long& l, const long& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx;
  long num = idxs[find_idx - curr_idx];

  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }

  return make_tuple(num, find_idx, idxs);
}

/*
 *  Retrieve all butterfly updates in this batch corresponding with active
 *  vertices in I.
 *  
 *  eti       : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite       : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current   : Array to keep track of active edges
 *  ps        : Holds all array space needed, to be reused between buckets
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  max_wedges: Max number of wedges in a batch
 *  curr_idx  : Starting index of this batch
 * 
 *  Returns a pair, the first of which indicates the number of updates in this
 *  batch and the second of which indicates the index starting the next
 *  batch.
 */
template<class Sequence>
pair<long, intT> getIntersectWedges(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, Sequence I, intT num_I,
                                    bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (curr_idx == 0){parallel_for(intT i=0; i < num_I; ++i){ current[I[i]] = 1; } }

  auto p = getNextIntersectWedgeIdx(eti, ite, I, num_I, GA, use_v, max_wedges, curr_idx);
  long num_wedges = _getIntersectWedges(ps, eti, ite, current, I, num_I, GA, use_v, get<0>(p), get<2>(p), 
					curr_idx, get<1>(p));
  free(get<2>(p));

  if (get<1>(p) == num_I) {
    parallel_for(intT i=0; i < num_I; ++i){
      edgesV[I[i]] = UINT_E_MAX; edgesU[ite[I[i]]] = UINT_E_MAX;
      current[I[i]] = 0;
    } }

  return make_pair(num_wedges, get<1>(p));
}

// Function that returns true if an element of an array is not max
template <class E>
struct oneMaxF { 
  uintE* same;
oneMaxF(uintE* _same) : same(_same) {}
  inline E operator() (const uintT& i) const {
    return (same[i] != UINT_E_MAX);
  }
};

/*
 *  Given a wedge, find the intersection of the neighborhoods of the two
 *  endpoints u and u2. This gives updated butterfly counts per edge
 *  in the discovered butterflies. Store these updated counts in a hash table.
 * 
 *  eti       : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite       : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current   : Array to keep track of active edges
 *  ps        : Holds all array space needed, to be reused between buckets
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  u         : Endpoint of wedge
 *  v         : Center of wedge
 *  u2        : Endpoint of wedge
 *  idx_vu    : Index of u in the adjacency list for v
 *  idx_vu2   : Index of u2 in the adjacency list for v
 */
void intersect_hash(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, bipartiteCSR& GA, bool use_v,
                    uintE u, uintE v, uintE u2, uintE idx_vu, uintE idx_vu2) {

  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  intT u2_offset = offsetsU[u2];
  intT u2_deg = offsetsU[u2+1] - u2_offset;

  intT u_offset = offsetsU[u];
  intT u_deg = offsetsU[u+1] - u_offset;
  bool* same = newA(bool, u2_deg);

  intT int_offset = 0;
  // Iterate through all neighbors v2 of u2
  parallel_for(intT j = 0; j < u2_deg; ++j) { same[j] = 0; }
  for(intT j = 0; j < u2_deg; ++j) {
    uintE v2 = edgesU[u2_offset + j];
    uintE idx_v2u2 = eti[u2_offset + j];
    // Find if v2 is also a neighbor of u
    if(v2 != UINT_E_MAX && v2 != v && (!current[idx_v2u2] || idx_v2u2 >= idx_vu)) {
      while(int_offset < u_deg && (edgesU[u_offset + int_offset] == UINT_E_MAX || edgesU[u_offset + int_offset] < v2)) {
        int_offset++;
      }
      if (int_offset < u_deg && edgesU[u_offset+int_offset] == v2 &&
	  (!current[eti[u_offset + int_offset]] || idx_vu < eti[u_offset + int_offset])) {
        // Note the butterfly on (v2, u2) and (v2, u)
        same[j] = 1;
        ps.update_hash.insert(make_pair(eti[u2_offset + j],1));
        ps.update_hash.insert(make_pair(eti[u_offset + int_offset],1));
      }
      else if(int_offset >= u_deg) break;
    }
  }

  // Note the butterfly on (v, u2)
  long num_same = sequence::sum(same, u2_deg);
  if (num_same > 0) { ps.update_hash.insert(make_pair(idx_vu2, num_same)); }
  free(same);
}

/*
 *  Find all wedges in this batch and use intersect_hash to obtain the requisite
 *  butterfly updates per edge.
 * 
 *  eti       : Array that converts edge indices from one orientation (as given by use_v) to the other
 *  ite       : Array that converts edge indices from the other orientation to the first, opposite of eti
 *  current   : Array to keep track of active edges
 *  ps        : Holds all array space needed, to be reused between buckets
 *  I         : Array that holds active vertices
 *  num_I     : Length of I
 *  GA        : Bipartite graph in CSR format
 *  use_v     : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *              produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 */
template<class Sequence>
void getIntersectWedgesHash(uintE* eti, uintE* ite, bool* current, PeelESpace& ps, Sequence I, intT num_I,
                            bipartiteCSR& GA, bool use_v) {
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  
  {parallel_for(intT i=0; i < num_I; ++i){ current[I[i]] = 1; } }

  auto wedges = ps.update_hash;
  parallel_for(intT i=0; i < num_I; ++i){
    // Set up for each active vertex
    uintE idx = I[i];
    uintE u = edgesV[idx];
    uintE v = edgesU[ite[idx]];
    intT v_offset = offsetsV[v];
    intT v_deg = offsetsV[v+1] - v_offset;
    // Iterate through each neighbor u2 of v
    parallel_for (intT k=0; k < v_deg; ++k) { 
      uintE u2 = edgesV[v_offset + k];
      // Find the intersection of N(u2) and N(u), and save the updated butterfly count in a hash table
      if (u2 != u && u2 != UINT_E_MAX && (!current[v_offset+k] || idx < v_offset + k)) {
        intersect_hash(eti, ite, current, ps, GA, use_v, u, v, u2, idx, v_offset+k);
      }
    }
  }

  // Reset current active set
  parallel_for(intT i=0; i < num_I; ++i){
    edgesV[I[i]] = UINT_E_MAX; edgesU[ite[I[i]]] = UINT_E_MAX;
    current[I[i]] = 0;
  }
}


//***************************************************************************************************
//***************************************************************************************************

/*
 *  Given vertices/edges with updated butterfly counts, find updated buckets
 *  in parallel.
 *  
 *  update_idxs: Vertices/edges with updated butterfly counts
 *  num_updates: Length of update_idxs
 *  butterflies: Butterfly counts on vertices/edges
 *  D          : Map of vertices to the buckets that they're in
 *  b          : Bucketing structure
 *  k          : ID of current bucket
 * 
 *  Returns an array with vertices/edges and their corresponding updated buckets,
 *  and the length of that array.
 */
pair<tuple<uintE,long>*,long> _updateBuckets(uintE* update_idxs, long num_updates, long* butterflies, 
                                             array_imap<long> D, buckets<array_imap<long>> b, long k) {
  using X = tuple<uintE,long>;
  X* update = newA(X,num_updates);
  const long eltsPerCacheLine = 64/sizeof(long);

  // Iterate through each vertex/edge with updated butterfly counts
  parallel_for(long i=0; i < num_updates; ++i) {
    const uintE u_idx = update_idxs[i];
    long old_b = D.s[u_idx];
    // Find updated buckets
    if (old_b > k) {
      long new_b = max(butterflies[eltsPerCacheLine*u_idx],k);
      D.s[u_idx] = new_b;
      long new_bkt = b.get_bucket(old_b, new_b);
      update[i] = make_tuple(u_idx, new_bkt);
    }
    else { update[i] = make_tuple(UINT_E_MAX, LONG_MAX); }
  }

  X* update_filter = newA(X, num_updates);
  auto nonMaxLF = [] (const X& r) { return get<0>(r) != UINT_E_MAX || get<1>(r) != LONG_MAX; };
  long num_updates_filter = sequence::filter(update, update_filter, num_updates, nonMaxLF);
  free(update);

  return make_pair(update_filter,  num_updates_filter);
}

/*
 *  Given vertices/edges with updated butterfly counts, find updated buckets
 *  sequentially.
 *  
 *  update_idxs: Vertices/edges with updated butterfly counts
 *  num_updates: Length of update_idxs
 *  butterflies: Butterfly counts on vertices/edges
 *  D          : Map of vertices to the buckets that they're in
 *  b          : Bucketing structure
 *  k          : ID of current bucket
 * 
 *  Returns an array with vertices/edges and their corresponding updated buckets,
 *  and the length of that array.
 */
pair<tuple<uintE,long>*,long> _updateBuckets_seq(uintE* update_idxs, long num_updates, long* butterflies, 
                                                 array_imap<long> D, buckets<array_imap<long>> b, long k) {
  using X = tuple<uintE, long>;
  X* update = newA(X, num_updates);
  const long eltsPerCacheLine = 64/sizeof(long);
  long idx = 0;

  // Iterate through each vertex/edge with updated butterfly counts
  for(long i=0; i < num_updates; ++i) {
    uintE u_idx = update_idxs[i];
    long old_b = D.s[u_idx];
    // Find updated buckets
    if(old_b > k) {
      long new_b = max(butterflies[eltsPerCacheLine*u_idx], k);
      D.s[u_idx] = new_b;
      long new_bkt = b.get_bucket(old_b, new_b);
      update[idx] = make_tuple(u_idx, new_bkt);
      ++idx;
    }
  }
  return make_pair(update, idx);
}

/*
 *  Given vertices/edges with updated butterfly counts, update the bucketing
 *  structure.
 *  
 *  D          : Map of vertices to the buckets that they're in
 *  b          : Bucketing structure
 *  k          : ID of current bucket
 *  butterflies: Butterfly counts on vertices/edges
 *  is_seq     : Bool indicating if there are enough buckets to be processed in parallel
 *  nu         : Number of vertices/edges that are being bucketed
 *  update_idxs: Vertices/edges with updated butterfly counts
 *  num_updates: Length of update_idxs
 */
void updateBuckets(array_imap<long>& D, buckets<array_imap<long>>& b, uintE k, long* butterflies, bool is_seq, long nu,
                   uintE* update_idxs, long num_updates) {
  pair<tuple<uintE,long>*,long> bucket_pair;
  if (is_seq) bucket_pair = _updateBuckets_seq(update_idxs, num_updates, butterflies, D, b, k);
  else bucket_pair = _updateBuckets(update_idxs, num_updates, butterflies, D, b, k);

  vertexSubsetData<long> moved = vertexSubsetData<long>(nu, bucket_pair.second, bucket_pair.first);
  b.update_buckets(moved.get_fn_repr(), moved.size());

  moved.del();
}

#endif

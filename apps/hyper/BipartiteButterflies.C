#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#define HYPER 1

#include "hygra.h"
#include "index_map.h"
#include "bucket.h"
#include "edgeMapReduce.h"

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "hypergraphIO.h"
#include "parseCommandLine.h"
#include "vertex.h"
#include "sequence.h"
#include "sparseSet.h"
#include "sampleSort.h"
#include "../../lib/histogram.h"
#include "../../lib/gbbs-histogram.h"

#include "butterfly_count.h"

using namespace std;

/*
 *  Retrieves all seagulls associated with the vertices in active (paths of length two with an endpoint
 *  in active), and stores the two endpoints using the constructor cons. The vertices in active are
 *  assumed to be a subset of U, and V represents the other bipartition of the bipartite graph.
 * 
 *  active: Set of active vertices (subset of U)
 *  U     : One bipartition of vertices
 *  V     : The other bipartition of vertices
 *  cons  : Constructor to create seagulls, given two endpoints (where the first is the active vertex, 
 *          and the second is the other endpoint)
 * 
 *  Returns: Array of seagulls and the number of seagulls
 */
template<class seagull, class vertex, class seagullCons>
pair<seagull*, long> getSeagulls(vertexSubset active, vertex* V, vertex* U, seagullCons cons) {
  // First, we must retrive the indices at which to store each seagull (so that storage can be parallelized)
  // Array of indices associated with seagulls for each active vertex
  long* sg_idxs = newA(long, active.size() + 1);
  sg_idxs[active.size()] = 0;

  // 2D array of indices associated with seagulls for each neighbor of each active vertex
  using T = long*;
  T* nbhd_idxs = newA(T, active.size()); 

  parallel_for(long i=0; i < active.size(); ++i) {
    // Set up for each active vertex
    uintE u_idx = active.vtx(i);
    const vertex u = U[u_idx];
    const uintE u_deg = u.getOutDegree();
    // Allocate space for indices associated with each neighbor of u
    nbhd_idxs[i] = newA(long, u_deg + 1);
    (nbhd_idxs[i])[u_deg] = 0;
    // Assign indices associated with each neighbor of u
    parallel_for(long j=0; j < u_deg; ++j) {
      (nbhd_idxs[i])[j] = V[u.getOutNeighbor(j)].getOutDegree()-1;
    }
    sequence::plusScan(nbhd_idxs[i], nbhd_idxs[i], u_deg+1);
    // Set up indices associated with u
    sg_idxs[i] = (nbhd_idxs[i])[u_deg];
  }
  // Assign indices associated with each active vertex
  sequence::plusScan(sg_idxs, sg_idxs, active.size() + 1);

  // Allocate space for seagull storage
  long num_sg = sg_idxs[active.size()];
  seagull* seagulls = newA(seagull, num_sg);

  // Store seagulls in parallel
  parallel_for(long i=0; i < active.size(); ++i) {
    uintE u_idx = active.vtx(i);
    const vertex u = U[u_idx];
    const uintE u_deg = u.getOutDegree();
    long sg_idx = sg_idxs[i];
    // Consider each neighbor v of active vertex u
    parallel_for(long j=0; j < u_deg; ++j) {
      const vertex v = V[u.getOutNeighbor(j)];
      const uintE v_deg = v.getOutDegree();
      long nbhd_idx = (nbhd_idxs[i])[j];
      long idx = 0;
      // Find neighbors (not equal to u) of v
      for (long k = 0; k < v_deg; ++k) {
        uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx != u_idx) {
          seagulls[sg_idx+nbhd_idx+idx] = cons(u_idx,u2_idx);
          ++idx;
        }
      }
    }
  }

  // Cleanup
  parallel_for(long i=0; i<active.size();++i) { free(nbhd_idxs[i]); }
  free(nbhd_idxs);
  free(sg_idxs);

  return make_pair(seagulls, num_sg);
}

//***************************************************************************************************
//***************************************************************************************************

/*
 *  Retrieves the number of butterflies remaining on each non-active vertex (where the number of 
 *  butterflies has changed), given an array of seagulls (constructed from some set of active vertices).
 *  Uses repeated sorting to obtain these counts.
 * 
 *  nu         : The number of vertices in U, which represents one bipartition of vertices, of which 
 *               the active vertices are taken to be a subset of
 *  seagulls   : Array of seagulls constructed from some set of active vertices. The active vertex that
 *               forms the endpoint of a seagull is assumed to be the first vertex in the unordered
 *               vertex pair, and the non-active vertex is assumed to be the second vertex.
 *  num_sgs    : The number of seagulls
 *  butterflies: An updated array containing the number of butterflies on each vertex in U (where
 *               the vertex is given by the index in butterflies)
 * 
 *  Returns: Array of butterfly counts to be changed with corresponding vertex, and the number of
 *           such changes
 */
pair<tuple<uintE, uintE>*, long> getSeagullFreqs(const long nu, UVertexPair* seagulls, long num_sgs, uintE* butterflies) {
  using X = tuple<uintE,uintE>;
  // Sort seagulls (considering both active + non-active endpoints), and retrieve frequency counts
  pair<uintE*, long> freq_pair = getFreqs(seagulls, num_sgs, UVertexPairCmp(nu), UVertexPairEq());
  long num_sg_freqs = freq_pair.second - 1;
  X* sg_freqs = newA(X, num_sg_freqs);
  // When retrieving frequency counts, store the frequency choose 2 with the non-active endpoint
  // This gives us the number of butterflies to be removed on the non-active endpoint
  parallel_for(long i=1; i < freq_pair.second; ++i) {
    uintE num = freq_pair.first[i] - freq_pair.first[i-1];
    uintE idx = seagulls[freq_pair.first[i-1]].v2;
    sg_freqs[i-1] = make_tuple(idx, num * (num-1) / 2;);
  }
  free(freq_pair.first);

  // Filter out any entries that have 0 butterflies to be removed (bucketing cannot handle these)
  num_sg_freqs = sequence::filter(sg_freqs, sg_freqs, num_sg_freqs, nonZeroF());
  
  // Now, collate all butterflies to be removed with the same non-active endpoint
  // Do this by sorting on the non-active endpoint, and summing the frequencies
  pair<uintE*, long> sg_freq_pair = getFreqs(sg_freqs, num_sg_freqs, uintETupleLt(), uintETupleEq());
  long num_sg_freqs_f = sg_freq_pair.second - 1;
  X* sg_freqs_f = newA(X, num_sg_freqs_f);
  parallel_for(long i=1; i < sg_freq_pair.second; ++i) {
    uintE num_freq = sg_freq_pair.first[i] - sg_freq_pair.first[i-1];
    // Reduce to sum the butterflies over the necessary range
    X sg_freq = sequence::reduce(&(sg_freqs[sg_freq_pair.first[i-1]]), num_freq, uintETupleAdd());
    uintE u_idx = get<0>(sg_freq);
    // Remove these butterflies from our array of butterfly counts
    butterflies[u_idx] -= get<1>(sg_freq);
    // Add the updated butterfly count to an update array, to be returned
    sg_freqs_f[i-1] = make_tuple(u_idx, butterflies[u_idx]);
  }
  free(sg_freq_pair.first);

  return make_pair(sg_freqs_f, num_sg_freqs_f);
}

/*
 *  Precisely getSeagullFreqs, but using histograms instead of repeated sortings.
 */
pair<tuple<uintE, uintE>*, long> getSeagullFreqsHist(const long nu, uintE* seagulls, long num_sgs, uintE* butterflies) {
  using X = tuple<uintE,uintE>;
  // TODO integrate sequence into histogram code (so we don't have to convert?)
  pbbsa::sequence<uintE> sgs_seq = pbbsa::sequence<uintE>(seagulls,num_sgs);

  // Place seagulls into a histogram to retrieve frequency counts (considering both active + non-active endpoints)
  tuple<size_t,X*> sgs_tuple = pbbsa::sparse_histogram<uintE,uintE>(sgs_seq,nu);
  X* sgs_freqs = get<1>(sgs_tuple);
  // Filter out any frequency count <= 1, since these won't contribute towards a butterfly
  size_t sgs_freqs_n = sequence::filter(sgs_freqs, sgs_freqs, get<0>(sgs_tuple), greaterOneF());
  parallel_for(long i=0; i < sgs_freqs_n; ++i) {
    uintE num = get<1>(sgs_freqs[i]);
    sgs_freqs[i] = make_tuple(get<0>(sgs_freqs[i]) % nu, num * (num-1) / 2);
  }

  // Now, we have to collate our seagulls again
  pbbsa::sequence<X> sgs_freqs_seq = pbbsa::sequence<X>(sgs_freqs,sgs_freqs_n);
  tuple<size_t,X*> sgs_freqs_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(sgs_freqs_seq,nu,chooseAdd<uintE,uintE>,chooseAddReduce<uintE,uintE>);
  X* sgs_freqs_f = get<1>(sgs_freqs_tuple);
  size_t num_sgs_freqs_f = get<0>(sgs_freqs_tuple);

  parallel_for(long i=0; i < num_sgs_freqs_f; ++i) {
    uintE u_idx = get<0>(sgs_freqs_f[i]);
    butterflies[u_idx] -= get<1>(sgs_freqs_f[i]);
    sgs_freqs_f[i] = make_tuple(u_idx, butterflies[u_idx]);
  }

  return make_pair(sgs_freqs_f, (long) num_sgs_freqs_f);
}

/*
 *  Assuming that the vertices in active are removed, computes the number of butterflies that are
 *  deleted by this removal from non-active vertices (where active and non-active vertices are in
 *  the same bipartition of the bipartite graph). Stores the number of removed butterflies in a 
 *  hash table seagulls_total, where the key is the non-active vertex. The vertices in active are 
 *  assumed to be a subset of U, and V represents the other bipartition of the bipartite graph.
 * 
 *  seagulls_total: Hash table to store the number of removed butterflies, where the key is the
 *                  non-active vertex
 *  active        : Set of active vertices (susbet of U)
 *  V             : One bipartition of vertices
 *  U             : The other bipartition of vertices
 *  nu            : The number of vertices in U
 * 
 *  Returns: None
 */
template<class vertex>
void getSeagullFreqsHash(sparseAdditiveSet<uintE>& seagulls_total, vertexSubset active, vertex* V, vertex* U, const long nu) {
  parallel_for (long i=0; i < active.size(); ++i) {
    // Set up for each active vertex
    uintE u_idx = active.vtx(i);
    const vertex u = U[u_idx];
    const uintE u_deg = u.getOutDegree();
    // Construct a (temporary) hash table to store seagulls with endpoint u (using key on non-active endpoint)
    float f = ((float) u_deg) / ((float) nu);
    sparseAdditiveSet<uintE> seagulls = sparseAdditiveSet<uintE>(nu,f,UINT_E_MAX);
    parallel_for (long j=0; j < u_deg; ++j ) {
        const vertex v = V[u.getOutNeighbor(j)];
        const uintE v_deg = v.getOutDegree();
        // Find all seagulls with center v
        parallel_for (long k=0; k < v_deg; ++k) {
          const uintE u2_idx = v.getOutNeighbor(k);
          if (u2_idx != u_idx) seagulls.insert(pair<uintE,uintE>(u2_idx, 1));
        }
    }
    _seq<pair<uintE,uintE>> sgs_freqs = seagulls.entries();
    // Accumulate in seagulls_total the removed butterflies, with active endpoint u (using key on non-active endpoint)
    parallel_for (long j=0; j < sgs_freqs.n; ++j) {
  	    pair<uintE,uintE> sgs_freq_pair = sgs_freqs.A[j];
  	    uintE num_butterflies = sgs_freq_pair.second;
        uintE u2_idx = sgs_freq_pair.first;
        // Removed butterflies are given by the number of seagulls choose 2; only insert non-zero
        // quantities of removed butterflies
        if (num_butterflies > 1) seagulls_total.insert(pair<uintE,uintE>(u2_idx, num_butterflies * (num_butterflies - 1)/2));
    }
    sgs_freqs.del();
    seagulls.del();
  }
}

//***************************************************************************************************
//***************************************************************************************************

/*
 *  Computes updated butterfly counts after active vertices are removed. The vertices in active are
 *  assumed to be a subset of U, and V represents the other bipartition of the bipartite graph.
 *  Uses sorting to obtain counts.
 * 
 *  active     : Set of active vertices (subset of U)
 *  butterflies: An updated array containing the number of butterflies on each vertex in U (where
 *               the vertex is given by the index in butterflies)
 *  V          : One bipartition of vertices
 *  U          : The other bipartition of vertices
 *  nu         : The number of vertices in U
 * 
 *  Returns: Array of butterfly counts to be changed with corresponding vertex, and the number of
 *           such changes
 */
template<class vertex>
pair<tuple<uintE, uintE>*, long> PeelSort(vertexSubset active, uintE* butterflies, vertex* V, vertex* U,const long nu) {
  // Retrieve all seagulls
  pair<UVertexPair*, long> sg_pair = getSeagulls<UVertexPair>(active, V, U, UVertexPairCons()); 
  // Compute updated butterfly counts
  return getSeagullFreqs(nu, sg_pair.first , sg_pair.second, butterflies);
}

/*
 *  Precisely PeelSort, but using histograms instead of repeated sortings.
 */
template<class vertex>
pair<tuple<uintE, uintE>*, long> PeelHist(vertexSubset active, uintE* butterflies, vertex* V, vertex* U,const long nu) {
  // Retrieve all seagulls
  pair<uintE*, long> sg_pair = getSeagulls<uintE>(active, V, U, UVertexPairIntCons(nu)); 
  // Compute updated butterfly counts
  return getSeagullFreqsHist(nu, sg_pair.first , sg_pair.second, butterflies);
}

/*
 *  Precisely PeelSort, but using hash tables instead of repeated sortings.
 */
template<class vertex>
pair<tuple<uintE, uintE>*, long> PeelHash(vertexSubset active, uintE* butterflies, vertex* V, vertex* U,const long nu) {
  // Compute number of butterflies to be removed
  sparseAdditiveSet<uintE> sgs_total = sparseAdditiveSet<uintE>(nu,(float) 1,UINT_E_MAX);
  getSeagullFreqsHash(sgs_total, active, V, U, nu);

  _seq<pair<uintE,uintE>> sgs_freqs = sgs_total.entries();
  long num_updates = sgs_freqs.n;
  using X = tuple<uintE,uintE>;
  X* update = newA(X, num_updates);

  // Compute the updated butterfly counts, given the number of butterflies to be removed
  parallel_for (long i=0; i < num_updates; ++i) {
  	pair<uintE,uintE> sgs_freq_pair = sgs_freqs.A[i];
    uintE u_idx = sgs_freq_pair.first;

    butterflies[u_idx] -= sgs_freq_pair.second;
    update[i] = make_tuple(u_idx, butterflies[u_idx]);
  }

  sgs_freqs.del();
  sgs_total.del();
  return make_pair(update, num_updates);
}

//***************************************************************************************************
//***************************************************************************************************

template <class vertex>
array_imap<uintE> Peel(bipartiteGraph<vertex>& GA, bool use_v, uintE* butterflies, long type=0, size_t num_buckets=128) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  vertex* V = use_v ? GA.V : GA.U;
  vertex* U = use_v ? GA.U : GA.V; // butterflies are on nu

  auto D = array_imap<uintE>(nu, [&] (size_t i) { return butterflies[i]; });

  auto b = make_buckets(nu, D, increasing, num_buckets);

  size_t finished = 0;
  while (finished != nu) {
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();

    pair<tuple<uintE, uintE>*, long> update_pair;
    if (type ==0) update_pair = PeelHash(active,butterflies,V,U,nu);
    else if(type==1) update_pair = PeelSort(active, butterflies, V,U,nu);
    else update_pair=PeelHist(active,butterflies,V,U,nu);

    vertexSubsetData<uintE> moved = vertexSubsetData<uintE>(nu, update_pair.second, update_pair.first);

    b.update_buckets(moved.get_fn_repr(), moved.size());
    
    moved.del(); active.del();
  }
  return D;
}

// Note: must be invoked with symmetricVertex
template <class vertex>
void Compute(hypergraph<vertex>& GA, commandLine P) {
  // Method type for counting + peeling
  long ty = P.getOptionLongValue("-t",0);
  long tp = P.getOptionLongValue("-tp",0);

  // Number of vertices if generating complete graph
  long nv = P.getOptionLongValue("-nv", 10);
  long nu = P.getOptionLongValue("-nu", 10);

  // 0 if using input file, 1 if using generated graph
  long gen = P.getOptionLongValue("-i",0);
    
	bipartiteGraph<symmetricVertex> G = (gen != 0) ? 
		bpGraphFromFile<symmetricVertex>(GA) : bpGraphComplete<symmetricVertex>(nv,nu);

  pair<bool,long> use_v_pair = cmpWedgeCounts(G);
  bool use_v = use_v_pair.first;
  long num_wedges = use_v_pair.second;

  timer t;
  t.start();
  uintE* butterflies = Count(G,use_v, num_wedges,ty);
	t.stop();

  if (ty==0) t.reportTotal("Sort:");
  else if (ty==1) t.reportTotal("SortCE:");
  else if (ty==2) t.reportTotal("Hash:");
  else if (ty==3) t.reportTotal("HashCE:");
  else if (ty==4) t.reportTotal("Hist:");
  else if (ty==5) t.reportTotal("HistNT:");
  else t.reportTotal("HistCE:"); 

  timer t2;
  t2.start();
  auto cores = Peel(G, use_v, butterflies, tp);
  t2.stop();
  if (tp ==0) t2.reportTotal("Hash Peel:");
  else if (tp==1) t2.reportTotal("Sort Peel:");
  else t2.reportTotal("Hist Peel:");

  long num_idxs = use_v ? G.nu : G.nv;
  uintE mc = 0;
  for (size_t i=0; i < num_idxs; i++) { mc = std::max(mc, cores[i]); }
  cout << "### Max core: " << mc << endl;

  free(butterflies);
	G.del();
}
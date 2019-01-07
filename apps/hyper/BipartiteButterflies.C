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
  // sg_idx is index of seagulls for each active vertex; nbhd_idx is index of seagulls for each neighbor of 
  // each active vertex, in a 2d nested array
  long* sg_idxs = newA(long, active.size() + 1);
  using T = long*;
  T* nbhd_idxs = newA(T, active.size());
  sg_idxs[active.size()] = 0;
  parallel_for(long i=0; i < active.size(); ++i) {
    uintE u_idx = active.vtx(i);
    const vertex u = U[u_idx];
    const uintE u_deg = u.getOutDegree();
    nbhd_idxs[i] = newA(long, u_deg + 1);
    (nbhd_idxs[i])[u_deg] = 0;
    parallel_for(long j=0; j < u_deg; ++j) {
      (nbhd_idxs[i])[j] = V[u.getOutNeighbor(j)].getOutDegree()-1;
    }
    sequence::plusScan(nbhd_idxs[i], nbhd_idxs[i], u_deg+1);
    sg_idxs[i] = (nbhd_idxs[i])[u_deg];
  }
  sequence::plusScan(sg_idxs, sg_idxs, active.size() + 1);
  long num_sg = sg_idxs[active.size()];
  seagull* seagulls = newA(seagull, num_sg);
  parallel_for(long i=0; i < active.size(); ++i) {
    uintE u_idx = active.vtx(i);
    const vertex u = U[u_idx];
    const uintE u_deg = u.getOutDegree();
    long sg_idx = sg_idxs[i];
    parallel_for(long j=0; j < u_deg; ++j) {
      const vertex v = V[u.getOutNeighbor(j)];
      const uintE v_deg = v.getOutDegree();
      long nbhd_idx = (nbhd_idxs[i])[j];
      long idx = 0;
      for (long k = 0; k < v_deg; ++k) {
        uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx != u_idx) {
          seagulls[sg_idx+nbhd_idx+idx] = cons(u_idx,u2_idx);
          ++idx;
        }
      }
    }
  }
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
  // Retrieve frequency counts on sorted seagulls
  pair<uintE*, long> freq_pair = getFreqs(seagulls, num_sgs, UVertexPairCmp(nu), UVertexPairEq());
  // This will hold our seagull frequencies choose 2, in (seagull end, frequency choose 2) tuples
  long num_sg_freqs = freq_pair.second - 1;
  X* sg_freqs = newA(X, num_sg_freqs);
  parallel_for(long i=1; i < freq_pair.second; ++i) {
    uintE num = freq_pair.first[i] - freq_pair.first[i-1];
    num = num * (num-1) / 2;
    uintE idx = seagulls[freq_pair.first[i-1]].v2;
    sg_freqs[i-1] = make_tuple(idx, num);
  }
  free(freq_pair.first);

  num_sg_freqs = sequence::filter(sg_freqs, sg_freqs, num_sg_freqs, nonZeroF());
  
  // Now, we have to collate our seagulls again
  pair<uintE*, long> sg_freq_pair = getFreqs(sg_freqs, num_sg_freqs, uintETupleLt(), uintETupleEq());
  long num_sg_freqs_f = sg_freq_pair.second - 1;
  X* sg_freqs_f = newA(X, num_sg_freqs_f);
  parallel_for(long i=1; i < sg_freq_pair.second; ++i) {
    uintE num_freq = sg_freq_pair.first[i] - sg_freq_pair.first[i-1];
    X sg_freq = sequence::reduce(&(sg_freqs[sg_freq_pair.first[i-1]]), num_freq, uintETupleAdd());
    uintE u_idx = get<0>(sg_freq);//get<0>(sg_freqs[sg_freq_pair.first[i-1]]); //
    butterflies[u_idx] -= get<1>(sg_freq);
    sg_freqs_f[i-1] = make_tuple(u_idx, butterflies[u_idx]);
  }
  free(sg_freq_pair.first);
  return make_pair(sg_freqs_f, num_sg_freqs_f);
}

pair<tuple<uintE, uintE>*, long> getSeagullFreqsHist(const long nu, uintE* seagulls, long num_sgs, uintE* butterflies) {
  using X = tuple<uintE,uintE>;
  pbbsa::sequence<uintE> sgs_seq = pbbsa::sequence<uintE>(seagulls,num_sgs);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.
  tuple<size_t,X*> sgs_tuple = pbbsa::sparse_histogram<uintE,uintE>(sgs_seq,nu);
  X* sgs_freqs = get<1>(sgs_tuple);
  size_t sgs_freqs_n = sequence::filter(sgs_freqs, sgs_freqs, get<0>(sgs_tuple), greaterOneF());
  //size_t sgs_freqs_n = get<0>(sgs_tuple);
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


template<class vertex>
void getSeagullFreqsHash(sparseAdditiveSet<uintE>& seagulls_total, vertexSubset active, vertex* V, vertex* U, const long nu) {
  parallel_for (long i=0; i < active.size(); ++i) {
      uintE u_idx = active.vtx(i);
      const vertex u = U[u_idx];
      const uintE u_deg = u.getOutDegree();
      float f = ((float) u_deg) / ((float) nu);
      sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(nu,f,UINT_E_MAX);
      parallel_for (long j=0; j < u_deg; ++j ) {
        const vertex v = V[u.getOutNeighbor(j)];
        const uintE v_deg = v.getOutDegree();
        parallel_for (long k=0; k < v_deg; ++k) {
          const uintE u2_idx = v.getOutNeighbor(k);
          if (u2_idx != u_idx) wedges.insert(pair<uintE,uintE>(u2_idx, 1));
        }
      }
      _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
      parallel_for (long j=0; j < wedge_freqs.n; ++j) {
  	    pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[j];
  	    uintE num_butterflies = wedge_freq_pair.second;
        uintE u2_idx = wedge_freq_pair.first;
        if (num_butterflies > 1) seagulls_total.insert(pair<uintE,uintE>(u2_idx, num_butterflies * (num_butterflies - 1)/2));
      }
      wedge_freqs.del();
      wedges.del();
    }
}

//***************************************************************************************************
//***************************************************************************************************

template<class vertex>
pair<tuple<uintE, uintE>*, long> PeelSort(vertexSubset active, uintE* butterflies, vertex* V, vertex* U,const long nu) {
  pair<UVertexPair*, long> sg_pair = getSeagulls<UVertexPair>(active, V, U, UVertexPairCons()); 
  return getSeagullFreqs(nu, sg_pair.first , sg_pair.second, butterflies);
}

template<class vertex>
pair<tuple<uintE, uintE>*, long> PeelHist(vertexSubset active, uintE* butterflies, vertex* V, vertex* U,const long nu) {
  pair<uintE*, long> sg_pair = getSeagulls<uintE>(active, V, U, UVertexPairIntCons(nu)); 
  return getSeagullFreqsHist(nu, sg_pair.first , sg_pair.second, butterflies);
}

template<class vertex>
pair<tuple<uintE, uintE>*, long> PeelHash(vertexSubset active, uintE* butterflies, vertex* V, vertex* U,const long nu) {
    sparseAdditiveSet<uintE> wedges_total = sparseAdditiveSet<uintE>(nu,(float) 1,UINT_E_MAX);
    getSeagullFreqsHash(wedges_total, active, V, U, nu);
    
    //go through wedges, and for each entry, we add to S a tuple with idx and val orig val - (n choose 2)
    _seq<pair<uintE,uintE>> wedge_freqs = wedges_total.entries();
    using X = tuple<uintE,uintE>;
    X* update = newA(X, wedge_freqs.n);

    parallel_for (long i=0; i < wedge_freqs.n; ++i) {
  	  pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
  	  uintE num_butterflies = wedge_freq_pair.second;
      uintE u = wedge_freq_pair.first;

      butterflies[u] -= num_butterflies;
      update[i] = make_tuple(u, butterflies[u]);
    }
    wedge_freqs.del();
    wedges_total.del();
    return make_pair(update, wedge_freqs.n);
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

template <class vertex>
void Compute(hypergraph<vertex>& GA, commandLine P) {
    std::string iFileConst = P.getOptionValue("-i", "");
    char* iFile = new char[iFileConst.length()+1];
    strcpy(iFile, iFileConst.c_str());
    long ty = P.getOptionLongValue("-t",0);
    long tp = P.getOptionLongValue("-tp",0);
  	long nv = P.getOptionLongValue("-nv", 10);
  	long nu = P.getOptionLongValue("-nu", 10);
	bipartiteGraph<symmetricVertex> G = (iFileConst.length() != 0) ? 
		bpGraphFromFile<symmetricVertex>(iFile) : bpGraphComplete<symmetricVertex>(nv,nu);
  delete [] iFile;

  uintE* butterflies;
  pair<bool,long> use_v_pair = cmpWedgeCounts(G);
  bool use_v = use_v_pair.first;
  long num_wedges = use_v_pair.second;

  timer t;
  t.start();
  butterflies=Count(G,use_v, num_wedges,ty);
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
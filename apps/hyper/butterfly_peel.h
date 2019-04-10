#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

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

#include "butterfly_utils.h"

/*
sparseAdditiveSet<uintE> getSeagullFreqsHash(sparseAdditiveSet<uintE>& seagulls, const long nu) {
  sparseAdditiveSet<uintE> update_hash = sparseAdditiveSet<uintE>(nu, (float) 1, UINT_E_MAX);
  _seq<pair<uintE,uintE>> sg_freqs = seagulls.entries();

  using T = pair<uintE,uintE>;

  //parallel_for (long j=0; j < sg_freqs.n; ++j) {
  granular_for ( j, 0, sg_freqs.n, (sg_freqs.n > 1000), {
    T sg_freq_pair = sg_freqs.A[j];
    uintE num = sg_freq_pair.second;
    if (num > 1) {update_hash.insert(T(sg_freq_pair.first % nu, num * (num - 1)/2));}
  });
  sg_freqs.del();
  return update_hash;
}

pair<uintE*, long> getUpdatesHash(sparseAdditiveSet<uintE>& update_hash, uintE* butterflies) {
  _seq<pair<uintE,uintE>> update_seq = update_hash.entries();
  long num_updates = update_seq.n;
  uintE* update = newA(uintE, num_updates);

  using T = pair<uintE,uintE>;

  //parallel_for (long i=0; i < num_updates; ++i) {
  granular_for(i, 0, num_updates, (num_updates>1000), {
    T update_pair = update_seq.A[i];
    uintE u_idx = update_pair.first;
    butterflies[u_idx] -= update_pair.second;
    update[i] = u_idx;
  });

  update_seq.del();
  return make_pair(update, num_updates);
}*/

template<class vertex>
pair<uintE*, long> PeelEHash(vertexSubset active, uintE* butterflies, vertex* V, vertex* U, const long nu, const long nv) {
  sparseAdditiveSet<uintE> update_hash = sparseAdditiveSet<uintE>(nv*nu+nv, (float) 1, UINT_E_MAX); //TODO find num updates lol
  parallel_for(long i=0; i < active.size(); ++i){
    // Set up for each active vertex
    uintE u_idx = active.vtx(i) % nu;
    uintE v_idx = active.vtx(i) / nu;
    //cout << active.vtx(i) << ", " << v_idx << ", " << u_idx << "\n";
    //fflush(stdout);
    vertex u = U[u_idx];
    vertex v = V[v_idx];
    uintE* v_neighbors = v.getOutNeighbors();

    parallel_for(long j=0; j < u.getOutDegree(); ++j) {
      uintE v_prime_idx = u.getOutNeighbor(j);
      if (v_prime_idx != v_idx) {
      vertex v_prime = V[v_prime_idx];
      uintE* v_prime_neighbors = v_prime.getOutNeighbors();
      tuple<long, uintE*> inter = intersect(v_neighbors, v_prime_neighbors, v.getOutDegree(), v_prime.getOutDegree());
      uintE num_butterflies = get<0>(inter) - 1;

      parallel_for(long k=0; k < get<0>(inter); ++k) {
        if (get<1>(inter)[k] == u_idx) {
        update_hash.insert(make_pair(get<1>(inter)[k] * nv + v_idx, num_butterflies));
        update_hash.insert(make_pair(get<1>(inter)[k] * nv + v_prime_idx, num_butterflies));
        }
      }
      free(get<1>(inter));
      }
    }
  }
  _seq<pair<uintE,uintE>> update_seq = update_hash.entries();
  long num_updates = update_seq.n;
  uintE* update = newA(uintE, num_updates);

  using T = pair<uintE,uintE>;

  parallel_for (long i=0; i < num_updates; ++i) {
    T update_pair = update_seq.A[i];
    uintE u_idx = update_pair.first;
    butterflies[u_idx] -= update_pair.second;
    update[i] = u_idx;
  }

  update_seq.del();
  update_hash.del();

  return make_pair(update, num_updates);
}

template <class vertex>
array_imap<uintE> PeelE(bipartiteGraph<vertex>& GA, bool use_v, uintE* butterflies, long type=0, size_t num_buckets=128) {
  // Butterflies are assumed to be stored on U
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  vertex* V = use_v ? GA.V : GA.U;
  vertex* U = use_v ? GA.U : GA.V;
  
  using X = tuple<uintE,uintE>;

  auto D = array_imap<uintE>(nu*(nv-1)+nu-1, [&] (size_t i) { return butterflies[i]; });

  auto b = make_buckets(nu*(nv-1)+nu-1, D, increasing, num_buckets);

  size_t finished = 0;
  while (finished != nu*(nv-1)+nu-1) {
    /*auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();

    bool is_seq = (active.size() < 1000);

    // Obtain butterfly updates based on desired method
    pair<uintE*, long> update_pair;
    if (type == 0) {
      update_pair = PeelEHash(active, butterflies, V, U, nu, nv);
    }
    //else if(type == 1){
    //  if (is_seq) update_pair = PeelSort_seq(active, butterflies, V, U, nu);
    //  else update_pair = PeelSort(active, butterflies, V, U, nu);
    //}
    //else {
    //  if (is_seq) update_pair = PeelHist_seq(active, butterflies, V, U, nu);
    //  else update_pair = PeelHist(active, butterflies, V, U, nu);
    //}

    pair<tuple<uintE,uintE>*,long> bucket_pair;
    if (is_seq) bucket_pair = updateBuckets_seq(update_pair.first, update_pair.second, butterflies, D, b, k);
    else bucket_pair = updateBuckets(update_pair.first, update_pair.second, butterflies, D, b, k);

    free(update_pair.first);

    vertexSubsetData<uintE> moved = vertexSubsetData<uintE>(nu, bucket_pair.second, bucket_pair.first);
    b.update_buckets(moved.get_fn_repr(), moved.size());

    moved.del(); active.del();*/
  }
  return D;
}
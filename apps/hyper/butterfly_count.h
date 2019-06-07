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

void countButterfliesSort(long* butterflies, UVertexPair* wedges, long* wedge_freqs_f, long num_wedge_freqs_f) {
  // Then, retrieve the frequency counts for each distinct key by taking the difference between
  // these indices
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // Store butterfly counts with each vertex in U
  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for (long i = 0; i < num_wedge_freqs_f-1; ++i) {
    long num_butterflies = wedge_freqs_f[i+1] - wedge_freqs_f[i];
    long wedge_idx = wedge_freqs_f[i];

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v1],num_butterflies); 
    writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v2],num_butterflies); 
  }
}

// This is the original compute function, without the more cache-efficient sorting method
intT CountSort(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  pair<long, intT> wedges_pair  = getWedges<UVertexPair>(cs.wedges_seq_uvp, GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges,
    wedge_idxs);

  UVertexPair* wedges = cs.wedges_seq_uvp.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UVertexPairCmp(), UVertexPairEq(), LONG_MAX, nonMaxLongF());

  countButterfliesSort(butterflies, wedges, freq_pair.first, freq_pair.second);

  free(freq_pair.first);
  return wedges_pair.second;
}

intT CountSortCE(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UVertexPair>(cs.wedges_seq_uvp, GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges, wedge_idxs);
  UVertexPair* wedges = cs.wedges_seq_uvp.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate by v1, v2
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UVertexPairCmp(), UVertexPairEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  using X = tuple<uintE, long>;
  if (cs.butterflies_seq_intt.n < 2*freq_pair.second-2) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*freq_pair.second-2);
    cs.butterflies_seq_intt.n = 2*freq_pair.second-2;
  }

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num = freq_arr[i+1] - freq_arr[i];
    num = num * (num-1) / 2;
    //parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
    long j = freq_arr[i];
    cs.butterflies_seq_intt.A[2*i] = make_tuple(wedges[j].v1, num);
    cs.butterflies_seq_intt.A[2*i+1] = make_tuple(wedges[j].v2, num);
    //}
  }

  free(freq_arr);

  // now, we need to collate by our indices
  pair<long*, long> b_freq_pair = getFreqs<long>(cs.butterflies_seq_intt.A, 2*freq_pair.second-2, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;
  const intT eltsPerCacheLine = 64/sizeof(long);

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    long num_freq = b_freq_arr[i] - b_freq_arr[i-1];
    // Reduce to sum the butterflies over the necessary range
    X reduce = sequence::reduce(&(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]), num_freq, tupleAdd<uintE,long>());
    // These are our butterfly counts
    butterflies[get<0>(reduce)*eltsPerCacheLine] += get<1>(reduce);
  }

  free(b_freq_arr);

  return wedges_pair.second;
}

intT CountSortCE(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  const intT eltsPerCacheLine = 64/sizeof(long);
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx, num_wedges,
    wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_f = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // We need to first collate
  auto freq_pair = getFreqs<long>(wedges, num_wedges_f, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  auto freq_arr = freq_pair.first;

  using X = tuple<uintE,long>;
  if (cs.butterflies_seq_intt.n < num_wedges_f) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, num_wedges_f);
    cs.butterflies_seq_intt.n = num_wedges_f;
  }

  // store these counts in another array so we can store in CE manner
  parallel_for(long i=0; i < freq_pair.second - 1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];
    if ((wedges[wedge_idx].v2 & 0b1) && num_butterflies > 1) {
    num_butterflies = num_butterflies * (num_butterflies - 1)/2;
    cs.butterflies_seq_intt.A[freq_arr[i]] = make_tuple(wedges[wedge_idx].v1, num_butterflies);
    cs.butterflies_seq_intt.A[freq_arr[i]+1] = make_tuple(wedges[wedge_idx].v2 >> 1, num_butterflies);
    parallel_for(long j=freq_arr[i]+2; j < freq_arr[i+1]; ++j) {cs.butterflies_seq_intt.A[j] = make_tuple(UINT_E_MAX, 0);}
    } else if (!(wedges[wedge_idx].v2 & 0b1) && num_butterflies > 1){
    parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
      //writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].u, num_butterflies - 1);
      cs.butterflies_seq_intt.A[j] = make_tuple(wedges[j].u, num_butterflies - 1);
    }
    }
    else {
      parallel_for(long j=freq_arr[i]; j < freq_arr[i+1]; ++j) {
        cs.butterflies_seq_intt.A[j] = make_tuple(UINT_E_MAX, 0);
      }
    }
  }

  free(freq_arr);

  // now, we need to collate by our indices
  auto b_freq_pair = getFreqs<long>(cs.butterflies_seq_intt.A, num_wedges_f, tupleLt<uintE,long>(), tupleEq<uintE,long>(), LONG_MAX, nonMaxLongF());
  auto b_freq_arr = b_freq_pair.first;

  parallel_for(long i=1; i < b_freq_pair.second; ++i) {
    if (get<0>(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]) != UINT_E_MAX) {
      long num_freq = b_freq_arr[i] - b_freq_arr[i-1];
      // Reduce to sum the butterflies over the necessary range
      X reduce = sequence::reduce(&(cs.butterflies_seq_intt.A[b_freq_arr[i-1]]), num_freq, tupleAdd<uintE,long>());
      // These are our butterfly counts
      butterflies[eltsPerCacheLine*get<0>(reduce)] += get<1>(reduce);
    }
  }

  free(b_freq_arr);
  return wedges_pair.second;
}

intT CountSort(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx, num_wedges,
    wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());

  const intT eltsPerCacheLine = 64/sizeof(long);
  parallel_for (long i = 0; i < freq_pair.second-1; ++i) {
    long num_butterflies = freq_pair.first[i+1] - freq_pair.first[i];
    long wedge_idx = freq_pair.first[i];

    if ((wedges[wedge_idx].v2 & 0b1) && (num_butterflies > 1)) {
      num_butterflies = (num_butterflies * (num_butterflies - 1))/2;
      writeAdd(&butterflies[eltsPerCacheLine*wedges[wedge_idx].v1],num_butterflies); 
      writeAdd(&butterflies[eltsPerCacheLine*(wedges[wedge_idx].v2 >> 1)],num_butterflies);
    }
    else if (!(wedges[wedge_idx].v2 & 0b1) && (num_butterflies > 1)) {
      parallel_for(long j=freq_pair.first[i]; j < freq_pair.first[i+1]; ++j) {
        writeAdd(&butterflies[eltsPerCacheLine*wedges[j].u], num_butterflies - 1);
      }
    }
  }

  free(freq_pair.first);
  return wedges_pair.second;
}

//********************************************************************************************
//********************************************************************************************

//cs, g, num_wedges, use_v, nv, nu, butterflies, max_wedges, wedge_idxs, ranks
intT CountHash(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
  intT curr_idx=0) {
  const intT eltsPerCacheLine = 64/sizeof(long);

  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);
  intT num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    if (num_butterflies > 1 && (wedge_freq_pair.first & 0b1)) {
      uintE u2 = ((wedge_freq_pair.first >> 1) % (GA.n));
      uintE u = ((wedge_freq_pair.first >> 1) / (GA.n));
      num_butterflies = (num_butterflies * (num_butterflies - 1))/2;
      writeAdd(&butterflies[eltsPerCacheLine*u], num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*u2], num_butterflies);
    }
    // if slow, could try to do an intersect here to find v
  }

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for (intT j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
      if (v > i && (GA.edges[u_offset+j] & 0b1)) {
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2 = GA.edges[v_offset+k] >> 1;
        if (u2 > i) {
          long to_find = (((long) i) *GA.n + (long) u2) << 1;
          long num_butterflies = cs.wedges_hash.find(to_find).second;
          if (num_butterflies > 1) writeAdd(&butterflies[eltsPerCacheLine*v], num_butterflies - 1);
        }
        else break;
      }
      }
    }
  }

  return next_idx;
}

intT CountHashCE(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
  intT curr_idx=0) {
  using T = pair<uintE,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);
  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);
  intT num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    if (num_butterflies > 1 && (wedge_freq_pair.first & 0b1)) {
      uintE u2 = (wedge_freq_pair.first >> 1) % (GA.n);
      uintE u = (wedge_freq_pair.first >> 1) / (GA.n);
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      cs.butterflies_hash.insert(T(u, num_butterflies));
      cs.butterflies_hash.insert(T(u2, num_butterflies));
    } 
  }
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for (intT j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
      if (v > i && (GA.edges[u_offset+j] & 0b1)) {
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2 = GA.edges[v_offset+k] >> 1;
        if (u2 > i) {
          long to_find = ((long)i*GA.n + u2) << 1;
          long num_butterflies = cs.wedges_hash.find(to_find).second;
          if (num_butterflies > 1) cs.butterflies_hash.insert(T(v, num_butterflies-1));
        }
        else break;
      }
      }
    }
  }

  num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.butterflies_seq_intp);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterfly_pair = cs.butterflies_seq_intp.A[i];
    butterflies[eltsPerCacheLine*butterfly_pair.first] += butterfly_pair.second;
  }

  return next_idx;
}

// TODO modularize this stuff better
intT CountHash(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
  long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  
  intT num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp); //TODO should be long num?

  // Retrieve count on each key; that number choose 2 is the number of butterflies  
  parallel_for (long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;

    if (num_butterflies > 1){
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      writeAdd(&butterflies[eltsPerCacheLine*v1],num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*v2],num_butterflies);
    }
  }

  return next_idx;
}

intT CountHashCE(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges,
  long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  //getWedgesFromHashTimer.start();
  intT num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);
  //getWedgesFromHashTimer.stop();
  //numButterfliesHashInsertTimer.start();
  using T = pair<uintE, long>;
  // Retrieve count on each key; that number choose 2 is the number of butterflies
  parallel_for (long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    uintE v1 = wedge_freq_pair.first / nu;
    uintE v2 = wedge_freq_pair.first % nu;

    if (num_butterflies > 1) {
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      cs.butterflies_hash.insert(T(v1, num_butterflies));
      cs.butterflies_hash.insert(T(v2, num_butterflies));
    }
  }
  num_wedges_seq = cs.butterflies_hash.entries_no_init(cs.butterflies_seq_intp);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto butterfly_pair = cs.butterflies_seq_intp.A[i];
    butterflies[eltsPerCacheLine*butterfly_pair.first] += butterfly_pair.second;
  }
  //numButterfliesHashInsertTimer.stop();  
  return next_idx;
}

//********************************************************************************************
//********************************************************************************************

intT CountHist(CountSpace& cs, graphCSR& GA, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, 
  intT curr_idx=0) {
  using X = tuple<long,uintE>;
  const intT eltsPerCacheLine = 64/sizeof(long);

  pair<long, intT> wedges_pair  = getWedges<long>(cs.wedges_seq_int, GA, UWedgeIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);
  intT next_idx = wedges_pair.second;

  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(cs.wedges_seq_int.A, wedges_pair.first);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  cs.wedges_hash.resize(wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    cs.wedges_hash.insert(make_pair(get<0>(wedge_freqs[i]), get<1>(wedge_freqs[i])));
  }

  intT num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num_butterflies = wedge_freq_pair.second;
    if (num_butterflies > 1 && (wedge_freq_pair.first & 0b1)) {
      uintE u2 = ((wedge_freq_pair.first >> 1) % (GA.n));
      uintE u = ((wedge_freq_pair.first >> 1) / (GA.n));
      num_butterflies = (num_butterflies * (num_butterflies - 1))/2;
      writeAdd(&butterflies[eltsPerCacheLine*u], num_butterflies);
      writeAdd(&butterflies[eltsPerCacheLine*u2], num_butterflies);
    }
    // if slow, could try to do an intersect here to find v
  }

  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    parallel_for (intT j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
      if (v > i && (GA.edges[u_offset+j] & 0b1)) {
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2 = GA.edges[v_offset+k] >> 1;
        if (u2 > i) {
          long to_find = (((long) i) *GA.n + (long) u2) << 1;
          long num_butterflies = cs.wedges_hash.find(to_find).second;
          if (num_butterflies > 1) writeAdd(&butterflies[eltsPerCacheLine*v], num_butterflies - 1);
        }
        else break;
      }
      }
    }
  }

  return next_idx;
}

intT CountHist(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  const intT eltsPerCacheLine = 64/sizeof(long);

  using T = tuple<long,uintE>;

  pair<long, intT> wedges_list_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_curr = wedges_list_pair.first;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_curr);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.

  tuple<size_t,T*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  T* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    auto wedge_freq_pair = wedge_freqs[i];
    long num_butterflies = get<1>(wedge_freq_pair);
    auto wedge_freq_pair_first = get<0>(wedge_freq_pair);

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
intT CountHistCE(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long* butterflies, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  const long nv = use_v ? GA.nv : GA.nu;
  const intT eltsPerCacheLine = 64/sizeof(long);
 
  using T = tuple<long, uintE>;
  using X = tuple<uintE, long>;

  pair<long, intT> wedges_list_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_curr = wedges_list_pair.first;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_curr);

  tuple<size_t,T*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nv*nu + nu, cs.tmp, cs.out);
  T* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  if (cs.butterflies_seq_intt.n < 2*wedge_freqs_n) {
    free(cs.butterflies_seq_intt.A);
    cs.butterflies_seq_intt.A = newA(X, 2*wedge_freqs_n);
    cs.butterflies_seq_intt.n = 2*wedge_freqs_n;
  }

  parallel_for(long i = 0; i < wedge_freqs_n; i++) {
    auto wedge_freq_pair = wedge_freqs[i];
    long num = get<1>(wedge_freq_pair);
    num = (num * (num-1)) / 2;
    auto wedge_num = get<0>(wedge_freq_pair);
    cs.butterflies_seq_intt.A[2*i] = make_tuple(wedge_num % nu, num);
    cs.butterflies_seq_intt.A[2*i + 1] = make_tuple(wedge_num / nu, num);
  }

  pbbsa::sequence<X> wedge_freqs_i_seq = pbbsa::sequence<X>(cs.butterflies_seq_intt.A,2*wedge_freqs_n);
  tuple<size_t, X*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,long>(wedge_freqs_i_seq,nu,getAdd<uintE,long>, getAddReduce<uintE,long>, cs.tmp_uint, cs.out_uint);
  X* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);

  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[eltsPerCacheLine*get<0>(butterflies_l[i])] += get<1>(butterflies_l[i]);
  }

  return wedges_list_pair.second;
}

//********************************************************************************************
//********************************************************************************************


void CountOrigCompactParallel(bipartiteCSR& GA, bool use_v, long* butterflies) {
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
      long shift = nu*(i-step*stepSize);
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
        writeAdd(&butterflies[i*eltsPerCacheLine],  (long)((long) wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
        writeAdd(&butterflies[u2_idx*eltsPerCacheLine], (long)((long) wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
        wedges[u2_idx+shift] = 0;
      }
    }
  }
  t2.reportTotal("main loop");

  free(wedges);
  free(used);
}

long* CountWorkEfficientParallel(graphCSR& GA, long* butterflies) {
  timer t1,t2;
  t1.start();

  long stepSize = getWorkers() * 7; //15 tunable parameter
  uintE* wedges = newA(uintE, GA.n*stepSize);
  uintE* used = newA(uintE, GA.n*stepSize);

  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);

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
	    if (wedges[shift+u2_idx] > 1) writeAdd(&butterflies[eltsPerCacheLine*v], (long)(wedges[shift+u2_idx]-1));
	  }
	  else break;
	}
      }

      for(long j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j] >> 1;
        if(used[shift+j] & 0b1) {
        writeAdd(&butterflies[i*eltsPerCacheLine],  (long)((long)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
        writeAdd(&butterflies[u2_idx*eltsPerCacheLine], (long)((long)wedges[shift+u2_idx]*(wedges[shift+u2_idx]-1) / 2));
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

long* CountRank(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long type,
  uintE* ranks, uintE* rankV, uintE* rankU) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  timer t_rank;
  t_rank.start();
  auto g = rankGraph(GA, use_v, ranks, rankV, rankU);
  free(ranks);
  t_rank.reportTotal("graph rank");

  timer t_time;
  t_time.start();

  const intT eltsPerCacheLine = 64/sizeof(long);
  long* rank_butterflies = newA(long,eltsPerCacheLine*g.n);
  granular_for(i,0,g.n,g.n > 10000, { rank_butterflies[eltsPerCacheLine*i] = 0; });

  long* wedge_idxs = (type == 11) ? nullptr : countWedgesScan(g);
  CountSpace cs = CountSpace(type, g.n, true);

  if (type == 11) CountWorkEfficientParallel(g, rank_butterflies);
  else {
  if (max_wedges >= num_wedges) {
    if (type == 0) CountSort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
    else if (type == 1) CountSortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
    else if (type == 2) CountHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
    else if (type == 3) CountHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
    else if (type == 4) CountHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs);
  }
  else {
  intT curr_idx = 0;
  while(curr_idx < g.n) {
    if (type ==0) curr_idx = CountSort(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type ==1) curr_idx = CountSortCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type ==2) curr_idx = CountHash(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type == 3) curr_idx = CountHashCE(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type == 4) curr_idx = CountHist(cs, g, num_wedges, rank_butterflies, max_wedges, wedge_idxs, curr_idx);
    cs.clear();
  }
  }
  }
  g.del();
  cs.del();
  if (type != 11) free(wedge_idxs);

  //uintE* rank_butterflies2 = newA(uintE,eltsPerCacheLine*g.n);
  //granular_for(i,0,g.n,g.n > 10000, { rank_butterflies2[eltsPerCacheLine*i] = 0; });
  //CountWorkEfficientParallel(g, rank_butterflies2);
  //for(long i=0; i < g.n; ++i) {assert(rank_butterflies2[eltsPerCacheLine*i] == rank_butterflies[eltsPerCacheLine*i]);}

  t_time.reportTotal("counting");

  timer t_convert;
  t_convert.start();

  long* butterflies = newA(long, eltsPerCacheLine*nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[eltsPerCacheLine * i] = 0; });

  auto rank_converter = use_v ? rankU : rankV;
  granular_for(i,0,nu,nu > 10000, { 
    butterflies[eltsPerCacheLine*i] = rank_butterflies[eltsPerCacheLine*rank_converter[i]];
  });
  free(rank_butterflies);
  free(rankU); free(rankV);

  t_convert.reportTotal("convert");

  return butterflies;
}

long* Count(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long type=0, long tw=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  // tw 0 is use side, tw 1 is use co core ranks, tw 2 is use approx co core ranks, tw 3 is use the better b/w approx co core and side (TODO put in deg for tw 3)
if (tw != 0) {
  timer t_rank;
  t_rank.start();
  //auto rank_tup = getDegRanks(GA);
  // auto rank_tup = getCoreRanks(GA);
  auto rank_tup = tw == 1 ? getCoCoreRanks(GA) : getApproxCoCoreRanks(GA);

    //long num_rwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup), get<1>(rank_tup)));
    //num_rwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup)));

    //long num_cwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup2), get<1>(rank_tup2)));
    //num_cwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup2), get<2>(rank_tup2)));

  long num_ccwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
    rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup), get<1>(rank_tup)));
  num_ccwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
    rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup)));
  
  t_rank.reportTotal("ranking");

  //cout << "Rank wedges: " << num_rwedges << "\n";
  cout << "Side wedges: " << num_wedges << "\n";
  //cout << "Core wedges: " << num_cwedges << "\n";
  cout << "Co Core wedges: " << num_ccwedges << "\n";

  if (num_ccwedges < num_wedges + 1000 || tw == 1 || tw == 2) return CountRank(GA, use_v, num_ccwedges, max_wedges, type, get<0>(rank_tup), get<1>(rank_tup), get<2>(rank_tup));
  free(get<0>(rank_tup)); free(get<1>(rank_tup)); free(get<2>(rank_tup));
}

  const intT eltsPerCacheLine = 64/sizeof(long);
  long* butterflies = newA(long, eltsPerCacheLine*nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[eltsPerCacheLine * i] = 0; });

  long* wedge_idxs = (type == 7 || type == 11) ? nullptr : countWedgesScan(GA, use_v, true);
  CountSpace cs = CountSpace(type, nu, false);

  if (type == 7) CountOrigCompactParallel(GA, use_v, butterflies);
  else if (type == 11) CountOrigCompactParallel(GA, use_v, butterflies);
  else {

  if (max_wedges >= num_wedges) {
    if (type == 0) CountSort(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else if (type == 1) CountSortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else if (type == 2) CountHash(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else if (type == 3) CountHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else if (type == 4) CountHist(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
    else if (type == 6) CountHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs);
  }
  else {
  intT curr_idx = 0;
  while(curr_idx < nu) {
    if (type == 0) CountSort(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type == 1) CountSortCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type == 2) CountHash(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type == 3) CountHashCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type == 4) CountHist(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    else if (type == 6) CountHistCE(cs, GA, use_v, num_wedges, butterflies, max_wedges, wedge_idxs, curr_idx);
    cs.clear();
  }
  }
  }

  if (type != 7 && type != 11) free(wedge_idxs);
  cs.del();
  return butterflies;

}

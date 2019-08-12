#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#define HYPER 1
#define MCX16 1
#define VERBOSE 1

#include "hypergraphIO.h"
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

#include "butterfly_count.h"
#include "butterfly_ecount.h"
#include "butterfly_peel.h"
#include "butterfly_epeel.h"
#include "butterfly_count_total.h"

#include <vector>

using namespace std;

//JS: to do: integrate into framework so that ranked versions work as well. also need one for edge counting.
void CountOrigCompactSerial(bipartiteCSR& GA, bool use_v) {
  timer t1,t2;
  t1.start();
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long results = 0;
  uintE* wedges = newA(uintE, nu);
  uintE* used = newA(uintE, nu);

  for(long i=0; i < nu; ++i) { wedges[i] = 0; }

  long* butterflies = newA(long,nu);
  for(long i=0; i < nu; ++i) { butterflies[i] = 0; }

  t1.reportTotal("preprocess");
  t2.start();

  for(intT i=0; i < nu; ++i){
    intT used_idx = 0;
    intT u_offset  = offsetsU[i];
    intT u_deg = offsetsU[i+1]-u_offset;
    for (intT j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1]-offsetsV[v];
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2_idx = edgesV[v_offset+k];
        if (u2_idx < i) {
          butterflies[i] += wedges[u2_idx];
          butterflies[u2_idx] += wedges[u2_idx];
          wedges[u2_idx]++;
          if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
        }
        else break;
      }
    }
    for(intT j=0; j < used_idx; ++j) { wedges[used[j]] = 0; }
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);

  for(long i=0;i<nu;i++) results += butterflies[i];
  free(butterflies);
  cout << "num: " << results/2 << "\n";
}

void CountOrigCompactSerialTotal(bipartiteCSR& GA, bool use_v) {
  timer t1,t2;
  t1.start();

  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long results = 0;
  uintE* wedges = newA(uintE, nu);
  uintE* used = newA(uintE, nu);

  for(long i=0; i < nu; ++i) { wedges[i] = 0; }
  long nb = 0;

  t1.reportTotal("preprocess");
  t2.start();

  for(intT i=0; i < nu; ++i){
    intT used_idx = 0;
    intT u_offset  = offsetsU[i];
    intT u_deg = offsetsU[i+1]-u_offset;
    for (intT j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1]-offsetsV[v];
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2_idx = edgesV[v_offset+k];
        if (u2_idx < i) {
          nb += wedges[u2_idx];
          wedges[u2_idx]++;
          if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
        }
        else break;
      }
    }
    for(intT j=0; j < used_idx; ++j) { wedges[used[j]] = 0; }
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);

  cout << "num: " << nb << "\n";
}

string CountTypeToStr(CountType ty) {
  switch(ty) {
    case ASORT: return "Sort"; break;
    case SORT: return "SortCE"; break;
    case AHASH: return "Hash"; break;
    case HASH: return "HashCE"; break;
    case AHIST: return "Hist"; break;
    case HIST: return "HistCE"; break;
    case BATCHS: return "Par"; break;
    case SERIAL: return "Serial"; break;
    case BATCHWA: return "WedgePar"; break;
    default: break;
  }
  return "";
}

void Compute(bipartiteCSR& GA, commandLine P) {
  const size_t eltsPerCacheLine = 64/sizeof(long);

  // Type for counting
  long count_type_long = P.getOptionLongValue("-t", LONG_MAX);
  string count_type_str = P.getOptionValue("-countType", "");
  CountType ty;
  if (count_type_long == LONG_MAX) {
    if (count_type_str == "ASORT") ty = ASORT;
    else if (count_type_str == "SORT") ty = SORT;
    else if (count_type_str == "AHASH") ty = AHASH;
    else if (count_type_str == "HASH") ty = HASH;
    else if (count_type_str == "AHIST") ty = AHIST;
    else if (count_type_str == "HIST") ty = HIST;
    else if (count_type_str == "BATCHS") ty = BATCHS;
    else if (count_type_str == "BATCHWA") ty = BATCHWA;
    else ty = SERIAL;
  }
  else {
    switch(count_type_long) {
      case 0: ty = ASORT; break;
      case 1: ty = SORT; break;
      case 2: ty = AHASH; break;
      case 3: ty = HASH; break;
      case 4: ty = AHIST; break;
      case 6: ty = HIST; break;
      case 5: case 7: case 11: ty = BATCHS; break;
      case 9: case 12: ty = SERIAL; break;
      case 8: ty = BATCHWA; break;
      default: break;
    }
  }

  // Type for ranking
  long rank_type_long = P.getOptionLongValue("-w", LONG_MAX);
  string rank_type_str = P.getOptionValue("-rankType", "");
  RankType tw;
  if (rank_type_long == LONG_MAX) {
    if (rank_type_str == "SIDE") tw = SIDE;
    else if (count_type_str == "COCORE") tw = COCORE;
    else if (count_type_str == "ACOCORE") tw = ACOCORE;
    else if (count_type_str == "DEG") tw = DEG;
    else tw = ADEG;
  }
  else {
    switch(rank_type_long) {
      case 0: tw = SIDE; break;
      case 1: tw = COCORE; break;
      case 2: tw = ACOCORE; break;
      case 3: tw = DEG; break;
      case 4: tw = ADEG; break;
      default: break;
    }
  }

  long tp = P.getOptionLongValue("-tp",0);
  long te = P.getOptionLongValue("-e",0);
  bool nopeel = P.getOptionValue("-nopeel");
  bool total = P.getOptionValue("-total");
  
  // # of max wedges
  long max_wedges = P.getOptionLongValue("-m",2577500000);
  long max_array_size = P.getOptionLongValue("-a",23090996160);

  long denom = P.getOptionLongValue("-d",25);
  long sparse = P.getOptionLongValue("-s",0); 
  if (sparse > 0) {total=true;}
  if (total) {nopeel = true;}

  cout << "count " << ty << ", " << "peel " << tp << ", " << "edge " << te << ", " << "rank " << tw << ", " << "peel " << nopeel << "\n";

  timer t1;
  t1.start();
  tuple<bool,long> use_v_tuple = cmpWedgeCounts(GA);
  bool use_v = get<0>(use_v_tuple);
  long num_wedges = get<1>(use_v_tuple);
  t1.reportTotal("preprocess (wedge counts)");

  // if you only want total counts
  if (total) {
    if (tw == SIDE && (ty == SERIAL)) CountOrigCompactSerialTotal(GA,use_v);
    if (tw == SIDE && (ty == SERIAL)) return;

    timer t;
    t.start();
    long num_butterflies = CountTotal(GA, use_v, num_wedges, max_wedges, max_array_size, ty, tw);
    t.stop();

    t.reportTotal(CountTypeToStr(ty));

    if (sparse == 1) num_butterflies *= pow(denom, 3);
    else if (sparse == 2) num_butterflies *= pow(denom, 4);
    cout << "number of butterflies: " << num_butterflies << "\n";
    return;
  }

 //TODO seq code integrate w/count
  if (te == 0) {
    if (tw == SIDE && ty == SERIAL) CountOrigCompactSerial(GA,use_v);
    if (tw == SIDE && ty == SERIAL) return;

    timer t;
    t.start();
    long* butterflies = Count(GA, use_v, num_wedges, max_wedges, max_array_size, ty, tw);
    t.stop();

    t.reportTotal(CountTypeToStr(ty));

    long num_idxs = use_v ? GA.nu : GA.nv;

    auto butterflies_extract_f = [&] (const long i) -> const long {
      return butterflies[eltsPerCacheLine*i];
    };
    
    long total = sequence::reduce<long>((long)0,(long)num_idxs,addF<long>(),butterflies_extract_f);
    
    cout << "number of butterflies: " << total/2 << "\n";

    if(!nopeel) {
      timer t2;
      t2.start();
      auto cores = Peel(GA, use_v, butterflies, max_wedges, tp, max_array_size);
      t2.stop();
      if (tp ==0) t2.reportTotal("Hash Peel:");
      else if (tp==1) t2.reportTotal("Sort Peel:");
      else if (tp==2) t2.reportTotal("Hist Peel:");
      else if (tp == 3) t2.reportTotal("Par Peel:");
      else if (tp == 5) t2.reportTotal("WedgePar Peel:");
    }
    free(butterflies);
  }
  else {
  
    timer t3,t4;
    t4.start();
    auto eti = edgeToIdxs(GA, use_v);
    t4.reportTotal("edgeToIdxs");
    
    t3.start();
    long* ebutterflies = CountE(eti, GA, use_v, num_wedges, max_wedges, max_array_size, ty, tw);
    t3.stop();
  
    t3.reportTotal(CountTypeToStr(ty));

    auto butterflies_extract_f = [&] (const long i) -> const long {
      return ebutterflies[eltsPerCacheLine*i];
    };
    
    long total = sequence::reduce<long>((long)0,(long)GA.numEdges,addF<long>(),butterflies_extract_f);
    cout << "number of edge butterflies: " << total/4 << "\n";

    if(!nopeel) {
      timer t2;
      t2.start();
      auto ite = idxsToEdge(GA, use_v);
      auto cores = PeelE(eti, ite, GA, use_v, ebutterflies, max_wedges, tp);
      free(ite);	    
      t2.stop();
      if (tp ==0) t2.reportTotal("Hash E Peel:");
      else if (tp==1) t2.reportTotal("Sort E Peel:");
      else if (tp == 2) t2.reportTotal("Hist E Peel:");
      else if (tp == 3) t2.reportTotal("Par E Peel:");
      else if (tp == 4) t2.reportTotal("NoUpdatePar E Peel:");
      else if (tp == 5) t2.reportTotal("WedgePar E Peel:");
    }
    free(eti);

    free(ebutterflies);
  }
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," <inFile>");
  char* iFile = P.getArgument(0);
  long rounds = P.getOptionLongValue("-r",3);
  long te = P.getOptionLongValue("-e",0);
  bool nopeel = P.getOptionValue("-nopeel");
  long denom = P.getOptionLongValue("-d",25);
  long sparse = P.getOptionLongValue("-s",0); // 0 for not sparse, 1 for clr, 2 for edge

  if ((te == 0 || nopeel) && sparse == 0) {
    bipartiteCSR G = readBipartite(iFile);

    Compute(G,P);
    for(int r=0;r<rounds;r++) {
      Compute(G,P);
    }
    G.del();
  }
  else {
  bipartiteCSR G;
  if (sparse == 0) G = readBipartite(iFile);
  else {
    auto tmp = readBipartite(iFile);
    G = sparse == 1 ? clrSparseBipartite(tmp, denom, 0) : eSparseBipartite(tmp, denom, 0);
    tmp.del();
  }
  Compute(G,P);
  G.del();

  for(int r=0;r<rounds;r++) {
    if (sparse == 0) G = readBipartite(iFile);
    else{
      auto tmp = readBipartite(iFile);
      G = sparse == 1 ? clrSparseBipartite(tmp, denom, (r+1)*(tmp.nv+tmp.nu)) : eSparseBipartite(tmp, denom, (r+1)*(tmp.nv+tmp.nu));
      tmp.del();
    }
    Compute(G,P);
    G.del();
  }
  }
}


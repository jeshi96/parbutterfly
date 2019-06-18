#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#define HYPER 1
//#define LONG 1
//#define EDGELONG 1
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
//#include "../../lib/gbbs-histogram.h"

#include "butterfly_count.h"
#include "butterfly_ecount.h"
#include "butterfly_peel.h"
#include "butterfly_epeel.h"
#include "butterfly_count_total.h"

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }

#include <vector>

using namespace std;

//JS: to do: integrate into framework so that ranked versions work as well. also need one for edge counting.
void CountOrigCompactSerial(bipartiteCSR& GA, bool use_v) {
  timer t1,t2;
  t1.start();
  //cout << GA.nv << " " << GA.nu << " " << GA.numEdges << endl;
  cout << "Original Serial (make sure running with CILK_NWORKERS=1)" << endl;  
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
          //results += wedges[u2_idx];
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

// Note: must be invoked with symmetricVertex
void Compute(bipartiteCSR& GA, commandLine P) {
  // Method type for counting + peeling
  long ty = P.getOptionLongValue("-t",0);
  long tp = P.getOptionLongValue("-tp",0);
  long te = P.getOptionLongValue("-e",0);
  long tw = P.getOptionLongValue("-w",0);
  bool nopeel = P.getOptionValue("-nopeel");
  bool total = P.getOptionValue("-total");
  
  // # of max wedges
  long max_wedges = P.getOptionLongValue("-m",2577500000);
  long max_array_size = P.getOptionLongValue("-a",23090996160);

  cout << "count " << ty << ", " << "peel " << tp << ", " << "edge " << te << ", " << "rank " << tw << ", " << "peel " << nopeel << "\n";

  timer t1;
  t1.start();
  tuple<bool,long> use_v_tuple = cmpWedgeCounts(GA);
  bool use_v = get<0>(use_v_tuple);
  long num_wedges = get<1>(use_v_tuple);
  t1.reportTotal("preprocess (wedge counts)");

  // if you only want total counts
  if (total) {
    timer t;
    t.start();
    long num_butterflies = CountTotal(GA, use_v, num_wedges, max_wedges, max_array_size, ty, tw);
    t.stop();

    if (ty==0) t.reportTotal("Sort:");
    else if (ty==2) t.reportTotal("Hash:");
    else if (ty==4) t.reportTotal("Hist:");
    else if (ty==7 || ty == 11) t.reportTotal("Par");
    else if (ty==8) t.reportTotal("WedgePar");
    cout << "number of butterflies: " << num_butterflies << "\n";
    return;
  }

 //TODO seq code integrate w/count
  if (te == 0) {
    // 12 for work efficient serial
    if (ty == 9) CountOrigCompactSerial(GA,use_v);

    if (ty == 9) return;
    const intT eltsPerCacheLine = 64/sizeof(long);

    timer t;
    t.start();
    long* butterflies = Count(GA, use_v, num_wedges, max_wedges, max_array_size, ty, tw);
    t.stop();

    if (ty==0) t.reportTotal("Sort:");
    else if (ty==1) t.reportTotal("SortCE:");
    else if (ty==2) t.reportTotal("Hash:");
    else if (ty==3) t.reportTotal("HashCE:");
    else if (ty==4) t.reportTotal("Hist:");
    else if (ty==6) t.reportTotal("HistCE:");
    else if (ty==7 || ty == 11) t.reportTotal("Par");
    else if (ty==8) t.reportTotal("WedgePar");
    else if (ty==12) t.reportTotal("Serial");

    if (ty==12) return;

  
    long num_idxs = use_v ? GA.nu : GA.nv;

    auto butterflies_extract_f = [&] (const long i) -> const long {
      return butterflies[i*eltsPerCacheLine];
    };
    
    long total = sequence::reduce<long>((long)0,(long)num_idxs,addF<long>(),butterflies_extract_f);
    
    cout << "number of butterflies: " << total/2 << "\n";
  
    //uintE* butterflies2 = Count(GA, use_v, num_wedges, max_wedges, 0, 0);
    //for (long i=0; i < num_idxs; ++i) { assertf(butterflies[eltsPerCacheLine*i] == butterflies2[eltsPerCacheLine*i], "%d, %d, %d", i, butterflies[eltsPerCacheLine*i], butterflies2[eltsPerCacheLine*i]); }
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

      /*long mc = 0;
      for (size_t i=0; i < num_idxs; i++) { mc = std::max(mc, cores[i]); }
      cout << "### Max core: " << mc << endl;*/
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
    if(ty==2) t3.reportTotal("E Hash:");
    else if (ty == 3) t3.reportTotal("E HashCE:");
    else if (ty == 0) t3.reportTotal("E Sort:");
    else if (ty==1) t3.reportTotal("E SortCE:");
    else if (ty==4) t3.reportTotal("E Hist:");
    else if (ty==6) t3.reportTotal("E HistCE:");
    else if (ty == 5 || ty == 11) t3.reportTotal("E Par:");
    else if (ty == 8) t3.reportTotal("E WedgePar:");

    const intT eltsPerCacheLine = 64/sizeof(long);

    auto butterflies_extract_f = [&] (const long i) -> const long {
      return ebutterflies[i*eltsPerCacheLine];
    };
    
    long total = sequence::reduce<long>((long)0,(long)GA.numEdges,addF<long>(),butterflies_extract_f);
    cout << "number of edge butterflies: " << total/4 << "\n";

    //uintE* butterflies2 = CountE(eti, GA, use_v, num_wedges, max_wedges, 0, 0);
    //for (long i=0; i < GA.numEdges; ++i) { assertf(ebutterflies[eltsPerCacheLine*i] == butterflies2[eltsPerCacheLine*i], "%d, %d, %d", i, ebutterflies[eltsPerCacheLine*i], butterflies2[eltsPerCacheLine*i]); }

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
  if (te == 0 || nopeel) {
  bipartiteCSR G = readBipartite(iFile);

  Compute(G,P);
  for(int r=0;r<rounds;r++) {
    Compute(G,P);
  }
  G.del();
  }
  else {
  bipartiteCSR G = readBipartite(iFile);

  Compute(G,P);
  G.del();
  for(int r=0;r<rounds;r++) {
    G = readBipartite(iFile);
    Compute(G,P);
    G.del();
  }
  }
}


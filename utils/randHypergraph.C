// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "parseCommandLine.h"
#include "graphIO.h"
#include "utils.h"
#include "parallel.h"
using namespace benchIO;
using namespace std;

template <class intT>
hyperedgeArray<intT> hyperedgeRandom(long n, long m, long cardinality) {
  long numEdges = m*cardinality;
  edge<intT> *VE = newA(edge<intT>,numEdges);
  edge<intT> *HE = newA(edge<intT>,numEdges);
  {parallel_for(long i=0;i<m;i++) {
      for(long j=0;j<cardinality;j++) {
	ulong offset = i*cardinality+j;
	HE[offset] = edge<intT>(i,hashInt(offset));
	VE[offset] = edge<intT>(hashInt(offset),i);
      }}}
  return hyperedgeArray<intT>(VE,HE,n,m,numEdges,numEdges);
}


//Generates a symmetrized hypergraph with n vertices and m hyperedges,
//each hyperedge containing c vertices.
int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-n <numvertices>] [-m <numhyperedges>] [-c <cardinality>] <outFile>");
  char* fname = P.getArgument(0);
  
  long n = P.getOptionLongValue("-n",100);
  long m = P.getOptionLongValue("-m", 10*n);  
  long cardinality = P.getOptionLongValue("-c", 3);

  hyperedgeArray<uintT> EA = hyperedgeRandom<uintT>(n, m, cardinality);
  hypergraph<uintT> G = hypergraphFromHyperedges<uintT>(EA);
  EA.del();
  writeHypergraphToFile<uintT>(G, fname);
  G.del();
}

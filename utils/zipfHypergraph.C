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

// int zipf(double alpha, int n)
// {
//   static int first = TRUE;      // Static first time flag
//   static double c = 0;          // Normalization constant
//   static double *sum_probs;     // Pre-calculated sum of probabilities
//   double z;                     // Uniform random number (0 < z < 1)
//   int zipf_value;               // Computed exponential value to be returned
//   int    i;                     // Loop counter
//   int low, high, mid;           // Binary-search bounds

//   // Compute normalization constant on first call only
//   if (first == TRUE)
//     {
//       for (i=1; i<=n; i++)
// 	c = c + (1.0 / pow((double) i, alpha));
//       c = 1.0 / c;

//       sum_probs = malloc((n+1)*sizeof(*sum_probs));
//       sum_probs[0] = 0;
//       for (i=1; i<=n; i++) {
// 	sum_probs[i] = sum_probs[i-1] + c / pow((double) i, alpha);
//       }
//       first = FALSE;
//     }

//   // Pull a uniform random number (0 < z < 1)
//     do
//       {
// 	z = rand_val(0);
//       }
//     while ((z == 0) || (z == 1));

//     // Map z to the value
//     low = 1, high = n, mid;
//     do {
//       mid = floor((low+high)/2);
//       if (sum_probs[mid] >= z && sum_probs[mid-1] < z) {
// 	zipf_value = mid;
// 	break;
//       } else if (sum_probs[mid] >= z) {
// 	high = mid-1;
//       } else {
// 	low = mid+1;
//       }
//     } while (low <= high);

//     // Assert that zipf_value is between 1 and N
//     assert((zipf_value >=1) && (zipf_value <= n));

//     return(zipf_value);
//}

struct edgeFirstCmp {
  bool operator() (edge<uintT> e1, edge<uintT> e2) {
    return e1.u < e2.u;
  }
};

struct edgeSecondCmp {
  bool operator() (edge<uintT> e1, edge<uintT> e2) {
    return e1.v < e2.v;
  }
};

template <class intT>
struct nonNeg {bool operator() (edge<intT> e) {return (e.u != UINT_T_MAX && e.v != UINT_T_MAX);}};

intT binSearch(double* sumProbs, double x, intT nv) {
  intT low = 0, high = nv-1, mid;
  do {
    mid = (low+high)/2;
    if(sumProbs[mid] >= x && sumProbs[mid-1] < x) {
      return mid-1;
    } else if(sumProbs[mid] >= x) {
      high = mid;
    } else {
      low = mid+1;
    }
  } while (low <= high);
  return low-1;
}

template <class intT>
hyperedgeArray<intT> hyperedgeZipf(long nv, long nh, long cardinality, double alpha) {
  long numEdges = nh*cardinality;
  edge<intT> *VE = newA(edge<intT>,numEdges);
  edge<intT> *HE = newA(edge<intT>,numEdges);

  double* sumProbs = newA(double,nv);
  {parallel_for(long i=0;i<nv;i++) sumProbs[i] = 1/pow((double)(1+i),alpha);}
  double c = 1/sequence::plusReduce(sumProbs,nv);
  {parallel_for(long i=0;i<nv;i++) sumProbs[i] *= c; }
  double total = sequence::plusScan(sumProbs,sumProbs,nv);
  //cout << c << " " << total << endl;

  {parallel_for(long i=0;i<nh;i++) {
      for(long j=0;j<cardinality;j++) {
	ulong offset = i*cardinality+j;
	double r = (double) hashInt(offset) / ULONG_MAX;
	intT ngh = binSearch(sumProbs,r,nv);
	HE[offset] = edge<intT>(i,ngh);
	VE[offset] = edge<intT>(ngh,i);
      }
      //need to remove duplicates
      quickSort(VE+i*cardinality,cardinality,edgeFirstCmp());
      quickSort(HE+i*cardinality,cardinality,edgeSecondCmp());
      intT curr = HE[i*cardinality].v;
      for(long j=1;j<cardinality;j++) {
	ulong offset = i*cardinality+j;
	if(HE[offset].v == curr) {HE[offset].v = UINT_T_MAX; VE[offset].u = UINT_T_MAX; }
	else curr = HE[offset].v;
      }
    }}
  free(sumProbs);
  //filter out -1's
  edge<intT> *HE2 = newA(edge<intT>,numEdges);
  intT mh = sequence::filter(HE,HE2,numEdges,nonNeg<intT>());
  free(HE);
  edge<intT> *VE2 = newA(edge<intT>,numEdges);
  intT mv = sequence::filter(VE,VE2,numEdges,nonNeg<intT>());
  free(VE);
  //cout << mv << " " << mh << endl;
  return hyperedgeArray<intT>(VE2,HE2,nv,nh,mv,mh);
}

//Generates a symmetrized hypergraph with n vertices and m hyperedges,
//each hyperedge containing c vertices drawn randomly from a Zipfian
//distribution. Duplicate vertices in a hyperedge are removed.
int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-nv <numvertices>] [-nh <numhyperedges>] [-c <cardinality>] [-a <alpha>] <outFile>");
  char* fname = P.getArgument(0);
  
  long nv = P.getOptionLongValue("-nv",100);
  long nh = P.getOptionLongValue("-nh", 10*nv);  
  long cardinality = P.getOptionLongValue("-c", 3);
  double alpha = P.getOptionDoubleValue("-a",2);
  
  hyperedgeArray<uintT> EA = hyperedgeZipf<uintT>(nv, nh, cardinality, alpha);
  hypergraph<uintT> G = hypergraphFromHyperedges<uintT>(EA);
  EA.del();
  writeHypergraphToFile<uintT>(G, fname);
  G.del();
}

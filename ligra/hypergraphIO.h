// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
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

#pragma once
#include "IO.h"
using namespace std;

template <class vertex>
hypergraph<vertex> readHypergraphFromFile(char* fname, bool isSymmetric, bool mmap) {
  words W;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    char *bytes = newA(char, S.n);
    // Cannot mutate the graph unless we copy.
    parallel_for(size_t i=0; i<S.n; i++) {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  } else {
    _seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n);
  }
#ifndef WEIGHTED
  if (W.Strings[0] != (string) "AdjacencyHypergraph") {
#else
  if (W.Strings[0] != (string) "WeightedAdjacencyHypergraph") {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m -1;
  long n_v = atol(W.Strings[1]);
  long m_v = atol(W.Strings[2]);
  long n_h = atol(W.Strings[3]);
  long m_h = atol(W.Strings[4]);

#ifndef WEIGHTED
  if (len != n_v + m_v + n_h + m_h + 4) {
#else
    if (len != n_v + 2*m_v + n_h + 2*m_h + 4) {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  uintT* offsetsV = newA(uintT,n_v);
  uintT* offsetsH = newA(uintT,n_h);
#ifndef WEIGHTED
  uintE* edgesV = newA(uintE,m_v);
  uintE* edgesH = newA(uintE,m_h);
#else
  intE* edgesV = newA(intE,2*m_v);
  intE* edgesH = newa(intE,2*m_h);
#endif

  {parallel_for(long i=0; i < n_v; i++) offsetsV[i] = atol(W.Strings[i + 5]);}
    {parallel_for(long i=0; i<m_v; i++) {
#ifndef WEIGHTED
      edgesV[i] = atol(W.Strings[i+n_v+5]);
#else
      edgesV[2*i] = atol(W.Strings[i+n_v+5]);
      edgesV[2*i+1] = atol(W.Strings[i+n_v+m_v+5]);
#endif
    }}

#ifndef WEIGHTED
    {parallel_for(long i=0; i < n_h; i++) offsetsH[i] = atol(W.Strings[i + n_v + m_v + 5]);}
#else
    {parallel_for(long i=0; i < n_h; i++) offsetsH[i] = atol(W.Strings[i + n_v + 2*m_v + 5]);}
#endif
    {parallel_for(long i=0; i<m_h; i++) {
#ifndef WEIGHTED
      edgesH[i] = atol(W.Strings[i+n_v+m_v+n_h+5]);
#else
      edgesH[2*i] = atol(W.Strings[i+n_v+2*m_v+n_h+5]);
      edgesH[2*i+1] = atol(W.Strings[i+n_v+2*m_v+n_h+m_h+5]);
#endif
    }}

  //W.del(); // to deal with performance bug in malloc

  vertex* v = newA(vertex,n_v);
  vertex* h = newA(vertex,n_h);

  //vertices
  {parallel_for (uintT i=0; i < n_v; i++) {
    uintT o = offsetsV[i];
    uintT l = ((i == n_v-1) ? m_v : offsetsV[i+1])-offsetsV[i];
    v[i].setOutDegree(l);
#ifndef WEIGHTED
    v[i].setOutNeighbors(edgesV+o);
#else
    v[i].setOutNeighbors(edgesV+2*o);
#endif
    }}

  //hyperedges
  {parallel_for (uintT i=0; i < n_h; i++) {
    uintT o = offsetsH[i];
    uintT l = ((i == n_h-1) ? m_h : offsetsH[i+1])-offsetsH[i];
    h[i].setOutDegree(l);
#ifndef WEIGHTED
    h[i].setOutNeighbors(edgesH+o);
#else
    h[i].setOutNeighbors(edgesH+2*o);
#endif
    }}

  if(!isSymmetric) {
    //in-edges for vertices obtained from out-edges for hyperedges
    uintT* tOffsets = newA(uintT,n_v);
    {parallel_for(long i=0;i<n_v;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair,m_h);
#else
    intTriple* temp = newA(intTriple,m_h);
#endif
    {parallel_for(long i=0;i<n_h;i++){
      uintT o = offsetsH[i];
      for(uintT j=0;j<h[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(h[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(h[i].getOutNeighbor(j),make_pair(i,h[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsetsH);

#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m_h,n_v+1,getFirst<uintE>());
#else
    quickSort(temp,m_h,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m_h,n_v+1,getFirst<intPair>());
#else
    quickSort(temp,m_h,pairFirstCmp<intPair>());
#endif
#endif

    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdgesV = newA(uintE,m_h);
    inEdgesV[0] = temp[0].second;
#else
    intE* inEdgesV = newA(intE,2*m_h);
    inEdgesV[0] = temp[0].second.first;
    inEdgesV[1] = temp[0].second.second;
#endif
    {parallel_for(long i=1;i<m_h;i++) {
#ifndef WEIGHTED
      inEdgesV[i] = temp[i].second;
#else
      inEdgesV[2*i] = temp[i].second.first;
      inEdgesV[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n_v,minF<uintT>(),(uintT)m_h);

    {parallel_for(long i=0;i<n_v;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n_v-1) ? m_h : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
#ifndef WEIGHTED
      v[i].setInNeighbors(inEdgesV+o);
#else
      v[i].setInNeighbors(inEdgesV+2*o);
#endif
      }}

    free(tOffsets);

    //in-edges for hyperedges obtained from out-edges for vertices
    tOffsets = newA(uintT,n_h);
    {parallel_for(long i=0;i<n_h;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    temp = newA(intPair,m_v);
#else
    temp = newA(intTriple,m_v);
#endif
    {parallel_for(long i=0;i<n_v;i++){
      uintT o = offsetsV[i];
      for(uintT j=0;j<v[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsetsV);

#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m_v,n_h+1,getFirst<uintE>());
#else
    quickSort(temp,m_v,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m_v,n_h+1,getFirst<intPair>());
#else
    quickSort(temp,m_v,pairFirstCmp<intPair>());
#endif
#endif

    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdgesH = newA(uintE,m_v);
    inEdgesH[0] = temp[0].second;
#else
    intE* inEdgesH = newA(intE,2*m_v);
    inEdgesH[0] = temp[0].second.first;
    inEdgesH[1] = temp[0].second.second;
#endif
    {parallel_for(long i=1;i<m_v;i++) {
#ifndef WEIGHTED
      inEdgesH[i] = temp[i].second;
#else
      inEdgesH[2*i] = temp[i].second.first;
      inEdgesH[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n_h,minF<uintT>(),(uintT)m_v);

    {parallel_for(long i=0;i<n_h;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n_h-1) ? m_v : tOffsets[i+1])-tOffsets[i];
      h[i].setInDegree(l);
#ifndef WEIGHTED
      h[i].setInNeighbors(inEdgesH+o);
#else
      h[i].setInNeighbors(inEdgesH+2*o);
#endif
      }}

    free(tOffsets);

    Uncompressed_Mem_Hypergraph<vertex>* mem =
      new Uncompressed_Mem_Hypergraph<vertex>(v,h,n_v,m_v,n_h,m_h,edgesV,edgesH,inEdgesV,inEdgesH);
    return hypergraph<vertex>(v,h,n_v,m_v,n_h,m_h,mem);
  }
  else {
    free(offsetsV); free(offsetsH);
    Uncompressed_Mem_Hypergraph<vertex>* mem =
      new Uncompressed_Mem_Hypergraph<vertex>(v,h,n_v,m_v,n_h,m_h,edgesV,edgesH);
    return hypergraph<vertex>(v,h,n_v,m_v,n_h,m_h,mem);
  }
}

/* template <class vertex> */
/* graph<vertex> readGraphFromBinary(char* iFile, bool isSymmetric) { */
/*   char* config = (char*) ".config"; */
/*   char* adj = (char*) ".adj"; */
/*   char* idx = (char*) ".idx"; */
/*   char configFile[strlen(iFile)+strlen(config)+1]; */
/*   char adjFile[strlen(iFile)+strlen(adj)+1]; */
/*   char idxFile[strlen(iFile)+strlen(idx)+1]; */
/*   *configFile = *adjFile = *idxFile = '\0'; */
/*   strcat(configFile,iFile); */
/*   strcat(adjFile,iFile); */
/*   strcat(idxFile,iFile); */
/*   strcat(configFile,config); */
/*   strcat(adjFile,adj); */
/*   strcat(idxFile,idx); */

/*   ifstream in(configFile, ifstream::in); */
/*   long n; */
/*   in >> n; */
/*   in.close(); */

/*   ifstream in2(adjFile,ifstream::in | ios::binary); //stored as uints */
/*   in2.seekg(0, ios::end); */
/*   long size = in2.tellg(); */
/*   in2.seekg(0); */
/* #ifdef WEIGHTED */
/*   long m = size/(2*sizeof(uint)); */
/* #else */
/*   long m = size/sizeof(uint); */
/* #endif */
/*   char* s = (char *) malloc(size); */
/*   in2.read(s,size); */
/*   in2.close(); */
/*   uintE* edges = (uintE*) s; */

/*   ifstream in3(idxFile,ifstream::in | ios::binary); //stored as longs */
/*   in3.seekg(0, ios::end); */
/*   size = in3.tellg(); */
/*   in3.seekg(0); */
/*   if(n != size/sizeof(intT)) { cout << "File size wrong\n"; abort(); } */

/*   char* t = (char *) malloc(size); */
/*   in3.read(t,size); */
/*   in3.close(); */
/*   uintT* offsets = (uintT*) t; */

/*   vertex* v = newA(vertex,n); */
/* #ifdef WEIGHTED */
/*   intE* edgesAndWeights = newA(intE,2*m); */
/*   {parallel_for(long i=0;i<m;i++) { */
/*     edgesAndWeights[2*i] = edges[i]; */
/*     edgesAndWeights[2*i+1] = edges[i+m]; */
/*     }} */
/*   //free(edges); */
/* #endif */
/*   {parallel_for(long i=0;i<n;i++) { */
/*     uintT o = offsets[i]; */
/*     uintT l = ((i==n-1) ? m : offsets[i+1])-offsets[i]; */
/*       v[i].setOutDegree(l); */
/* #ifndef WEIGHTED */
/*       v[i].setOutNeighbors((uintE*)edges+o); */
/* #else */
/*       v[i].setOutNeighbors(edgesAndWeights+2*o); */
/* #endif */
/*     }} */

/*   if(!isSymmetric) { */
/*     uintT* tOffsets = newA(uintT,n); */
/*     {parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;} */
/* #ifndef WEIGHTED */
/*     intPair* temp = newA(intPair,m); */
/* #else */
/*     intTriple* temp = newA(intTriple,m); */
/* #endif */
/*     {parallel_for(intT i=0;i<n;i++){ */
/*       uintT o = offsets[i]; */
/*       for(uintT j=0;j<v[i].getOutDegree();j++){ */
/* #ifndef WEIGHTED */
/* 	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i); */
/* #else */
/* 	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j))); */
/* #endif */
/*       } */
/*       }} */
/*     free(offsets); */
/* #ifndef WEIGHTED */
/* #ifndef LOWMEM */
/*     intSort::iSort(temp,m,n+1,getFirst<uintE>()); */
/* #else */
/*     quickSort(temp,m,pairFirstCmp<uintE>()); */
/* #endif */
/* #else */
/* #ifndef LOWMEM */
/*     intSort::iSort(temp,m,n+1,getFirst<intPair>()); */
/* #else */
/*     quickSort(temp,m,pairFirstCmp<intPair>()); */
/* #endif */
/* #endif */
/*     tOffsets[temp[0].first] = 0; */
/* #ifndef WEIGHTED */
/*     uintE* inEdges = newA(uintE,m); */
/*     inEdges[0] = temp[0].second; */
/* #else */
/*     intE* inEdges = newA(intE,2*m); */
/*     inEdges[0] = temp[0].second.first; */
/*     inEdges[1] = temp[0].second.second; */
/* #endif */
/*     {parallel_for(long i=1;i<m;i++) { */
/* #ifndef WEIGHTED */
/*       inEdges[i] = temp[i].second; */
/* #else */
/*       inEdges[2*i] = temp[i].second.first; */
/*       inEdges[2*i+1] = temp[i].second.second; */
/* #endif */
/*       if(temp[i].first != temp[i-1].first) { */
/* 	tOffsets[temp[i].first] = i; */
/*       } */
/*       }} */
/*     free(temp); */
/*     //fill in offsets of degree 0 vertices by taking closest non-zero */
/*     //offset to the right */
/*     sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m); */
/*     {parallel_for(long i=0;i<n;i++){ */
/*       uintT o = tOffsets[i]; */
/*       uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i]; */
/*       v[i].setInDegree(l); */
/* #ifndef WEIGHTED */
/*       v[i].setInNeighbors((uintE*)inEdges+o); */
/* #else */
/*       v[i].setInNeighbors((intE*)(inEdges+2*o)); */
/* #endif */
/*       }} */
/*     free(tOffsets); */
/* #ifndef WEIGHTED */
/*     Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges,inEdges); */
/*     return graph<vertex>(v,n,m,mem); */
/* #else */
/*     Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights,inEdges); */
/*     return graph<vertex>(v,n,m,mem); */
/* #endif */
/*   } */
/*   free(offsets); */
/* #ifndef WEIGHTED */
/*   Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges); */
/*   return graph<vertex>(v,n,m,mem); */
/* #else */
/*   Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights); */
/*   return graph<vertex>(v,n,m,mem); */
/* #endif */
/* } */

template <class vertex>
hypergraph<vertex> readHypergraph(char* iFile, bool compressed, bool symmetric, bool binary, bool mmap) {
  //if(binary) return readGraphFromBinary<vertex>(iFile,symmetric);
  //else
  return readHypergraphFromFile<vertex>(iFile,symmetric,mmap);
}

template <class vertex>
hypergraph<vertex> readCompressedHypergraph(char* iFile, bool symmetric, bool mmap) {
  //if(binary) return readGraphFromBinary<vertex>(iFile,symmetric);
  //else
  return readHypergraphFromFile<vertex>(iFile,symmetric,mmap);
}

/* template <class vertex> */
/* graph<vertex> readCompressedGraph(char* fname, bool isSymmetric, bool mmap) { */
/*   char* s; */
/*   if (mmap) { */
/*     _seq<char> S = mmapStringFromFile(fname); */
/*     // Cannot mutate graph unless we copy. */
/*     char *bytes = newA(char, S.n); */
/*     parallel_for(size_t i=0; i<S.n; i++) { */
/*       bytes[i] = S.A[i]; */
/*     } */
/*     if (munmap(S.A, S.n) == -1) { */
/*       perror("munmap"); */
/*       exit(-1); */
/*     } */
/*     s = bytes; */
/*   } else { */
/*     ifstream in(fname,ifstream::in |ios::binary); */
/*     in.seekg(0,ios::end); */
/*     long size = in.tellg(); */
/*     in.seekg(0); */
/*     cout << "size = " << size << endl; */
/*     s = (char*) malloc(size); */
/*     in.read(s,size); */
/*     in.close(); */
/*   } */

/*   long* sizes = (long*) s; */
/*   long n = sizes[0], m = sizes[1], totalSpace = sizes[2]; */

/*   cout << "n = "<<n<<" m = "<<m<<" totalSpace = "<<totalSpace<<endl; */
/*   cout << "reading file..."<<endl; */

/*   uintT* offsets = (uintT*) (s+3*sizeof(long)); */
/*   long skip = 3*sizeof(long) + (n+1)*sizeof(intT); */
/*   uintE* Degrees = (uintE*) (s+skip); */
/*   skip+= n*sizeof(intE); */
/*   uchar* edges = (uchar*)(s+skip); */

/*   uintT* inOffsets; */
/*   uchar* inEdges; */
/*   uintE* inDegrees; */
/*   if(!isSymmetric){ */
/*     skip += totalSpace; */
/*     uchar* inData = (uchar*)(s + skip); */
/*     sizes = (long*) inData; */
/*     long inTotalSpace = sizes[0]; */
/*     cout << "inTotalSpace = "<<inTotalSpace<<endl; */
/*     skip += sizeof(long); */
/*     inOffsets = (uintT*) (s + skip); */
/*     skip += (n+1)*sizeof(uintT); */
/*     inDegrees = (uintE*)(s+skip); */
/*     skip += n*sizeof(uintE); */
/*     inEdges = (uchar*)(s + skip); */
/*   } else { */
/*     inOffsets = offsets; */
/*     inEdges = edges; */
/*     inDegrees = Degrees; */
/*   } */


/*   vertex *V = newA(vertex,n); */
/*   parallel_for(long i=0;i<n;i++) { */
/*     long o = offsets[i]; */
/*     uintT d = Degrees[i]; */
/*     V[i].setOutDegree(d); */
/*     V[i].setOutNeighbors(edges+o); */
/*   } */

/*   if(sizeof(vertex) == sizeof(compressedAsymmetricVertex)){ */
/*     parallel_for(long i=0;i<n;i++) { */
/*       long o = inOffsets[i]; */
/*       uintT d = inDegrees[i]; */
/*       V[i].setInDegree(d); */
/*       V[i].setInNeighbors(inEdges+o); */
/*     } */
/*   } */

/*   cout << "creating graph..."<<endl; */
/*   Compressed_Mem<vertex>* mem = new Compressed_Mem<vertex>(V, s); */

/*   graph<vertex> G(V,n,m,mem); */
/*   return G; */
/* } */

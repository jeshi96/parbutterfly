ParButterfly: Parallel algorithms for butterfly computations
======================

Organization
--------

The code for ParButterfly is located in the apps/hyper/ directory, which is
also where compilation should be performed. Graph utilities, particularly
for reading and storing graphs into the proper format, are provided in the
utils/ directory.

Compilation
--------

Compilation is done from within the apps/hyper/ directory.

Compilers

* g++ &gt;= 5.3.0 with support for Cilk Plus
* g++ &gt;= 5.3.0 with OpenMP
* Intel icpc compiler

To compile with g++ using Cilk Plus, define the environment variable
CILK. To compile with icpc, define the environment variable MKLROOT
and make sure CILK is not defined. Using Cilk Plus seems to give the best
parallel performance in our experience.

After the appropriate environment variables are set, to compile,
simply run

```
$ make -j  #compiles with all threads
```

The following commands cleans the directory:
```
$ make clean #removes all executables
$ make cleansrc #removes all executables and linked files from the ligra/ directory
```

Running code
-------
The application BipartiteButterflies in apps/hyper/ takes as input a bipartite
graph in Ligra format, which is described in more detail in the next section.
It also takes parameters for counting and peeling, as follows:

| Parameter&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| Default     | Description                                      |
| ---------    | -------     | ------------------------------------------------ |
| `-countType` | SERIAL      | Wedge/butterfly aggregation type for counting (ASORT, SORT, AHASH, HASH, AHIST, HIST, BATCHS, BATCHWA, or SERIAL) |
| `-rankType`  | ADEG        | Vertex ordering (SIDE, COCORE, ACOCORE, DEG, or ADEG). The prefix A means approximate.                |
| `-peelType`  | NONE        | Butterfly aggregation type for peeling (SORT, HASH, HIST, BATCHS, BATCHWA, or NONE). NONE means don't execute peeling. |
| `-per`       | VERT        | Denotes whether to count per vertex (VERT), per edge (EDGE), or in total (TOTAL). If TOTAL is selected, then -peelType must be NONE. |
| `-m`         | 2577500000  | Maximum number of wedges that the system can hold in memory (for wedge aggregation methods apart from batching) |
| `-a`         | 23090996160 | Maximum number of vertices that the system can hold in memory (for wedge aggregation by batching) |
| `-sparseType`| NONE        | Sparsification options (EDGE, COLOR, NONE) |
| `-d`         | 25          | If -sparseType is EDGE, then 1/d is the probability of keeping an edge. If -sparseType is COLOR, then d is the number of colors. |
| `-r`         | 3           | Number of rounds to repeat butterfly counting/peeling. |

For example:

```
$ ./BipartiteButterflies -countType BATCHS -rankType ACOCORE -peelType BATCHS -per EDGE <input file>
``` 

On NUMA machines, adding the command "numactl -i all " when running
the program may improve performance for large graphs. For example:

```
$ numactl -i all ./BipartiteButterflies <input file>
``` 


Input Format for ParButterfly
-----------
The input format of bipartite graphs is based on the adjacenecy graph format
from the Problem Based Benchmark Suite and Ligra. The adjacency graph format
starts with a sequence of offsets one for each vertex in one bipartition V,
 followed by a sequence of directed edges ordered by their source
 vertex. The offset for a vertex i refers to the location of the start
 of a contiguous block of out edges for vertex i in the sequence of
 edges. The block continues until the offset of the next vertex, or
 the end if i is the last vertex. All vertices and offsets are 0 based
 and represented in decimal. This then repeats for the second bipartition U.

 Let nv and nu denote the number of vertices in bipartitions V and U respectively. 
 Let mv and mu denote the number of edges (mv = mu). 
 The specific format is as follows:

AdjacencyHypergraph  
&lt;nv>  
&lt;mv>  
&lt;nu>  
&lt;mu>  
&lt;offsetv(0)>   
&lt;offsetv(1)>  
...  
&lt;offsetv(nv-1)>  
&lt;edgev(0)>  
&lt;edgev(1)>  
...  
&lt;edgev(mv-1)>   
&lt;offsetu(0)>  
&lt;offsetu(1)>  
...  
&lt;offsetu(nv-1)>  
&lt;edgeu(0)>  
&lt;edgeu(1)>  
...  
&lt;edgeu(mv-1)> 

Graph Utilities
---------

Several graph utilities are provided in the utils/ directory and can
be compiled using "make". Importantly, ParButterfly works only with 
hypergraph formats as described in the previous section.

**KONECTtoHyperAdj** converts a graph in [KONECT
format](http://konect.cc/) and converts it to Ligra's
adjacency hypergraph format. The first required parameter is the input
(KONECT) file name and second required parameter is the output (Ligra)
file name.

Optimizations
---------
Independently of our work, [Wang et 
al.](http://www.vldb.org/pvldb/vol12/p1139-wang.pdf) described an 
algorithm for butterfly counting using degree ordering, and also 
proposed a cache-aware wedge processing optimization that changes 
the direction in which wedges are inspected. 
Since our initial publication and inspired by their proposed optimization, we have 
incorporated the cache-aware wedge processing optimization of Wang et al.
into our framework for all orderings
(previously we had only done for the SIDE ordering just by chance). 
This optimization can be turned on using
```
#define INVERSE 1
```
If INVERSE is not defined, then the optimization will be turned off.

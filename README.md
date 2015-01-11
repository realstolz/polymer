ABOUT POLYMER
=======

Polymer is a NUMA-aware Graph-structured Analytics Framework written in C++. It adopts a general design principle for NUMA machines by co-locating graph data and computation within NUMA-nodes as much as possible, with the goal of reducing remote memory accesses and balancing cross-node interconnect bandwidth. 

Specifically, Polymer differentially allocates and places graph data according to their access patterns. For graph topology data such as vertices and edges that are always accessed only their own threads, Polymer uses corporative allocation by letting an accessing thread to allocate the memory in its local memory node, and thus eliminating remote accesses. 

Second, for application-defined data (like the ranks of a web page in PageRank) whose memory locations are static but data will be updated dynamically, though corporative allocation may eliminate a lot of remote accesses, there are still inevitable remote accesses due to frequent exchanges of application-defined data during computation. Hence, Polymer allocates such data with contiguous virtual addresses, but distributes actual physical memory frames to the NUMA-node of the owning thread. This makes it seamless to access cross-node data. For mutable graph runtime states such as the current active vertices, as they are dynamically allocated in each iteration, Polymer allocates and updates such data in a distributed way but accesses it through a global lookup table to avoid contention.

Currently, Polymer has implemeted six popular graph algorithms files, including PageRank (numa-PageRank.C), SpMV (numa-SPMV.C), BP (numa-BP.C), BFS (numa-BFS.C), CC (numa-Componenet.C) and SSSP (numa-BellmanFord.C). To compile all of them, simply run "make" with the appropriate environment variables set as described above. The results of the computation are not used, but the code can be easily modified to output the results to a file.

To develop a new implementation, simply include "polymer.h" in the implementation files. When finished, one may add it to the ALL variable in Makefile.


LICENSE
=======

Polymwer is released under the Apache 2 license. Please see the file called LICENSE.

If you use GraphLab PowerGraph in your research, please cite our paper:
Kaiyuan Zhang, Rong Chen and Haibo Chen. NUMA-aware Graph-structured Analytics. ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP'15), Bay Area, California, USA, February, 2015. 


BUILD & RUN
=======

Polymer compiles with g++ version 4.8.0 or higher with support for Cilk+. To compile with g++ using Cilk, define the environment variable CILK. To compile with g++ with no parallel support, make sure CILK is not defined.

With correct version of g++ installed (Cilk+ recommended), use
below command to compile all alogrithms.
```
make
```

After that use the following command to run the algorithm:
```
PageRank: ./numa-PageRank [graph file] [maximum iteration]
SPMV: ./numa-SPMV [graph file] [maximum iteration]
BP: ./numa-BP [graph file] [maximum iteration]
BFS: ./numa-BFS [graph file] [start vertex number]
BellmanFord: ./numa-BellmanFord [graph file] [start vertex number]
ConnectedComponents: ./numa-Components [graph file]
```      

INPUT FORMAT
=======

The input format of an unweighted graphs should be in following format.

The adjacency graph format from the Problem Based Benchmark Suite (http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html). The adjacency graph format starts with a sequence of offsets one for each vertex, followed by a sequence of directed edges ordered by their source vertex. The offset for a vertex i refers to the location of the start of a contiguous block of out edges for vertex i in the sequence of edges. The block continues until the offset of the next vertex, or the end if i is the last vertex. All vertices and offsets are 0 based and represented in decimal. The specific format is as follows (NOTE: This file is represented as plain text.):
```
 AdjacencyGraph
 <n>
 <m>
 <o0>
 <o1>
 ...
 <o(n-1)>
 <e0>
 <e1>
 ...
 <e(m-1)>
```

CONTACT
=======

Website: http://ipads.se.sjtu.edu.cn/projects/polymer.html

Kaiyuan Zhang <johnzh@sjtu.edu.cn>
Rong Chen <rongchen@sjtu.edu.cn>
Haibo Chen <haibochen@sjtu.edu.cn>


ACKNOWLEDGEMENT
=======

Currently, Polymer mostly follows the interface from [Ligra](https://github.com/jshun/ligra). 


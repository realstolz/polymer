cmake_minimum_required(VERSION 3.6)
project(polymer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    micro-bench/barrier-bench.cc
    micro-bench/rw-cycle-bench.cc
    micro-bench/test-prefetch.cc
    micro-bench/two-thread-read.cc
    micro-bench/two-thread-write.cc
    ConvertTmp.C
    ConvertToBinary.C
    ConvertToJSON.C
    custom-barrier.h
    DegreeCount.C
    gettime.h
    graph.h
    IO-numa.h
    IO.h
    ligra-numa.h
    ligra.h
    numa-BellmanFord.C
    numa-BFS-async-pipe.C
    numa-BFS.C
    numa-BP.C
    numa-Components.C
    numa-PageRank-bin.C
    numa-PageRank-pull.C
    numa-PageRank-write.C
    numa-PageRank.C
    numa-PageRankDelta.C
    numa-SPMV.C
    parallel.h
    polymer-wgh.h
    polymer.h
    quickSort.h
    utils-numa.h
    utils.h)

add_executable(polymer ${SOURCE_FILES})
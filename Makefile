ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

# # no compare and swap!
# ifdef OPENMP
# PCC = g++
# PCFLAGS = -fopenmp -mcx16 -O3 -DOPENMP $(INTT) $(INTE)

ifdef CILK
PCC = g++
#-cilk
PCFLAGS = -fcilkplus -lcilkrts -O2 -DCILK $(INTT) $(INTE)
PLFLAGS = -fcilkplus -lcilkrts

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -O3 -DCILKP $(INTT) $(INTE)

else
PCC = g++
PCFLAGS = -O2 $(INTT) $(INTE)
endif

PCFLAGS = -fcilkplus -lcilkrts -O2 -DCILK $(INTT) $(INTE)
PLFLAGS = -fcilkplus -lcilkrts

COMMON= ligra.h ligra-rewrite.h ligra-rewrite-wgh.h graph.h utils.h IO.h parallel.h gettime.h quickSort.h

ALL= BFS BC Components Radii PageRank PageRankDelta BellmanFord PageRank-Pull 

MYAPPS= DegreeCount SPMV BP numa-BP numa-PageRank numa-PageRank-pull numa-PageRank-write numa-PageRankDelta numa-Components numa-BFS numa-SPMV numa-BellmanFord ConvertToJSON
MYHEADER= ligra-rewrite.h ligra-numa.h
LIBS_I_NEED= -pthread -lnuma

all: $(ALL) $(MYAPPS)

debug: PCFLAGS += -g
debug: all

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $< $(LIBS_I_NEED)

.PHONY : clean

clean :
	rm -f *.o $(ALL) $(MYAPPS)


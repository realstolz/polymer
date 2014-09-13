
CC = g++

LIBS = -pthread -lnuma
BENCHMARKS = two-thread-read two-thread-write rw-cycle-bench barrier-bench test-prefetch

all: $(BENCHMARKS)

% : %.cc
	$(CC) -o $@ $< $(LIBS)

.PHONY: clean

clean:
	rm -f *.o $(BENCHMARKS)

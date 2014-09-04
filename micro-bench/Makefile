
all:
	g++ -o read_bench two-thread-read.cc -pthread -lnuma
	g++ -o write_bench two-thread-write.cc -pthread -lnuma
	g++ -o rw-cycle-bench rw-cycle-bench.cc -pthread -lnuma
	g++ -o barrier-bench barrier-bench.cc -pthread -lnuma

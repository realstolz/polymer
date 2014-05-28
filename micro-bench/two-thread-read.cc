#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <numa.h>
#include <time.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

#define PAGE_SIZE (4 * 1024)
#define NUM_OF_PAGES (1024*100)
int REGION_SIZE = (PAGE_SIZE * NUM_OF_PAGES);

void *mem_pos[2]; //pos of two mem
char* thread_1_nodestring;
char* thread_2_nodestring;
bool randomAccess = false;

struct my_args {
    int my_tid;
    pthread_mutex_t *mut1;
    pthread_mutex_t *mut2;
};

void *sayHello(void *arg_ptr) {
    struct my_args *arg = (struct my_args *)arg_ptr;
    pthread_detach(pthread_self());
    pthread_mutex_t *mut_to_lock;
    pthread_mutex_t *another_one;
    char *cpuString;
    void **myPos;
    void **otherPos;
    int *buf;
    if (arg->my_tid == 0) {
	mut_to_lock = arg->mut1;
	another_one = arg->mut2;
	cpuString = thread_1_nodestring;
	myPos = &(mem_pos[0]);
	otherPos = &(mem_pos[1]);
    } else {
	mut_to_lock = arg->mut2;
	another_one = arg->mut1;
	cpuString = thread_2_nodestring;
	myPos = &(mem_pos[1]);
	otherPos = &(mem_pos[0]);
    }

    struct bitmask *cpumask = numa_parse_nodestring(cpuString);
    numa_bind(cpumask);
    srand(time(NULL));

    struct timezone tz = {0, 0};

    // allocate mem, bind thread to cpu
    int *ptr = (int *)numa_alloc_local(REGION_SIZE);
    int i;
    int numOfInt = REGION_SIZE / sizeof(int);
    for (i = 0; i < (REGION_SIZE / sizeof(int)); i++) {
	ptr[i] = rand() % numOfInt;
    }
    *myPos = (char *)ptr;
    pthread_mutex_unlock(mut_to_lock);
    printf("%d: start computing\n", arg->my_tid);
    pthread_mutex_lock(another_one);
    struct timeval start, end;
    if (arg->my_tid == 0) {
	gettimeofday(&start, &tz);
	// acquire lock to make sure mem been initialized
	ptr = (int *)*otherPos;
	int bufSum = 0;
	for (i = 0; i < (REGION_SIZE / sizeof(int)); i++) {
	    int idx = randomAccess?ptr[i]:i;
	    int val = ptr[idx];
	    bufSum += val;
	}
	gettimeofday(&end, &tz);
	//printf("%d\n", bufSum);
    } else {
	gettimeofday(&start, &tz);
	ptr = (int *)*myPos;
	int bufSum = 0;
	for (i = 0; i < (REGION_SIZE / sizeof(int)); i++) {
	    int idx = randomAccess?ptr[i]:i;
	    int val = ptr[idx];
	    bufSum += val;
	}
	//memcpy((void *)*otherPos, (void *)buf, REGION_SIZE);
	gettimeofday(&end, &tz);
	//printf("%d\n", bufSum);
    }
    double timeStart = ((double)start.tv_sec) + ((double)start.tv_usec) / 1000000.0;
    double timeEnd = ((double)end.tv_sec) + ((double)end.tv_usec) / 1000000.0;
    pthread_mutex_unlock(another_one);
    printf("%d: %lf\n", arg->my_tid, timeEnd - timeStart);
    pthread_exit(0);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
	printf("usage: write_bench [thread_1 node] [thread_2 node] [write_size]\n");
	return 0;
    }
    thread_1_nodestring = argv[1];
    thread_2_nodestring = argv[2];
    REGION_SIZE = atoi(argv[3]);
    randomAccess = false;
    if(argc > 4) if((std::string) argv[4] == (std::string) "-random") randomAccess = true;
    
    int i;
    pthread_mutex_t mut1;
    pthread_mutex_t mut2;
    pthread_mutex_init(&mut1, NULL);
    pthread_mutex_init(&mut2, NULL);
    pthread_mutex_lock(&mut1);
    pthread_mutex_lock(&mut2);
    mem_pos[0] = NULL;
    mem_pos[1] = NULL;
    for (i = 0; i < 2; i++) {
        pthread_t tid;
	struct my_args *arg = (struct my_args *)malloc(sizeof(struct my_args));
	arg->my_tid = i;
	arg->mut1 = &mut1;
	arg->mut2 = &mut2;
        pthread_create(&tid, NULL, sayHello, (void *)arg);
    }
    pthread_exit(NULL);
}

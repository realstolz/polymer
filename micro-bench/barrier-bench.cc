#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <pthread.h>
#include <numa.h>
#include <sys/time.h>
#include "../custom-barrier.h"

int socketIdx[8] = {0, 2, 3, 1, 5, 4, 6, 7};
int numOfNode = 0;
int CORES_PER_NODE = 0;

volatile int global_counter = 0;
volatile int global_toggle = 0;

volatile int global_counter2 = 0;
volatile int global_toggle2 = 0;

Custom_barrier *barriers;

pthread_barrier_t global_barr;
pthread_barrier_t submaster_barr;
pthread_barrier_t *local_barr_list;

void printTime(int subTid, struct timeval startT, struct timeval endT) {
    if (subTid != 0) {
	return;
    }
    double timeStart = (double)startT.tv_usec;
    double timeEnd = (double)endT.tv_usec;
    double barrTime = timeEnd - timeStart;
    printf("barrier time: %lf\n", barrTime);
}

void *threadSubFunc(void *arg) {
    int subTid = *(int *)arg;
    //printf("subTid: %d %d\n", subTid, subTid / CORES_PER_NODE);
    pthread_barrier_t *local_barr = &local_barr_list[subTid / CORES_PER_NODE];
    
    Custom_barrier global_custom(&global_counter, &global_toggle, numOfNode * CORES_PER_NODE);
    Custom_barrier submaster_custom(&global_counter2, &global_toggle2, numOfNode);
    Custom_barrier local_custom = barriers[subTid];

    pthread_barrier_wait(&global_barr); 

    if (subTid == 0) {
	sleep(3);
    }

    struct timeval startT, endT;
    struct timezone tz = {0, 0};
    gettimeofday(&startT, &tz);
    
    pthread_barrier_wait(&global_barr);
    
    gettimeofday(&endT, &tz);
    printTime(subTid, startT, endT);
    
    if (subTid == 0) {
	sleep(3);
    }

    gettimeofday(&startT, &tz);
    pthread_barrier_wait(local_barr);
    if (subTid % CORES_PER_NODE == 0) {
	pthread_barrier_wait(&submaster_barr);
    }
    pthread_barrier_wait(local_barr);
    
    gettimeofday(&endT, &tz);
    printTime(subTid, startT, endT);

    if (subTid == 0) {
	sleep(3);
    }

    gettimeofday(&startT, &tz);
    
    global_custom.wait();
    
    gettimeofday(&endT, &tz);
    printTime(subTid, startT, endT);

    if (subTid == 0) {
	sleep(3);
    }

    gettimeofday(&startT, &tz);
    
    local_custom.wait();
    if (subTid % CORES_PER_NODE == 0) {
	submaster_custom.wait();
    }
    local_custom.wait();
    
    gettimeofday(&endT, &tz);
    printTime(subTid, startT, endT);
    

    return NULL;
}

void *threadFunc(void *arg) {
    int tid = *(int *)arg;
    int node = socketIdx[tid];
    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    printf("%d ok %s\n", tid, nodeString);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);

    pthread_barrier_wait(&submaster_barr);

    pthread_barrier_init(&(local_barr_list[tid]), NULL, CORES_PER_NODE);

    pthread_barrier_wait(&submaster_barr);

    volatile int localCounter;
    volatile int localToggle;

    int sub_args[CORES_PER_NODE]; 
    pthread_t tids[CORES_PER_NODE];
    for (int i = 0; i < CORES_PER_NODE; i++) {
	sub_args[i] = i + tid * CORES_PER_NODE;
	Custom_barrier barr(&localCounter, &localToggle, CORES_PER_NODE);
	barriers[sub_args[i]] = barr;
	pthread_create(&tids[i], NULL, threadSubFunc, (void *)(&sub_args[i]));
    }

    for (int i = 0; i < CORES_PER_NODE; i++) {
	pthread_join(tids[i], NULL);
    }

    return NULL;
}

int main(int argc, char* argv[]) {
    numOfNode = numa_num_configured_nodes();
    CORES_PER_NODE = numa_num_configured_cpus() / numOfNode;
    if (argc > 1) {
	numOfNode = atoi(argv[1]);
    }

    pthread_barrier_init(&global_barr, NULL, numOfNode * CORES_PER_NODE);
    pthread_barrier_init(&submaster_barr, NULL, numOfNode);
    local_barr_list = (pthread_barrier_t *)malloc(sizeof(pthread_barrier_t) * numOfNode);
    
    barriers = (Custom_barrier *)malloc(sizeof(Custom_barrier) * numOfNode * CORES_PER_NODE);
    
    int args[numOfNode]; 
    pthread_t tids[numOfNode];
    for (int i = 0; i < numOfNode; i++) {
	args[i] = i;
	pthread_create(&tids[i], NULL, threadFunc, (void *)(&args[i]));
	printf("%p %d\n", &args[i], args[i]);
    }
    
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    return 0;
}

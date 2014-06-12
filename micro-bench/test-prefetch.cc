#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

int PREFETCH_DIST = 1;

int main(int argc, char *argv[]) {
    int size = atoi(argv[1]);
    if (argc > 2) {
	PREFETCH_DIST = atoi(argv[2]);
    }
    int *arr = (int *)malloc((size+1) * sizeof(int));
    srand(time(NULL));
    for (int i = 0; i < size; i++) {
	arr[i] = rand() % size;
    }

    //start random access
    struct timezone tz = {0, 0};
    struct timeval start, end;
    gettimeofday(&start, &tz);

    int sum = 0;
    for (int i = 0; i < size / 4; i++) {
	arr[(int)arr[i]]++;
	__builtin_prefetch(&arr[(int)arr[i+PREFETCH_DIST]], 1, 3);
	i++;
	arr[(int)arr[i]]++;
	__builtin_prefetch(&arr[(int)arr[i+PREFETCH_DIST]], 1, 3);
	i++;
	arr[(int)arr[i]]++;
	__builtin_prefetch(&arr[(int)arr[i+PREFETCH_DIST]], 1, 3);
	i++;
	arr[(int)arr[i]]++;
	__builtin_prefetch(&arr[(int)arr[i+PREFETCH_DIST]], 1, 3);
    }

    gettimeofday(&end, &tz);
    double timeStart = ((double)start.tv_sec) + ((double)start.tv_usec) / 1000000.0;
    double timeEnd = ((double)end.tv_sec) + ((double)end.tv_usec) / 1000000.0;
    printf("Time used: %lf\n", timeEnd - timeStart);
}

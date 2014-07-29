#include <stdio.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

__inline__ __attribute__((always_inline)) uint64_t rdtsc() {
    uint64_t a;
    __asm__ __volatile__ ("cpuid\n\t"
			  "rdtsc\n\t"
			  : "=A" (a));
    return a;
}

int main() {
    int i = 0;
    uint64_t s, e;
    s = rdtsc();
    printf("first: %" PRIu64 "\n", s);
    i++;
    e = rdtsc();
    printf("rdtsc: %" PRIu64 " %" PRIu64 "\n", e, (int)e -(int)s);
    return i;
}

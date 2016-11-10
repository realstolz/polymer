#ifndef CUSTOM_BARRIER
#define CUSTOM_BARRIER

#include <iostream>
#include <pthread.h>
#include <numa.h>

template<class ET>
inline void toggle(volatile ET *var) {
    *var = 1 - *var;
}

template<class ET>
inline void spinUntilNotEq(volatile ET *var, ET val) {
    while (*var == val) {
        __asm__ __volatile__ ("pause\n\t":: :"memory");
    }
}

struct Custom_barrier {
    volatile int *counter;
    volatile int *toggleVar;
    int count;

    Custom_barrier() {}

    Custom_barrier(volatile int *_counter, volatile int *_toggle, int _count) : counter(_counter), toggleVar(_toggle),
                                                                                count(_count) {}

    inline void wait() {
        int oldToggle = *toggleVar;
        int oldVal = __sync_add_and_fetch(counter, 1);
        if (oldVal >= count) {
            *counter = 0;
            toggle(toggleVar);
        } else {
            spinUntilNotEq(toggleVar, oldToggle);
        }
    }
};

#endif

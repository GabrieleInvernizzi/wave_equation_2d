#ifdef COLLECT_TIMINGS
#include "timings.h"

#include <bits/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


typedef struct {
    char* label;
    struct timespec init_time;
    int is_running;
} Timer;

static Timer t = { 0 };

void start_timer(char* label, int rank) {
    if (t.is_running) {
        fprintf(stderr, "Can't run multiple timers concurrently.\n");
        abort();
    }

    t.is_running = 1;
    t.label = label;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t.init_time);
}

void end_timer() {
    struct timespec end_time;

    if (!t.is_running) {
        fprintf(stderr, "Can't end stopped timer.\n");
        abort();
    }

    t.is_running = 0;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    long diff = (end_time.tv_sec - t.init_time.tv_sec) * (long)1e9 + (end_time.tv_nsec - t.init_time.tv_nsec);
    fprintf(stderr, "%s, %lu\n", t.label, diff);
}

#endif

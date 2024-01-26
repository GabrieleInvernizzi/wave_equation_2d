#pragma once

#ifdef COLLECT_TIMINGS

void start_timer(char *label, int rank);
void end_timer();
void save_timings();

#define START_TIMER(lab, rank) start_timer(lab, rank)
#define END_TIMER end_timer()

#else

#define START_TIMER(lab, rank)
#define END_TIMER

#endif

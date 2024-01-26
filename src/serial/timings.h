#pragma once

#ifdef COLLECT_TIMINGS

void start_timer(char *label);
void end_timer();
void save_timings();

#define START_TIMER(lab) start_timer(lab)
#define END_TIMER end_timer()
#define SAVE_TIMINGS save_timings()

#else

#define START_TIMER(lab)
#define END_TIMER
#define SAVE_TIMINGS

#endif

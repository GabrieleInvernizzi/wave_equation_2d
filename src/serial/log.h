#pragma once

#include <stdio.h>

#define LOG_FILE_DESC (stdout)

#ifdef NOLOG

#define LOGF(str, ...)

#else

#define LOGF(str, ...) fprintf(LOG_FILE_DESC, str "\n", ##__VA_ARGS__)

#endif


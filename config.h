#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "coefs.h"

float* tested_set = bruteforced_coefs;

#define MAX_THREADS (10)
#define MULTITHREADED (1)

#define SAVE_TO_FILE (1)

#define BRUTEFORCE (1)

#include "bruteforce.h"

#define RUN_BRUTEFORCE(body) \
{ \
    BRUTEFORCE_LOOP(0, 2, -10, 10, \
    BRUTEFORCE_LOOP(1, 2, -10, 10, \
    BRUTEFORCE_LOOP(2, 2, -10, 10, \
    BRUTEFORCE_LOOP(3, 2, -10, 10, \
    BRUTEFORCE_LOOP(4, 2, -10, 10, \
    BRUTEFORCE_LOOP(5, 2, -10, 10, \
    { \
        body \
    } \
    )))))) \
}

#endif
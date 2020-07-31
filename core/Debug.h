#ifndef SRC_SHARED_DEBUG_H
#define SRC_SHARED_DEBUG_H

//#define DEBUG

#ifdef DEBUG

#include <cstdio>
#define TRACE(...) printf(__VA_ARGS__)

#else

#define TRACE(...) (void)0

#endif /* DEBUG */

#endif /* SRC_SHARED_DEBUG_H */
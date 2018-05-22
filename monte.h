/* type definitions for monte carlo simulation */

/*
 * HISTORY
 * 16-Oct-88 Jeffrey Trull (jt1j) at Carnegie-Mellon University
 *       Created
 */

/* note: in actual fact I am typing this in, 30 years later, from a paper copy,
 * and being a bit inexact about it
 */

#ifndef MONTE_HDR
#define MONTE_HDR

typedef struct vector_str *vector_ptr;
typedef struct vector_str {
    float x, y, z;
} vector;

#endif // MONTE_HDR

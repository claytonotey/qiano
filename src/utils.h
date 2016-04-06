#ifndef UTILS_H
#define UTILS_H

#ifndef __has_extension
#define __has_extension(x) 0
#endif
#define vImage_Utilities_h
#define vImage_CVUtilities_h

#include <Accelerate/Accelerate.h>

//#define logf(format, ...) fprintf(stderr,format,__VA_ARGS__)
#define logf(format, ...) 

#include "immintrin.h"

typedef float vec4 __attribute__ ((vector_size (sizeof(float) * 4)));
typedef float vec8 __attribute__ ((vector_size (sizeof(float) * 8)));
typedef float vec16 __attribute__ ((vector_size (sizeof(float) * 16)));

#define HALFPI 1.5707963267948966192313216916398
#define PI 3.1415926535897932384626433832795
#define TWOPI 6.283185307179586476925286766559

std::ostream& operator<<(std::ostream& os, const vec4 &v);

float dot(int N, float *A, float *B);
float sum8(__m256 x);
float sse_dot(int N, float *A, float *B);
float dsp_dot(int N, float *A, float *B);
float sum4(vec4 x);
void ms4(float *x, float *y, float *z, int N);
void diff4(float *x, float *y, int N);

static inline float square(float x)
{
  return x * x;
}
 
static inline float cube(float x)
{
  return x * x * x;
}

#endif

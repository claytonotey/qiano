#ifndef FILTER_H
#define FILTER_H

#include "types.h"
#include <string.h>
#include "sse.h"

#ifndef __has_extension
#define __has_extension(x) 0
#endif
#define vImage_Utilities_h
#define vImage_CVUtilities_h

#include <Accelerate/Accelerate.h>

enum {
  MaxFilterUpsample = 8
};

void sse_dotN(int N, float *A, float *B, float *C);
float sse_dot(int N, float *A, float *B);
float dsp_dot(int N, float *A, float *B);
float dot(int N, float *A, float *B);

class Filter {
 public:
  Filter(int n);
  ~Filter();
  float phasedelay(float omega);
  float groupdelay(float omega);
  void merge(const Filter &c1, const Filter &c2);
  inline float filter(float in);
  void init(int upsample = 1);

  float *x;
  float *xc;
  float *xend;
  float *b;
  float *bend;
  int n;
  int nmax;
  int xstep;
  int xskip;
  int upsample;
};


inline float Filter :: filter(float in) 
{
  float *b = this->b;
  float *x = this->xc;
  float out = *(b) * in;  
  b+=2;
  x+=xstep;

  while(b <= bend) {
    if(x>xend) x -= xskip;
    out += *(b) * *(x) - *(b+1) * *(x+1);
    b+=2;
    x+=xstep;
  }
  
  x = this->xc;
  *(x) = in;
  *(x+1) = out;
  x-=2; if(x<this->x) x = xend; this->xc = x;

  return out;
}

class BiquadHP : public Filter 
{
 public:
 BiquadHP() : Filter(2) {}
  void create(float omega, float Q);
};

class DWGResonator {
 public:
  DWGResonator();
  void create(float omega, float gamma);
  float go(float in);

  // variables
  float g;
  float c;
  float b1t;

  // state
  float x1;
  float x2;
};

class ConvolutionResonator {
public:
  ConvolutionResonator();
  ~ConvolutionResonator();

  void create(float omega, float gamma);
  float go(float in);

  float *b;
  float *x;
  float *bend;
  float *xc;
  float *xend;
  int size;
};

class Thiran : public Filter {
 public:
 Thiran(int n) : Filter(n) {}
  void create(float D, int N, int upsample = 1);
};

class ThiranDispersion : public Thiran {
 public:
 ThiranDispersion() : Thiran(2) {}
  void create(float B, float f, int M, int upsample = 1);
};

class Loss : public Filter
{
 public:
 Loss() : Filter(1) {}
  float filter(float in) {
    float out = b[0] * in - b[3] * x[1];
    x[1] = out;
    return out;
  }

  void create(float f0, float c1, float c3, int upsample = 1) {
    f0 /= upsample;
    c1 /= upsample;
    c3 *= upsample;
    float g = 1.0 - c1/f0; 
    float c = 4.0*c3 + f0;
    float a1 = (-c+sqrt(c*c-16.0*c3*c3))/(4.0*c3);
    b[0] = g*(1+a1);
    b[2] = 0;
    b[1] = 1;    
    b[3] = a1;     
    
    n = 1;
    init(1);
   }
};

template <int size>
class Delay
{
 public:
  int di; 
  int d1;
  int cursor;
  int probeOffset;

  float x[size] __attribute__((aligned(32)));
  bool bFull;

  enum {
    mask = size-1 };

  Delay() {
    bFull = false;
    cursor = 0;
  }

  void setDelay(int di) {
    this->di = di;    
    if(bFull) {
      d1 = (cursor+size-di)&(mask);
    } else {
      d1 = cursor - di;
    }
  }

  void clear() {
    memset(x,0,size*sizeof(float));
  }

  float probe() {
    return x[(d1 + size)&mask];
  }

  float goDelay(float in) {
    float y0;
    x[cursor] = in;
    cursor++;
    if(d1 < 0) {
      y0 = 0;
      d1++;
    } else {
      y0 = x[d1];
      d1 = (d1 + 1) & mask;
    }
    if(cursor & size) {
      bFull = true;
      cursor = 0;
    }
    return y0; 
  }
};


class MSDFilter : public Filter {
public:
  MSDFilter() : Filter(2) {}
  void create(float Fs, float m, float k, float mu, float RT);
};

class MSD2Filter {
public:
 MSD2Filter() : f11(4), f12(4), f21(4), f22(4) {}
  void filter(float in[2], float out[2]);
  void create(float Fs,
              float m1, float k1, float R1, 
              float m2, float k2, float R2,
              float R12, float k12, float Z); 
    
  /*  Filter f11;
  Filter f12;
  Filter f21;
  Filter f22;
  */
  float f11,f12,f21,f22;
};

enum {
  DownSampleFilterSize = 64
};

class DownSampleFIR {
 public:
  int getDelay();
  void create(int DownSample);
  float filter(float in);
  float b[DownSampleFilterSize];
  float x[DownSampleFilterSize];
  float *xc;
  float *xend;
  float *bend;
};

#endif

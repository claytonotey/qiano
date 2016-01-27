#include "filter.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "types.h"
#include "utils.h"

#include "filter.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "types.h"


float dot(int N, float *A, float *B) {
  float dot = 0;
  for(int i=0; i<N; i++) {
    dot += A[i] * B[i];
  }
  return dot;
}

#define V8

#ifdef V8
enum {
VSIZE = 8
};
#define BROADCAST(v) _mm256_broadcast_ss(v)
#define FMADD(a,b,c) _mm256_fmadd_ps(a,b,c)
#else
enum {
VSIZE = 4
};
#define BROADCAST(v) _mm_broadcast_ss(v)
#define FMADD(a,b,c) _mm_fmadd_ps(a,b,c)
#endif

typedef float vec __attribute__ ((vector_size (sizeof(float) * VSIZE)));

void sse_dotN(int N, float *A, float *B, float *C) {
  vec temp = {0};

  //  vec *Av = (vec *)A;
  vec *Bv = (vec *)B;
  vec *Cv = (vec *)C;
  
  for(int i = 0; i < N; i++) {
    vec Av =  BROADCAST(A);
    vec Bv0 = *Bv;
    temp = FMADD( Av,Bv0,temp);

    A+=1;
    Bv+=1;
  }
  
  *Cv = temp;
  
}

float sse_dot(int N, float *A, float *B) {


    vec temp0 = {0};
    vec temp1 = {0};
    vec temp2 = {0};
    vec temp3 = {0};
    
    vec *Av = (vec *)A;
    vec *Bv = (vec *)B;

    vec Bv0;
    vec Bv1;
    vec Bv2;
    vec Bv3;


    int N32 = N / (4 * VSIZE);
    N -= 4*VSIZE*N32;
    int N8 =  N / VSIZE;
    N -= VSIZE*N8;

    for(int i = 0; i < N32; i++) {

      Bv0 = *(Bv);
      Bv1 = *(Bv+1);
      Bv2 = *(Bv+2);
      Bv3 = *(Bv+3);
      temp0 = FMADD(*Av,Bv0,temp0);
      temp1 = FMADD(*(Av+1),Bv1,temp1);
      temp2 = FMADD(*(Av+2),Bv2,temp2);
      temp3 = FMADD(*(Av+3),Bv3,temp3);

      Av+=4;
      Bv+=4;
    }

    for(int i = 0; i < N8; i++) {
      temp0 = FMADD(*Av,*Bv,temp0);
      Av++;
      Bv++;
    }

    A = (float*)Av;
    B = (float*)Bv;

    float dot = 0;
    for(int i=0; i<N; i++) {
      dot += A[i] * B[i];
    }

    union {
      vec tempv;
      float tempf[VSIZE];
    };
    
    
    if(N32 || N8) {
      tempv = temp0 + temp1 + temp2 + temp3;

      for(int i = 0; i < VSIZE; ++i) {
        dot += tempf[i];
      }
    }

    return dot;
}

void *_aligned_malloc(size_t size, size_t align)
{
  void *p;
  posix_memalign(&p,align,size);
  return p;
}

void _aligned_free(void *p)
{
  free(p);
}

float *filter_malloc(int size)
{
  size = 4*(size/4+1);
  return (float*)_aligned_malloc(size*sizeof(float),16);
}

Filter :: Filter(int nmax)
{
  this->nmax = nmax;
  b = (float*)malloc(2*(nmax+1)*sizeof(float));	
  x = (float*)malloc(2*(nmax+1)*sizeof(float));
  memset(x,0,2*(nmax+1)*sizeof(float));
  memset(b,0,2*(nmax+1)*sizeof(float));
  xc = x;
  xend = x+2*nmax;
}

void Filter :: init()
{
  bend = b + (n << 1);
}

float Db(float B, float f, int M) 
{
  float C1,C2,k1,k2,k3;
  if(M==4) {
    C1 = .069618;
    C2 = 2.0427;
    k1 = -.00050469;
    k2 = -.0064264;
    k3 = -2.8743;
  } else {
    C1 = .071089;
    C2 = 2.1074;
    k1 = -.0026580;
    k2 = -.014811;
    k3 = -2.9018;
  }

  float logB = log(B);
  float kd = exp(k1*logB*logB + k2*logB + k3);
  float Cd = exp(C1*logB+C2);
  float halfstep = pow(2.0,1.0/12.0);
  float Ikey = log(f*halfstep/27.5) / log(halfstep);
  float D = exp(Cd - Ikey*kd);
  return D;
}


void complex_divide(float Hn[2], float Hd[2], float H[2]) 
{
	float d2 = Hd[0]*Hd[0] + Hd[1]*Hd[1];
	H[0] = (Hn[0]*Hd[0] + Hn[1]*Hd[1])/d2;
  H[1] = (Hn[1]*Hd[0] - Hn[0]*Hd[1])/d2;
}


Filter :: ~Filter()
{
  _aligned_free(b);
  _aligned_free(x);
}

float Filter :: phasedelay(float omega) 
{
  float Hn[2];
  float Hd[2];
  float H[2];

  Hn[0] = 0.0; Hn[1] = 0.0;
  Hd[0] = 0.0; Hd[1] = 0.0;

   for(int k=0;k<=n;k++) {
    int k2 = (k<<1);
    float c = cos(k * omega);
    float s = sin(k * omega);
    Hn[0] += c * b[k2];
    Hn[1] += s * b[k2];
    Hd[0] += c * b[k2+1];
    Hd[1] += s * b[k2+1];
  }
  complex_divide(Hn,Hd,H);
  float arg = atan2(H[1],H[0]);
  if(arg<0) arg = arg + 2*PI;
  
  return arg/omega;
}


float Filter :: groupdelay(float omega)
{
  float dw = .001;
  float omega2 = omega + dw;
  float omega1 = omega - dw;
  return (omega2*phasedelay(omega2) - omega1*phasedelay(omega1))/(omega2-omega1);
}


void Filter :: merge(const Filter &c1, const Filter &c2)
{
	int n = c1.n + c2.n;
	this->n = n;
	for(int j=0;j<=n;j++) {
    int j2 = (j<<1);
		b[j2] = 0;
    b[j2+1] = 0;
	}
	for(int j=0;j<=c1.n;j++) {
    int j2 = (j<<1);
		for(int k=0;k<=c2.n;k++) {
      int k2 = (k<<1);
			b[j2+k2] += c1.b[j2]*c2.b[k2];
			b[j2+k2+1] += c1.b[j2+1]*c2.b[k2+1];
		}
	}
	
	init();
}

void Thiran :: create(float D, int N) 
{
  if(N < 1) {
    n = 0;
    b[0] = 1;
    b[1] = 0;
    init();
    return;
  }

  
  int choose = 1;
  for(int k=0;k<=N;k++) {
    float ak = choose;
    for(int n=0;n<=N;n++) {
      ak *= ((float)D-(float)(N-n));
      ak /= ((float)D-(float)(N-k-n));
    }
    b[(k<<1)+1] = ak;
    b[(N-k)<<1] = ak;
    choose = (-choose * (N-k)) / (k+1); 
  }  

  n = N;
  init();
}

void ThiranDispersion :: create(float B, float f, int M)
{
  int N = 2;
  float D;
  D = Db(B,f,M);

  if(D<=1.0) {
    n = 2;	
    b[1] = 1;
    b[3] = 0;
    b[5] = 0;
    b[0] = 1;
    b[2] = 0;
    b[4] = 0;	
    init();
  } else {
    Thiran :: create(D,N);
  }
}

void BiquadHP :: create(float omega, float Q)
{
  float a = 1.0/(2.0*tan(0.5*omega));
  float a2 = a*a;
  float aoQ = a/Q;
  float d = (4*a2+2*aoQ+1);

  b[1] = 1.0;
  b[3] = -(8*a2-2) / d;
  b[5] = (4*a2 - 2*aoQ + 1) / d;

  b[0] = 4*a2/d;
  b[2] = -8*a2/d;
  b[4] = 4*a2/d;

  n = 2;
  init();
}

DWGResonator :: DWGResonator() 
{
  x1 = 0;
  x2 = 0;
}

float DWGResonator :: go(float in)
{
  float x1t = g * x1;
  float v = c * (x1t + x2);
  x1 = v - x2 + b1t * in;
  x2 = x1t + v;
  
  //return in;
  return x2;
}

void DWGResonator :: create(float omega, float gamma)
{
  g = exp(-2*gamma);
  c = sqrt(1/(1+(square(tan(omega) * (1+g)) + square(1-g)) / (4*g)));
  if(omega > HALFPI) c = -c;    
  b1t = sqrt((1-c)/(1+c));
  //x1 = 0;
  //x2 = 0;
}


ConvolutionResonator :: ConvolutionResonator()
{
  b = NULL;
  x = NULL;
}

ConvolutionResonator :: ~ConvolutionResonator()
{
  if(b) delete [] b;
  if(x) delete [] x;
}

void ConvolutionResonator :: create(float omega, float gamma)
{ 
  static const float logeps = 5;
  size = lrintf(logeps / gamma);
  if(b) delete [] b;
  if(x) delete [] x;
  b = new float[size];
  x = new float[size];
  xc = x;
  xend = x + size;
  bend = b + size;
  float scale = exp(-gamma);
  float e = 1.0f;
  float ph = 0.0f;
  for(int k=0; k<size; k++) {
    b[k] = e * sin(ph);
    ph += omega;
    while(ph>TWOPI) ph -= TWOPI;
    e *= scale;
  }
}

float ConvolutionResonator :: go(float in)
{
  float *b = this->b;
  float *x = this->xc;
  float out = *(b) * in;  
  b++;
  x++;

  while(b <= bend) {
    if(x>xend) x = this->x;
    out += *(b) * *(x);
    b++;
    x++;
  }
  x = this->xc;
  *(x) = in;
  x--; if(x<this->x) x = xend; this->xc = x;

  return out;
}

void MSDFilter :: create(float Fs, float m, float k, float mu, float RT) 
{
  float alpha = 2 * Fs;
  float a2m = alpha * alpha * m;
  float a0 = k + a2m + alpha*(mu+RT);

  b[0] = alpha / a0;
  b[2] = 0;
  b[4] = -b[0];
  b[1] = 1;
  b[3] = (2*k-2*a2m)/a0;
  b[5] = (k + a2m - alpha*(mu+RT)) / a0;

  n = 2;
  init();
}

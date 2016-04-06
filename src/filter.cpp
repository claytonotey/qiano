#include "filter.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "AudioFFT.h"
#include "filter.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

inline void complex_divide(float Hn[2], float Hd[2], float H[2]) 
{
	float d2 = Hd[0]*Hd[0] + Hd[1]*Hd[1];
	H[0] = (Hn[0]*Hd[0] + Hn[1]*Hd[1])/d2;
  H[1] = (Hn[1]*Hd[0] - Hn[0]*Hd[1])/d2;
}

float Filter :: phasedelay(float omega) 
{
  float Hn[2];
  float Hd[2];
  float H[2];

  Hn[0] = 0.0; Hn[1] = 0.0;
  Hd[0] = 0.0; Hd[1] = 0.0;

  for(int k=0;k<=n;k++) {
    int j = n - k;
    float c = cos(k * omega);
    float s = sin(k * omega);
    Hn[0] += c * b[j];
    Hn[1] += s * b[j];
    Hd[0] -= c * a[j];
    Hd[1] -= s * a[j];
  }
  complex_divide(Hn,Hd,H);
  float arg = atan2(H[1],H[0]);
  if(arg<0) arg = arg + 2*PI;
  
  return arg/omega;
}

void Filter :: merge(const Filter &f)
{
  int n1 = n;
  n = n1 + f.n;
  
  float *aa = (float*)malloc((n1+1)*sizeof(float));
  float *bb = (float*)malloc((n1+1)*sizeof(float));
  memcpy(aa,a,(n1+1)*sizeof(float));
  memcpy(bb,b,(n1+1)*sizeof(float));
  memset(a,0,(n+1)*sizeof(float));
  memset(b,0,(n+1)*sizeof(float));

	for(int j=0;j<=n1;j++) {
		for(int k=0;k<=f.n;k++) {
			b[j+k] += bb[n1-j]*f.b[f.n-k];
			a[j+k] -= aa[n1-j]*f.a[f.n-k];
		}
	}
	
  free(aa);
  free(bb);
	init(upsample);
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

Filter :: ~Filter()
{

}

Filter :: Filter(int nmax)
{
  this->nmax = nmax;
  int n4 = nmax / 4;
  posix_memalign((void**)&b,32,(n4+3)*sizeof(vec4));
  posix_memalign((void**)&a,32,(n4+3)*sizeof(vec4));
  posix_memalign((void**)&x,32,2*(n4+3)*sizeof(vec4));
  posix_memalign((void**)&y,32,2*(n4+3)*sizeof(vec4));

  memset(x,0,2*(n4+3)*sizeof(vec4));
  memset(y,0,2*(n4+3)*sizeof(vec4));
  memset(b,0,(n4+3)*sizeof(vec4));
  memset(a,0,(n4+3)*sizeof(vec4));

  int xsize = 2*4*(n4 + 2);
  xend = x + xsize;
  yend = y + xsize;
  xc = x;
  yc = y;

  xskip = xsize;
}

void Filter :: init(int upsample)
{
  this->upsample = upsample;
  bend = b + n;
  aend = a + n - 1;
  aend4 = a + n - 4;

  vec4 v0 = {a[3], a[3], a[3], 0};
  a0 = v0;

  vec4 v1 = {a[2], a[2], 0, a[1]};
  a1 = v1;
  
  vec4 v2 = {a[1], 0, a[1], a[2] + a[1] * a[1]};
  a2 = v2;

  vec4 v3 = {0, a[1], a[1] * a[1] + a[2], a[1] * (a[1] * a[1] + 2 * a[2]) + a[3]};
  a3 = v3;

  // reverse order
  for(int i=0; i<=n/2; i++) {
    float tmp = a[i];
    a[i] = a[n-i];
    a[n-i] = tmp;
    
    tmp = b[i];
    b[i] = b[n-i];
    b[n-i] = tmp;
  }

  a[n+1] = 0;
  b[n+1] = 0;

}

float Filter :: filter(float in)
{
  float *b = this->b;
  float *x = this->xc - n;
  if(x < this->x) x += xskip;

  *xc = in;
  xc++; 
  if(xc>=xend) xc = this->x;


  float out = 0;
  while(b <= bend) {
    if(x >= xend) x -= xskip;
    out += *b * *x;
    b += upsample;
    x += upsample;
  }

  float *a = this->a;
  float *y = this->yc - n;
  if(y < this->y) y += xskip;
  while(a <= aend) {
    if(y >= yend) y -= xskip;
    out += *a * *y;
    a += upsample;
    y += upsample;
  }

  *yc = out;
  yc++;
  if(yc>=yend) yc = this->y;

  return out;
}

vec4 Filter :: filter4(vec4 in)
{
  float *b = this->b;
  float *x = this->xc - n;
  if(x < this->x) x += xskip;

  _mm_store_ps(xc, in);
  if(xc == this->x) _mm_store_ps(xend, in);
  xc+=4; 
  if(xc>=xend) xc = this->x;
  
  vec4 out = {0};

  while(b <= bend) {
    if(x >= xend) x = this->x;
    out = _mm_fmadd_ps(_mm_broadcast_ss(b), _mm_loadu_ps(x), out);
    b+=upsample;
    x+=upsample;
  }
  
  vec4 outa = {0};
  float *a = this->a;
  float *y = this->yc - n;
  if(y < this->y) y += xskip;
  
  while(a <= aend4) {
    if(y >= yend) y = this->y;
    outa = _mm_fmadd_ps(_mm_broadcast_ss(a), _mm_loadu_ps(y), outa);
    a+=upsample;
    y+=upsample;
  }
  
  out += outa;
  
  y = yc - 4;
  if(y < this->y) y += xskip;
  vec4 y4 = _mm_loadu_ps(y);  
  
  out = _mm_fmadd_ps(_mm_shuffle_ps(y4, y4, _MM_SHUFFLE(0,3,2,1)), a0, out);
  out = _mm_fmadd_ps(_mm_shuffle_ps(y4, out, _MM_SHUFFLE(2,0,3,2)), a1, out);
  out = _mm_fmadd_ps(_mm_shuffle_ps(y4, out, _MM_SHUFFLE(1,1,0,3)), a2, out);
  out = _mm_fmadd_ps(_mm_shuffle_ps(out, out, _MM_SHUFFLE(0,0,0,0)), a3, out);
  
  _mm_store_ps(yc, out);
  if(yc == this->y) _mm_store_ps(yend, out);
  yc+=4; 
  if(yc>=yend) yc = this->y;
    
  return out;
}



void Thiran :: create(float D, int N, int upsample)
{
  if(N < 1) {
    n = 0;
    a[0] = -1;
    b[0] = 1;
    init();
    return;
  }
  n = N*upsample;
  memset(a,0,(n+1)*sizeof(float));
  memset(b,0,(n+1)*sizeof(float));

  int choose = 1;
  for(int k=0;k<=N;k++) {
    float ak = choose;
    for(int n=0;n<=N;n++) {
      ak *= ((float)D-(float)(N-n));
      ak /= ((float)D-(float)(N-k-n));
    }
    a[upsample*(k)] = -ak;
    b[upsample*(N-k)] = ak;
    choose = (-choose * (N-k)) / (k+1); 
  }  

  init(upsample);
}

void ThiranDispersion :: create(float B, float f, int M, int downsample, int upsample)
{
  int N = 2;
  float D; 
  D = Db(B,f,M);
  D /= downsample;

  if(D<=1.0) {
    n = 2*upsample;	
    a[0] = -1;
    a[upsample] = 0;
    a[2*upsample] = 0;
    b[0] = 1;
    b[upsample] = 0;
    b[2*upsample] = 0;	
    init(upsample);
  } else {
    Thiran :: create(D,N,upsample);
  }
}

void BiquadHP :: create(float omega, float Q)
{
  float A = 1.0/(2.0*tan(0.5*omega));
  float A2 = A*A;
  float AoQ = A/Q;
  float d = (4*A2+2*AoQ+1);

  a[0] = -1.0;
  a[1] = (8*A2-2) / d;
  a[2] = -(4*A2 - 2*AoQ + 1) / d;

  b[0] = 4*A2/d;
  b[1] = -8*A2/d;
  b[2] = 4*A2/d;

  n = 2;
  init();
}

DWGResonator :: DWGResonator() 
{
  x1 = 0;
  x2 = 0;
}

vec4 DWGResonator :: go4(vec4 vin)
{
  float out[4] __attribute__((aligned(32)));
  float in[4] __attribute__((aligned(32)));
  _mm_store_ps(in, vin);

  for(int i=0; i<4; i++) {
    float x1t = g * x1;
    float v = c * (x1t + x2);
    x1 = v - x2 + b1t * in[i];
    x2 = x1t + v;
    out[i] = x2;
  }
   
  return _mm_load_ps(out);
}

float DWGResonator :: go(float in)
{
  float x1t = g * x1;
  float v = c * (x1t + x2);
  x1 = v - x2 + b1t * in;
  x2 = x1t + v;
  
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


void MSD2Filter :: filter(float in[2], float out[2])
{
  out[0] = f11 * in[0] + f12 * in[1];
  out[1] = f21 * in[0] + f22 * in[1];
}

void MSD2Filter :: filter4(vec4 in[2], vec4 out[2])
{
  out[0] = _mm_broadcast_ss(&f11) * in[0] + _mm_broadcast_ss(&f12) * in[1];
  out[1] = _mm_broadcast_ss(&f21) * in[0] + _mm_broadcast_ss(&f22) * in[1];
}

void MSD2Filter :: create(float Fs,
                          float m1, float k1, float R1, 
                          float m2, float k2, float R2,
                          float R12, float k12, float Zn, float Z) 
{
  float det = (R1 + Zn) * (R2 + Zn) - R12 * R12;
  f11 = 2.0 * Z * (R2 + Zn) / det;
  f12 = -2.0 * Z * R12 / det;
  f21 = -2.0 * Z * R12 / det;
  f22 = 2.0 * Z * (R1 + Zn) / det;
}

vec8 ResampleFIR :: filter8(vec8 in)
{
  float *b = this->b;
  float *x = this->xc - (ResampleFilterSize - 1);
  if(x < this->x) x += xsize;

  _mm256_store_ps(xc, in);
  if(xc == this->x) _mm256_store_ps(xend, in);
  xc+=8; 
  if(xc>=xend) xc = this->x;
  
  vec8 out = {0};
 
  while(b <= bend) {
    if(x >= xend) x = this->x;
    out = _mm256_fmadd_ps(_mm256_broadcast_ss(b), _mm256_loadu_ps(x), out);
    b++;
    x++;
  }

  return out;
}

ResampleFIR :: ResampleFIR()
{
  xsize = ResampleFilterSize * 4;
  memset(x,0,(xsize+16)*sizeof(float));
  xc = x;

  bend = this->b + ResampleFilterSize - 1;
  xend = this->x + xsize;
  bInit = false;
}

int ResampleFIR :: getDelay()
{
  return ResampleFilterSize / 2;
}

bool ResampleFIR :: isCreated()
{
  return bInit;
}

void ResampleFIR :: create(int resample)
{
  bInit = true;
  float f[ResampleFilterSize];
  float im[ResampleFilterSize];
  
  memset(f,0,ResampleFilterSize*sizeof(float));
  memset(im,0,ResampleFilterSize*sizeof(float));

  for(int i=0; i<ResampleFilterSize/resample/2; i++) {
    f[i] = 1.0;
    f[ResampleFilterSize - i - 1] = 1.0;
  }
  audiofft::AudioFFT fft;
  fft.init(ResampleFilterSize);
  fft.ifft(b, f, im);
  // fftshift
  int shift = ResampleFilterSize / 2;
  for(int i=0; i<shift; i++) {
    float tmp  = b[i];
    b[i] = b[i+shift];
    b[i+shift] = tmp;
  }

  for(int i=0; i<shift; i++) {
    float tmp  = b[i];
    b[i] = b[ResampleFilterSize-1-i];
    b[ResampleFilterSize-1-i] = tmp;
  }

  
}

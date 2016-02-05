#include "dwgs.h"
#include "types.h"
#include <math.h>
#include <stdio.h>
#include "utils.h"
#include <stdlib.h>


/*
  For transverse waves:
 There is a total delay of D+d.  D is an integer, d is fractional part.  
 Store an array of N=2*ceil((D+d)/2)
 The delay is split between elft and right going sections.
 The delay lines, and filters will access this at different points.
 It has an intepretation in x-space.  If the total filter is H(z), then 
 H(z)^(1/N)

For longitudinal waves, store another array of size N.
sTHere is just a leftgoing and righthoing delay line.
Every time step, transverse waves force every point in left, and rightgoing delay lines,
The l-waves move faster, so step the delay line multiple times
It might makes since to limit this to an integer delay
 */

  /*
     
 pin - 0-|--1 -- 2 bridge - soundboard
       | 
       3 
     hammer

  The hammer is imparted with an initial velocity.  
  The hammer velocity becomes an input load using the impedances at the junction.
  The 0 & 1 strings carry the velocity wave.
  The 0 string is terminated at an infinite impedance pin.
  The 1 string interacts with the soundboard by inducing a load
  There are no additional ingoing waves into the strings from either the hammer or the bridge

*/

//#define STRING_DEBUG 1

dwgs :: dwgs(float Fs) :
  fracDelayTop(8), fracDelayBottom(8), hammerDelay(8)

{
	this->Fs = Fs;
  memset(modeTable,0,nMaxLongModes*sizeof(float*));

  a0_1 = 0.0f;
  a0_2 = 0.0f;
  a0_3 = 0.0f;
  a0_4 = 0.0f;
  a1_1 = 0.0f;
  a1_2 = 0.0f;
  a1_3 = 0.0f;
  a1_4 = 0.0f;
  a0_5 = 0.0f;
}

void dwgs :: set(int upsample, float f, float c1, float c3, float B, float L, float longFreq1, float gammaL, float gammaL2, float inpos, float Z) 
{
  this->upsample = upsample;

  this->L = L;
  this->f = f;
  this->inpos = inpos;
  
	if(f > 400) {
		M = 1;		
	} else {
		M = 4;
	}

  for(int m=0;m<M;m++) {		
    dispersion[m].create(B,f,M,upsample);
  }
  this->B = B;
  this->longFreq1 = longFreq1;


  float deltot = Fs/f*upsample; 
  this->omega = TWOPI / deltot;
  float dispersiondelay = M*dispersion[0].phasedelay(omega);  

  del0 = (int)(0.5 * (inpos*deltot));
  del1 = (int)((inpos*deltot) - 6.0);
  if(del1 < 1) abort();
  float dHammer = inpos * deltot - del1 - 1;
  hammerDelay.create(dHammer,(int)dHammer);
  del5 = (int)dHammer;

  // XXX resize arrays
  del3 = (int)(0.5*(deltot - inpos * deltot) - dispersiondelay - 6.0);
  del2 = (int)(0.5*(deltot - inpos * deltot) - 5.0);

  if(del2 < 1) abort();
  if(del3 < 1) abort();

  float delHalf = 0.5f * deltot;
 
  d0.setDelay(1);
  d1.setDelay(del1-1);
  d2.setDelay(del2-1);
  d3.setDelay(del3-1);

  float delHammerHalf = 0.5 * inpos * deltot;

  float dTop = delHalf - delHammerHalf - del2;

  
  del4 = (int)dTop;
  
  fracDelayTop.create(dTop,(int)dTop);

  dBottomAndLoss = delHalf - delHammerHalf - del3 - dispersiondelay;

  delTab = del0 + del2 + del4;
  float delta = dTop - (int)dTop;
  float delta2 = delHalf - delTab - delta;
  
  posix_memalign((void**)&wave,32,delTab*sizeof(float));
  posix_memalign((void**)&wave0,32,delTab*sizeof(float));
  posix_memalign((void**)&wave1,32,delTab*sizeof(float));
  posix_memalign((void**)&wave2,32,delTab*sizeof(float));
  posix_memalign((void**)&dwave,32,delTab*sizeof(float));
  posix_memalign((void**)&Fl,32,delTab*sizeof(float));
  posix_memalign((void**)&dFl,32,delTab*sizeof(float));
  memset(wave,0,delTab*sizeof(float));
  memset(wave0,0,delTab*sizeof(float));
  memset(wave1,0,delTab*sizeof(float));
  memset(wave2,0,delTab*sizeof(float));
  
  //fprintf(stderr,"dwgs top %d %d %d %d %g %g %g %g %g %g %g %g\n",del0,del1,del2,del3,dHammer+1+del2+dTop,del1+del3+dispersiondelay+lowpassdelay+dBottom, dTop, dBottom, dHammer, inpos*deltot, lowpassdelay, dispersiondelay);
  
  nLongModes = (int)(0.5f * Fs / longFreq1 - 0.5);
  //fprintf(stderr,"long modes %d \n",nLongModes);
  if(nLongModes >= nMaxLongModes) abort();
  
  for(int k=1; k <= nLongModes; k++) {
    float omegak = TWOPI * k * longFreq1 / Fs / upsample;
    float gammak = gammaL * (1.0 + gammaL2 * k * k);
    fLong[k] = omegak / TWOPI;
    longModeResonator[k].create(omegak, gammak);
    
    //fprintf(stderr,"fL%d=%g L%d/fT = %g\n",k,k*longFreq1,k,k*longFreq1/f );
    if(delTab) {
      float n = PI * k / delHalf;
      if(modeTable[k]) delete modeTable[k];
      posix_memalign((void**)&modeTable[k],32,delTab*sizeof(float));
      for(int i=0; i<delTab; i++) {
        float d = i + delta;
        modeTable[k][i] = sin(d*n);
      }
      modeTable[k][0] *= (0.5 + delta);
      modeTable[k][delTab-1] *= (0.5 + delta2);
    }
  }


  damper(c1,c3,gammaL,gammaL2);
}

void dwgs :: damper(float c1, float c3, float gammaL, float gammaL2) 
{
	loss.create(f,c1,c3,upsample);
  float lowpassdelay = loss.phasedelay(omega);
  float dBottom = dBottomAndLoss - lowpassdelay;
  fracDelayBottom.create(dBottom,5);
}

dwgs :: ~dwgs()
{
}

int dwgs :: getMaxDecimation(float f, int upsample, float magic) 
{
  //float magic = 120;
  int decimation = 1;
  do {
    float deltot = Fs*upsample / f;
    //printf("%g %g\n",deltot, decimation * magic);
    if(deltot > decimation * magic) {
      decimation *= 2;
    } else {
      break;
    }
  } while(true);

  return decimation;
}

int dwgs :: getMinUpsample(float f, float inpos)
{
  int upsample = 1; 
  
  do {
    float deltot = upsample * Fs/f; 
    del1 = (int)((inpos*deltot) - 6.0);
    if(del1 < 4) {
      upsample *= 2;
    } else {
      break;
    }
  } while(true);
  
  return upsample;
}

float dwgs :: input_velocity() {
  return a0_3 + a1_2;
}

float dwgs :: next_input_velocity() {
  return d2.probe() + d1.probe();
}

// returns input to soundboard
float dwgs :: go_string()
{
  float a;
  a = d0.goDelay(a0_2);
  a = hammerDelay.filter(a);
  a1_2 = d1.goDelay(-a);
  a0_3 = d2.goDelay(a0_4);
  a = d3.goDelay(a1_3);
  for(int m=0;m<M;m++) {
    a = dispersion[m].filter(a);
  }
  a = loss.filter(a);
  a1_4 = fracDelayBottom.filter(a);
  
  return a1_4;
}

float dwgs :: longTran(dwgs *top, dwgs *bottom) {
  return square(a0_5 - a1_4);
}


inline float hermite(float beta0, float beta3, float m0, float m1, float t) {
  if(t == 0.0) {
    return beta0;
  } else if(t == 1.0) {
    return beta3;
  }

  float beta1 = beta0 + 0.33333333333333333333 * m0;
  float beta2 = beta3 - 0.33333333333333333333 * m1;
  
  float t1 = 1.0 - t;
  beta0 = t1 * beta0 + t * beta1;
  beta1 = t1 * beta1 + t * beta2;
  beta2 = t1 * beta2 + t * beta3;

  beta0 = t1 * beta0 + t * beta1;
  beta1 = t1 * beta1 + t * beta2;

  beta0 = t1 * beta0 + t * beta1;
  
  return beta0;
}

//0 2
//1 3
/* wave is the difference of top and bottom strings 
   wave[0] is on the right side, so samples are delayed by frac(topDelay)
*/
// d/dx (dy/dx)^2
// dy/dx = dy/dt / dx/dt = v / vTran

void dwgs :: tran2long(dwgs *top, dwgs *bottom, int offset, int delay)
{
  int del2 = top->del2;
  int del0 = top->del0;
    
  float *x = top->d2.x;
  int cur = (top->d2.cursor + DelaySize - delay + del4) % DelaySize;

  int n = del2 + top->del4;
  if(n <= cur) {
    float *x1 = x + cur;
    for(int i=0; i<n; i++) {
      wave1[i] = x1[-i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=0; i<=cur; i++) {
      wave1[i] = x1[-i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=cur+1; i<n; i++) {
      wave1[i] = x1[-i];
    }
  }

  x = top->d0.x;
  cur = (top->d0.cursor + DelaySize - delay) % DelaySize;
  float *wave10 = wave1 + del2 + del4;

  if(del0 <= cur) {
    float *x1 = x + cur;
    for(int i=0; i<del0; i++) {
      wave10[i] = x1[-i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=0; i<=cur; i++) {
      wave10[i] = x1[-i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=cur+1; i<del0; i++) {
      wave10[i] = x1[-i];
    }
  }

#ifdef STRING_DEBUG
  for(int i=0; i<delTab; i++) {
    printf("%g ",wave1[i]);
  }
  printf("\n");
#endif


  /********* bottom *********/
  
  x = bottom->d1.x;
  cur = (bottom->d1.cursor + DelaySize - delay + del5 + 1 - del0) % DelaySize;
  wave10 = wave1 + delTab;
  if(del0 < cur) {
    float *x1 = x + cur;
    for(int i=-1; i>=-del0; i--) {
      wave10[i] -= x1[i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=-1; i>=-cur; i--) {
      wave10[i] -= x1[i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=-cur-1; i>=-del0; i--) {
      wave10[i] -= x1[i];
    }
  }

  x = bottom->d3.x;
  cur = (bottom->d3.cursor + DelaySize - delay) % DelaySize;
  n = del2 + top->del4;
  wave10 = wave1 + n;
  if(n < cur) {
    float *x1 = x + cur;
    for(int i=-1; i>=-n; i--) {
      wave10[i] -= x1[i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=-1; i>=-cur; i--) {
      wave10[i] -= x1[i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=-cur-1; i>=-n; i--) {
      wave10[i] -= x1[i];
    }
  }

#ifdef STRING_DEBUG
  for(int i=0; i<delTab; i++) {
    printf("%g ",wave1[i]);
  }
  printf("\n");
#endif

  for(int i=0; i<delTab; i++) {
    wave0[i] = wave2[i];
    wave2[i] = wave[i];
    dwave[i] = 0.5 * (wave1[i] - wave0[i]);
    wave[i] = wave1[i];
  }

  int last = delTab - 1;
  if(offset == 0) {
    Fl[0] = 2.0 * (wave[1] - wave[0]) * wave[0];
    for(int i=1; i<last; i++) {
      Fl[i] = (wave[i+1] - wave[i-1]) * wave[i];
    }
    Fl[last] = 2.0 * (wave[last] - wave[last-1]) * wave[last];

    for(int k=1; k<=nLongModes; k++) {
      float *tab = modeTable[k];
      hermite_p0[k] = hermite_p1[k];
      //hermite_p1[k] = sse_dot(delTab,tab,Fl);
      hermite_p1[k] = dsp_dot(delTab,tab,Fl);
    }
  } else if(offset == 1) {
    dFl[0] = 2.0 * ((dwave[1] - dwave[0]) * wave[0] + (wave[1] - wave[0]) * dwave[0]);
    for(int i=1; i<last; i++) {
      dFl[i] = (dwave[i+1] - dwave[i-1]) * wave[i] + (wave[i+1] - wave[i-1]) * dwave[i];
    }
    dFl[last] = 2.0 * ((dwave[last] - dwave[last-1]) * wave[last] + (wave[last] - wave[last-1]) * dwave[last]);
     
    for(int k=1; k<=nLongModes; k++) {
      float *tab = modeTable[k];
      hermite_m0[k] = hermite_m1[k];
      //hermite_m1[k] = sse_dot(delTab,tab,dFl);      
      hermite_m1[k] = dsp_dot(delTab,tab,dFl);      
    }
  }
}

void dwgs :: longForce(float *Ft, int n, int decimation) {
  float t = 0;
  float dt = 1.0 / decimation;
  float s = 1.0;
  for(int i=0; i<n; i++) {
    float FbL = 0;
    for(int k=1; k<=nLongModes; k++) {
      float F = hermite(hermite_p0[k], hermite_p1[k], decimation*hermite_m0[k], decimation*hermite_m1[k] , t);
      float Fl = longModeResonator[k].go(F);
      FbL += Fl;
      //printf("%g %g %g %g %g %g ", F, Fl, hermite_p0[k], hermite_p1[k] , hermite_m0[k], hermite_m1[k]);

    }
    //printf("\n");
    Ft[i] = FbL;
    t += dt;
  }
  for(int k=1; k<=nLongModes; k++) {
    //printf("%g %g %g %g ", hermite_p0[k], hermite_p1[k] , hermite_m0[k], hermite_m1[k]);
  }
  //printf("\n");
}

float dwgs :: go_soundboard(float load_h, float load_sb) 
{
  a0_2 = a0_3 + load_h;
  a1_3 = a1_2 + load_h;
  a0_5 = load_sb - a1_4;
  a0_4 = fracDelayTop.filter(a0_5);
  
  return load_sb - 2 * a1_4;
}


/*

| a0_1 --- del0 ---                 a0_2 | a0_3 --- del2 --- a0_4 --- fracDelayTop                            |   
|                                        H                                                                    | 
| a1_1 --- hammerDelay --- del1 --- a1_2 | a1_3  ... del3 ... dispersion + loss + fracDelayBottom --- a1_4 |
 */

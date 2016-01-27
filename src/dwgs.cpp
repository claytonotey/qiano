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

dwgs :: dwgs(float Fs, int type) :
  fracdelay(8)

{
	this->Fs = Fs;
  this->type = type;

  decimation = 8;
  upsample = 0;
  memset(modeTable,0,nMaxLongModes*sizeof(float*));

  a0_1 = 0.0f;
  a0_2 = 0.0f;
  a0_3 = 0.0f;
  a0_4 = 0.0f;
  a1_1 = 0.0f;
  a1_2 = 0.0f;
  a1_3 = 0.0f;
  a1_4 = 0.0f;
  
}

void dwgs :: set(int upsample, float f, float c1, float c3, float B, float L, float longFreq1, float gammaL, float gammaL2, float inpos, float Z, int hammerNutDelay) 
{
  if(this->upsample != upsample) {
    switch(upsample) {
    case 1:
      loss = new Loss<1>();
      break;
    case 2:
      loss = new Loss<2>();
      break;
    case 3:
      loss = new Loss<3>();
      break;
    case 4:
      loss = new Loss<4>();
      break;
    }
  } 
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
    dispersion[m].create(B,f,M);
  }
  this->B = B;
  this->longFreq1 = longFreq1;

  float deltot = Fs/f;
  if(hammerNutDelay > 0) {
    del1 = hammerNutDelay;
  } else {
    del1 = (int)(0.5f*inpos*deltot);
    if(del1 < 2) {
      del1 = 1;
    }
  }
  damper(c1,c3,gammaL,gammaL2);
  
  t = 0;
}

int dwgs :: getHammerNutDelay()
{
  return del1;
}

void dwgs :: damper(float c1, float c3, float gammaL, float gammaL2) 
{

	loss->create(f,Fs,c1,c3);  
  float deltot = Fs/f; 
  float omega = TWOPI * f / Fs;
  float lowpassdelay = loss->phasedelay(omega);
  float dispersiondelay = M*dispersion[0].phasedelay(omega);  

  // XXX resize arrays
  if(type == NoTopDispersion) {
    del3 = (int)(0.5*(deltot-del1 - del1)-dispersiondelay-lowpassdelay-5.0);
    del2 = (int)(0.5*(deltot-del1 - del1));
    int del12 = del1 + del2;
    posix_memalign((void**)&wave,32,(del12+1)*sizeof(float));
    posix_memalign((void**)&wave1,32,(del12+1)*sizeof(float));
    posix_memalign((void**)&dwave,32,(del12+1)*sizeof(float));
    posix_memalign((void**)&Fl,32,(del12+1)*sizeof(float));
    posix_memalign((void**)&dFl,32,(del12+1)*sizeof(float));
    memset(wave,0,(del12+1)*sizeof(float));
  } else if(type == NoBottomDispersion) {
    del3 = (int)(0.5*(deltot-del1 - del1));
    del2 = (int)(0.5*(deltot-del1 - del1)-dispersiondelay-lowpassdelay-5.0);
    int del13 = del1 + del3;
    posix_memalign((void**)&wave,32,(del13+1)*sizeof(float));
    posix_memalign((void**)&wave1,32,(del13+1)*sizeof(float));
    posix_memalign((void**)&dwave,32,(del13+1)*sizeof(float));
    posix_memalign((void**)&Fl,32,(del13+1)*sizeof(float));
    posix_memalign((void**)&dFl,32,(del13+1)*sizeof(float));
    memset(wave,0,(del13+1)*sizeof(float));
  } else {
    del3 = (int)(0.5*(deltot-del1 - del1)-dispersiondelay);
    del2 = (int)(0.5*(deltot-del1 - del1)-lowpassdelay-5.0);
  }

  if(del2 < 2)
    del2 = 1;
  if(del3 < 2)
    del3 = 1;

  float D = (deltot-(float)(del1 + del1+del2+del3)-dispersiondelay-lowpassdelay);
  if(D<1) D=1;


  fracdelay.create(D,(int)D);

  /*
  if(type == NoTopDispersion) {
    fprintf(stderr,"dwgs top %d %d %d %d %g %g %g\n",del1,del2,del3,del1+del2, del1+del3+dispersiondelay+lowpassdelay+D,lowpassdelay,dispersiondelay);
  } else if(type == NoBottomDispersion) {
    fprintf(stderr,"dwgs bot %d %d %d %d %g %g %g\n",del1,del2,del3,del1+del3, del1+del2+dispersiondelay+lowpassdelay+D,lowpassdelay,dispersiondelay);
  } else {
    fprintf(stderr,"dwgs meh %d %d %d %g %g %g %g\n",del1,del2,del3,del1+del2+lowpassdelay+D, del1+del3+dispersiondelay,lowpassdelay,dispersiondelay);
  }
  */

  d0.setDelay(del1-1);
  d1.setDelay(del1-1);
  d2.setDelay(del2-1);
  d3.setDelay(del3-1);

  int delTab;
  if(type == NoTopDispersion) {
    delTab = del1 + del2;
  } else if(type == NoBottomDispersion) {
    delTab = del1 + del3;
  } else {
    delTab = 0;
  }

  float delHalf = 0.5f * deltot;
  float delta = delTab - delHalf;
  iMax = delTab + 1;

  nLongModes = (int)(0.5f * Fs / longFreq1 - 0.5);
  if(nLongModes >= nMaxLongModes) abort();

  for(int k=1; k <= nLongModes; k++) {
    float omegak = TWOPI * k * longFreq1 / Fs;
    float gammak = gammaL * (1.0 + gammaL2 * k * k);
    fLong[k] = omegak / TWOPI;
    longModeResonator[k].create(omegak, gammak);

    //fprintf(stderr,"fL%d=%g L%d/fT = %g\n",k,k*longFreq1,k,k*longFreq1/f );
    if(delTab) {
      float n = PI * k / delHalf;
      if(modeTable[k]) delete modeTable[k];
      posix_memalign((void**)&modeTable[k],32,(delTab+1)*sizeof(float));
      for(int i=0; i<iMax; i++) {
        float d = i - delta;
        if(d < 0 ) d = 0;
        modeTable[k][i] = sin(d*n);
      }
    }
  }


  //fprintf(stderr,"%d %d\n",nTranModes,nLongModes);
}

dwgs :: ~dwgs()
{
}

float dwgs :: input_velocity() {
  return a0_3 + a1_2;
}


//0 2
//1 3

// returns input to soundboard
float dwgs :: go_string()
{
  a0_1 = d0.goDelay(a0_2);
  a1_2 = d1.goDelay(a1_1);
  a0_3 = d2.goDelay(a0_4);
  a1_4 = d3.goDelay(a1_3);
  return a1_4;
}

float dwgs :: longTran(dwgs *top, dwgs *bottom) {
  return square(top->a0_4 - bottom->a1_4);
}

// d/dx (dy/dx)^2
// dy/dx = dy/dt / dx/dt = v / vTran

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

void dwgs :: tran2long(dwgs *top, dwgs *bottom, bool bD)
{
  int del2 = top->del2;
  int del3 = bottom->del3;
  int del12 = top->del1 + top->del2;
  int del13 = bottom->del1 + bottom->del3;
  
  wave1[0] = top->a0_4;
  float *x = top->d2.x;
  int cur = top->d2.cursor;

  if(del2 <= cur) {
    float *x1 = x + cur;
    for(int i=1; i<del2; i++) {
      wave1[i] = x1[-i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=1; i<=cur; i++) {
      wave1[i] = x1[-i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=cur+1; i<del2; i++) {
      wave1[i] = x1[-i];
    }
  }
  wave1[del2] = top->a0_2;

  x = top->d0.x;
  cur = top->d0.cursor;
  float *wave10 = wave1 + del2;
  if(del1 <= cur) {
    float *x1 = x + cur;
    for(int i=1; i<del1; i++) {
      wave10[i] = x1[-i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=1; i<=cur; i++) {
      wave10[i] = x1[-i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=cur+1; i<del1; i++) {
      wave10[i] = x1[-i];
    }
  }
  wave1[del12] = top->a0_1;

  x = bottom->d1.x;
  cur = bottom->d1.cursor;

  wave1[del13] -= bottom->a1_1;
  wave10 = wave1 + del13;
  if(del1 <= cur) {
    float *x1 = x + cur;
    for(int i=-1; i>-del1; i--) {
      wave10[i] -= x1[i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=-1; i>=-cur; i--) {
      wave10[i] -= x1[i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=-cur-1; i>-del1; i--) {
      wave10[i] -= x1[i];
    }
  }
  wave1[del3] -= bottom->a1_2;

  x = bottom->d3.x;
  cur = bottom->d3.cursor;
  wave10 = wave1 + del3;
  if(del3 <= cur) {
    float *x1 = x + cur;
    for(int i=-1; i>-del3; i--) {
      wave10[i] -= x1[i];
    }
  } else {
    float *x1 = x + cur;
    for(int i=-1; i>=-cur; i--) {
      wave10[i] -= x1[i];
    }
    int cur2 = cur + DelaySize;
    x1 = x + cur2;
    for(int i=-cur-1; i>-del3; i--) {
      wave10[i] -= x1[i];
    }
  }
  wave1[0] -= bottom->a1_4;


  for(int i=0; i<=del12; i++) {
    dwave[i] = wave1[i] - wave[i];
    wave[i] = wave1[i];
  }


  Fl[0] = 2.0 * (wave[1] - wave[0]) * wave[0];
  for(int i=1; i<del12; i++) {
    Fl[i] = (wave[i+1] - wave[i-1]) * wave[i];
  }
  Fl[del12] = 2.0 * (wave[del12] - wave[del12-1]) * wave[del12];
  
  dFl[0] = 2.0 * ((dwave[1] - dwave[0]) * wave[0] + (wave[1] - wave[0]) * dwave[0]);
  for(int i=1; i<del12; i++) {
    dFl[i] = (dwave[i+1] - dwave[i-1]) * wave[i] + (wave[i+1] - wave[i-1]) * dwave[i];
  }
  dFl[del12] = 2.0 * ((dwave[del12] - dwave[del12-1]) * wave[del12] + (wave[del12] - wave[del12-1]) * dwave[del12]);


  for(int k=1; k<=nLongModes; k++) {
    float *tab = modeTable[k];
    if(bD) {

    } else {
      hermite_m0[k] = hermite_m1[k];
      hermite_m1[k] = sse_dot(iMax,tab,dFl);
      hermite_p0[k] = hermite_p1[k];
      hermite_p1[k] = sse_dot(iMax,tab,Fl);
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
      float F = hermite(hermite_p0[k], hermite_p1[k], hermite_m0[k], hermite_m1[k] , t);
      FbL += longModeResonator[k].go(F);  
    }
    Ft[i] = FbL;
    t += dt;
  }
}

float dwgs :: go_soundboard(float load_h, float load_sb) 
{
  float a;

  a1_1 = -a0_1;
  a0_2 = a0_3 + load_h;

  // bottom -> right
  a = a1_2 + load_h;
  if(type == BothDispersion) {
    for(int m=0;m<M;m++) {
      a = dispersion[m].filter(a);
    }
  } else if(type == NoTopDispersion) {
    for(int m=0;m<M;m++) {
      a = dispersion[m].filter(a);
    }
    a = loss->filter(a);
    a = fracdelay.filter(a);
  }
  a1_3 = a;

  a = load_sb - a1_4;
  if(type == BothDispersion) {
    a = loss->filter(a);
    a = fracdelay.filter(a);
  } else if(type == NoBottomDispersion) {
    for(int m=0;m<M;m++) {
      a = dispersion[m].filter(a);
    } 
    a = loss->filter(a);
    a = fracdelay.filter(a);
  }
  a0_4 = a;
  
  return load_sb - 2 * a1_4;
  

}


/*
both
| a0_1 --- del0 --- a0_2 | a0_3 --- del2 --- a0_4 ...   merged ...   |   
|                        H                                           | 
| a1_1 --- del1 --- a1_2 | ... dispersion ... a1_3 --- del3 --- a1_4 |

notop
| a0_1 --- del0 --- a0_2 | a0_3 --- del2 --- a0_4                             |   
|                        H                                                    | 
| a1_1 --- del1 --- a1_2 | ... dispersion + merged ... a1_3 --- del3 --- a1_4 |

nobottom
| a0_1 --- del0 --- a0_2 | a0_3 --- del2 --- a0_4 ... dispersion + merged ...   |   
|                        H                                                      | 
| a1_1 --- del1 --- a1_2 |  a1_3 --- del3 --- a1_4                              |
 */

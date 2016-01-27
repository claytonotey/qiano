#ifndef DWGS_H
#define DWGS_H

#include "filter.h"

class dwgs;

enum {
  DelaySize = 4096,
  nMaxLongModes = 32
};


enum {
  NoBottomDispersion = 0,
  NoTopDispersion,
  BothDispersion
};

class dwgs {
 public:
  dwgs(float Fs, int type);
  ~dwgs();

  int getHammerNutDelay();
  void set(int upsample, float f, float c1, float c3, float B, float L, float longFreq1, float gammaL, float gammaL2, float inpos, float Z, int hammerNutDelay);
  float input_velocity();
  float go_string();
  float go_soundboard(float hload, float sbload);
  void damper(float c1, float c3, float gammaL, float gammaL2);

  float longTran(dwgs *top, dwgs *bottom);
  void tran2long(dwgs *top, dwgs *bottom, bool bD);
  void longForce(float *F, int n, int decimation);

  int iMax;

  float *wave;
  float *wave1;
  float *dwave;
  float *Fl;
  float *dFl;

  float L;
  int type;
  float f,Fs;
  float inpos;
  float B;
  float longFreq1;
  int nLongModes;
  int upsample;

  Loss_ *loss;
  Thiran fracdelay;
  ThiranDispersion dispersion[4];

  float hermite_p0[nMaxLongModes];
  float hermite_p1[nMaxLongModes];
  float hermite_m0[nMaxLongModes];
  float hermite_m1[nMaxLongModes];
  float *modeTable[nMaxLongModes];
  float fLong[nMaxLongModes];
  DWGResonator longModeResonator[nMaxLongModes];
  int decimation;
  int t;

  float a0_1, a0_2, a0_3, a0_4;
  float a1_1, a1_2, a1_3, a1_4;
  int del1, del2, del3;
  Delay<DelaySize> d0;
  Delay<DelaySize> d1;
  Delay<DelaySize> d2;
  Delay<DelaySize> d3;
  
 
  int M;
};

#endif

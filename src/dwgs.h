#ifndef DWGS_H
#define DWGS_H

#include "filter.h"

class dwgs;

enum {
  DelaySize = 4096,
  nMaxLongModes = 128
};

class dwgs {
 public:
  dwgs(float Fs);
  ~dwgs();
  
  int getMaxDecimation(float f, int decimation, float magic);
  int getMinUpsample(float f, float inpos);
  void set(int upsample, float f, float c1, float c3, float B, float L, float longFreq1, float gammaL, float gammaL2, float inpos, float Z);
  float input_velocity();
  float next_input_velocity();
  float go_string();
  float go_soundboard(float hload, float sbload);
  void damper(float c1, float c3, float gammaL, float gammaL2);

  float longTran(dwgs *top, dwgs *bottom);
  void tran2long(dwgs *top, dwgs *bottom, int offset, int delay);
  void longForce(float *F, int n, int decimation);

  int delTab;
  float *wave;
  float *wave0;
  float *wave2;
  float *wave1;
  float *dwave;
  float *Fl;
  float *dFl;

  float L;
  float omega;
  float f,Fs;
  float inpos;
  float B;
  float longFreq1;
  int nLongModes;
  int upsample;

  Loss loss;
  Thiran fracDelayTop;
  Thiran fracDelayBottom;
  Thiran hammerDelay;
  ThiranDispersion dispersion[4];

  float hermite_p0[nMaxLongModes];
  float hermite_p1[nMaxLongModes];
  float hermite_m0[nMaxLongModes];
  float hermite_m1[nMaxLongModes];
  float *modeTable[nMaxLongModes];
  float fLong[nMaxLongModes];
  DWGResonator longModeResonator[nMaxLongModes];

  float dBottomAndLoss;
  float a0_1, a0_2, a0_3, a0_4, a0_5;
  float a1_1, a1_2, a1_3, a1_4;
  int del0, del1, del2, del3, del4, del5;
  Delay<DelaySize> d0;
  Delay<DelaySize> d1;
  Delay<DelaySize> d2;
  Delay<DelaySize> d3;
  
 
  int M;
};

#endif

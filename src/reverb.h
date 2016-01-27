#ifndef REVERB_H
#define REVERB_H
#include <stdio.h>
#include "filter.h"
#include "FFTConvolver.h"
#include "TwoStageFFTConvolver.h"
using namespace fftconvolver;

#include <algorithm>
using namespace std;

enum { revSize = 65536 };

template<int size>
class ConvolveReverb {

  enum {
    delaySize = 2*size,
    resSize = size,
    tailBlockSize = 512
  };

  int headBlockSize;

 public:
  static float res[resSize] __attribute__ ((aligned(32)));
  int k;
  Delay<delaySize> d;
  //TwoStageFFTConvolver fftConvolver;
  FFTConvolver fftConvolver;


 ConvolveReverb(int blockSize) : k(0) {
    headBlockSize = min(blockSize, (int)tailBlockSize);
    //fftConvolver.init(headBlockSize, tailBlockSize, res, size);
    fftConvolver.init(headBlockSize, res, size);
  } 

  void fft_conv(float *in, float *out, int N) {
    //fprintf(stderr,"%d %d\n",headBlockSize,N); 
    fftConvolver.process(in,out,N);
  }

};


enum {
  NumLengths = 18
};


class Reverb {
public:
  Reverb(float Fs);
  void set(float size, float c1, float c3);
  float reverb(float in); 
  float probe();

 protected:
  int getLength(int k);
  static int allLengths[NumLengths];
  int lengths[8];
  float Fs;
  Delay<1024> d[8];
  float o[8];
  float b[8];
  float c[8];
  Loss<1> decay[8];
  float out;
};

#endif

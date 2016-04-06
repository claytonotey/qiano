#include "qiano.h"
#include "vst.h"
#include "mex.h"
#include "matrix.h"


#define X plhs[0]

#define NOTE prhs[0]
#define LENGTH prhs[1]
#define VEL prhs[2]
#define TUNE prhs[3]
#define ARGS prhs[4]

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  if(nrhs < 3) { 
    mexErrMsgTxt("Too few input arguments."); 
  } else if(nrhs > 5) { 
    mexErrMsgTxt("Too many input arguments."); 
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments."); 
  }

  int note = (int) mxGetScalar(NOTE);
  int N = (int) mxGetScalar(LENGTH);
  float velocity = mxGetScalar(VEL);

  double *tune = mxGetPr(TUNE);

  X = mxCreateDoubleMatrix(N,1,mxREAL);
  double *x = mxGetPr(X);

  float ftune[128*3];
  float *fx = new float[N];
  if(tune) {
    for(int i=0; i<128*3; i++) {
      ftune[i] = tune[i];
    }
  }

  int blockSize = 512;
  int Fs = 44100;
  Piano *effect = (Piano*)createEffectInstance(vst2xPluginHostCallback);
  effect->init(Fs,blockSize);

  if(nrhs == 5) {
    if(!mxIsCell(ARGS)) {
      mexErrMsgTxt("5th argument must be a cell array of key-value pairs.");
    }
    int nargs = mxGetNumberOfElements(ARGS);
    if(nargs%2) {
      mexErrMsgTxt("5th argument must be a cell array of key-value pairs.");
    }
    
    for(int i=0; i<nargs; i+=2) {
      mxArray *keyArr = mxGetCell(ARGS, i);
      if(!mxIsChar(keyArr)) {
        mexErrMsgTxt("5th argument must be a cell array of key-value pairs.");
      }
      int buflen =  mxGetNumberOfElements(keyArr)+1;
      char *key = (char*)mxCalloc(buflen, sizeof(char));
      int status = mxGetString(keyArr, key, buflen); 
      
      fprintf(stderr,"len %d\n",buflen);
      mxArray *valArr = mxGetCell(ARGS, i+1);
      float val = mxGetScalar(valArr);

      int index = getParameterIndex(key);
      if(index >= 0) {
        effect->setParameterLiteral(index,val);
      }
      fprintf(stderr,"%s = %g %d\n",key,val,index);
    }
  }

  float *ftune2 = ftune+3*(note-1);
  effect->triggerOn(note,velocity,ftune2);


  for(int i=0; i<N; i += blockSize) {
    int block = min(blockSize,N - i);
    effect->process(fx+i,block);
  }

  for(int i=0; i<N; i++) {
    x[i] = fx[i];
  }
}

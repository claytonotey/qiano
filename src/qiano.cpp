#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "qiano.h"
#include <string.h>
#include "midi.h"
#include "vst.h"
#include "utils.h"
#include <xmmintrin.h>

// The dwgs uses velocity waves dy/dt,  which reflect with -1
// These wave directly interacts with the hammer
// which has incoming string impedance Z = sqrt(T mu)
// velocity waves are converted to slope waves as
// s = dy/dx = dy/dt / sqrt(T/mu)
// Total velocity = sTop + sBottom
// Total slope = sTop - sBottom

//XXX delay clear opt

float TUNE[3][3] = {{1.0, 0.0, 0.0 },
                    {0.9997, 1.0003, 0.0},
                    { 1.0001, 1.0003, 0.9996 }};    


Param params[NumParams] = {
  { pYoungsModulus, "E", "GPa"},
  { pStringDensity, "rho", "kg / m^3" },
  { pHammerMass, "m", "kg" },
  { pStringLength, "L", "m" },
  { pStringRadius, "r", "m" },
  { pHammerCompliance, "p", "" },
  { pHammerSpringConstant, "K", "kg / s^2" },
  { pHammerHysteresis, "alpha", "" },
  { pBridgeImpedance, "Zb", "" },
  { pBridgeHorizontalImpedance, "ZbH", "" },
  { pVerticalHorizontalImpedance, "Zvh", "" },
  { pHammerPosition, "pos", "" },
  { pSoundboardSize, "size", "" },
  { pStringDecay, "c1", "" },
  { pStringLopass, "c3", "" },
  { pDampedStringDecay, "d1", "" },
  { pDampedStringLopass, "d3", "" },
  { pSoundboardDecay, "s1", "" },
  { pSoundboardLopass, "s3", "" },
  { pLongitudinalGamma, "gammaL", "" },
  { pLongitudinalGammaQuadratic, "gammaL2", "" },
  { pLongitudinalGammaDamped, "gammaLDamped", "" },
  { pLongitudinalGammaQuadraticDamped, "gammaL2Damped", "" },
  { pLongitudinalMix, "lmix", "" },
  { pVolume, "volume", "" },
  { pMaxVelocity, "maxv", "m/s" },
  { pStringDetuning, "detune", "%" },
  { pBridgeMass, "mb", "kg" },
  { pBridgeSpring, "kb", "kg/s^2" },
  { pDecimation, "decimation", "" }
};

int getParameterIndex(const char *key) {
  for(int i=0; i<NumParams; i++) {
    if(!strcmp(key,params[i].name)) return i;
  }
  return -1;
}

double PianoNote :: freqTable[NUM_NOTES];

void PianoNote :: fillFrequencyTable() {
  double NOTE_UP_SCALAR = pow(2.0,1.0/12.0); 
  double A = 6.875;	// A
  A *= NOTE_UP_SCALAR;	// A#
  A *= NOTE_UP_SCALAR;	// B
  A *= NOTE_UP_SCALAR;	// C, frequency of midi note 0
  for (int i = 0; (i < NUM_NOTES); i++)	// 128 midi notes
    {
      freqTable[i] = A;
      A *= NOTE_UP_SCALAR;
    }
}

float PianoNote :: goUpsampled()
{
  float downSampleOutput;
  bool bLoop;

  while(downSampleDelayNeeded) {
    float output = go();
    downSampleFilter.filter(output);
    downSampleDelayNeeded--;
  }

  for(int s = 0; s<upsample; s++) {
    float output = go();
    if(s == 0) {
      downSampleOutput = downSampleFilter.filter(output);
    } else {
      downSampleFilter.filter(output);
    }
  }
  
  
  //printf("%g\n", downSampleOutput);
  return downSampleOutput;
}

float PianoNote :: go()
{
  float output;
  Value *v = piano->vals;

  while(nSamplesReady == 0) {
    float vstringT0 = 0.0;
    float vstringT1 = 0.0;
    
    for(int k=0;k<nstrings;k++) {
      vstringT0 += stringT[k]->input_velocity();
      vstringT1 += stringT[k]->next_input_velocity();
    }
    float vin0 = vstringT0 * nstringsi;
    float vin1 = vstringT1 * nstringsi;
    float hload = hammer->load(vin0,vin1)*Z2i;

    /* The note's string outputs sum to give a load at the bridge
       The bridge is treated seperately for each note,
       the sum of all incoming waves at the bridge becomes the 
       load for the soundobard.
       The magnitude of the tarnsverse oscillations is set by K
       The output of the soundboard reflects back into the bridge junction,
       and the soundboard output is also the output of the qiano 
    */
    
    float sbloadT = 0.0;
    float sbloadHT = 0.0;
    
    for(int k=0;k<nstrings;k++) {
      sbloadT += stringT[k]->go_string();
      sbloadHT += stringHT[k]->go_string();
    }
      
    float in[2] = {sbloadT, sbloadHT};
    float out[2];
    bridge.filter(in, out);
    
    out[0] *= 2.0 * Z;
    out[1] *= 2.0 * Z;
    
    //printf("%g %g %g %g\n",in[0], out[0], in[1], out[1]);
    //out[0] = sbloadT * alphasb;
    //out[1] = 0; 

    float tranForce = 0.0f;
    float tranForceH = 0.0f;
    float longTranForce = 0.0f;
    for(int k=0;k<nstrings;k++) {
      tranForce += stringT[k]->go_soundboard(hload,out[0]);
      tranForceH += stringHT[k]->go_soundboard(0,out[1]);
      longTranForce += stringT[k]->longTran(stringT[k],stringT[k]);
    }
    
    longTranForces[tTran] = longTranForce; 
    tranForces[tTran] = tranForce;
    tranForcesH[tTran] = tranForceH;
    tTran = (tTran + 1)%TranBufferSize;
    
    bool bInterpReady = false;
    if(tLong >= 0) { 
      for(int k=0;k<nstrings;k++) {
        // XXX ignoring horizontal -> longitudinal
        // XXX the derivative is wrong for the first sample  (should be forward derivative not centered)
        

        if(tLong % decimation == 0) {
          stringT[k]->tran2long(stringT[k],stringT[k], 0, longDelay);        
        } else if((tLong-1) % decimation == 0) {
          stringT[k]->tran2long(stringT[k],stringT[k], 1, longDelay);           
        } else if((tLong+1) % decimation == 0) {
          stringT[k]->tran2long(stringT[k],stringT[k], -1, longDelay);
        }

      }
    }

    tLong++;

    if(tLong == nLongNeeded) {
      int nSamples = decimation;
      memset(longForces,0,nSamples*sizeof(float));
      for(int k=0;k<nstrings;k++) {
        stringT[k]->longForce(longForcesK, nSamples, decimation);
        for(int i=0; i<nSamples; i++) {
          longForces[i] += longForcesK[i];
        }
      }
      nSamplesReady = nSamples;
      tLong = 0;
      if(hammer->isEscaped()) {
        decimation = decimationNoHammer;
      }
      nLongNeeded = decimation;
    }
  }
  //printf("%g %g %g %g\n", tranForces[t], tranForcesH[t], longTranForces[t], longForces[t]);
      
  //longTranBridgeForce = 0;
  output = tranForces[tTranRead] * tranBridgeForce + tranForcesH[tTranRead] * tranBridgeForce + longHP.filter(longTranForces[tTranRead] * longTranBridgeForce + longForces[tLongRead] * longBridgeForce) * v[pLongitudinalMix];
  //output = tranForces[t] * tranBridgeForce;
  float delayed = outputdelay.goDelay(output);
  energy = energy + output*output - delayed*delayed;
  if(energy>maxEnergy)
    maxEnergy = energy;

  tTranRead = (tTranRead + 1)%TranBufferSize;
  tLongRead++;
  if(tLongRead == nSamplesReady) {
    tLongRead = 0;
    nSamplesReady = 0;
  }

  return output;
}

void Piano :: addVoice(PianoNote *v) 
{ 
  if(voiceList) {
    v->next = voiceList;
    v->prev = voiceList->prev;
    voiceList->prev->next = v;
    voiceList->prev = v;
  } else {
    voiceList = v;
    v->prev = v;
    v->next = v;
  }
}

void Piano :: removeVoice(PianoNote *v)
{	
	if(v == voiceList) {
		if(v == v->next)
			voiceList = NULL;
		else
			voiceList = v->next;
	}
	PianoNote *p = v->prev;
	PianoNote *n = v->next;  
	p->next = n;
	n->prev = p;
}

void Piano :: init(float Fs, int blockSize) 
{
  this->Fs = Fs;
  this->blockSize = blockSize;
  voiceList = NULL;
  PianoNote :: fillFrequencyTable(); 
  for(int k=PIANO_MIN_NOTE; k<=PIANO_MAX_NOTE;k++) {
	  if(noteArray[k]) delete noteArray[k];
    noteArray[k] = new PianoNote(k,Fs,this);
  }
  
  if(input) delete input;
  input = new float[blockSize];

  
  if(soundboard) delete soundboard;
#ifdef FDN_REVERB
  soundboard = new Reverb(Fs);
#else
  soundboard = new ConvolveReverb<revSize>(blockSize);
#endif

  setParameter(pYoungsModulus, 0.5);
  setParameter(pStringDensity, 0.5);
  setParameter(pHammerMass, 0.5);
  setParameter(pStringLength, 0.25);
  setParameter(pStringRadius, 0.25);
  setParameter(pHammerCompliance, 0.5);
  setParameter(pHammerSpringConstant, 0.5);
  setParameter(pHammerHysteresis, 0.5);
  setParameter(pHammerPosition, 2.0/7.0);
  setParameter(pBridgeImpedance, 0.5);
  setParameter(pBridgeHorizontalImpedance, 0.5);
  setParameter(pVerticalHorizontalImpedance, 0.5);
  setParameter(pSoundboardSize, 0.0);
  setParameter(pStringDecay, 0.25);
  setParameter(pStringLopass, 0.5);
  setParameter(pDampedStringDecay, 0.5);
  setParameter(pDampedStringLopass, 0.5);
  setParameter(pSoundboardDecay, 0.5);
  setParameter(pSoundboardLopass, 0.5);
  setParameter(pLongitudinalGamma, 0.5);
  setParameter(pLongitudinalGammaQuadratic, 0.0);
  setParameter(pLongitudinalGammaDamped, 0.5);
  setParameter(pLongitudinalGammaQuadraticDamped, 0.0);
  setParameter(pLongitudinalMix, 0.5);
  setParameter(pVolume, 0.5);
  setParameter(pMaxVelocity, 0.5);
  setParameter(pStringDetuning, 0.5);
  setParameter(pBridgeMass, 0.5);
  setParameter(pBridgeSpring, 0.5);
  setParameter(pDecimation, 0.5);
}

Piano :: Piano(audioMasterCallback audioMaster, int parameters) : VstEffect(audioMaster, parameters){
  
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  setNumInputs (0);		// stereo in
  setNumOutputs (2);		// stereo out
  setUniqueID (8881887);	// identify
  isSynth(true);
  for(int k=PIANO_MIN_NOTE; k<=PIANO_MAX_NOTE;k++) {
    noteArray[k] = NULL;
  }
  input = NULL;
  soundboard = NULL;
}

PianoNote :: PianoNote(int note, int Fs, Piano *piano) {
  this->Fs = Fs;
  this->note = note;
  this->f = freqTable[note];
  this->piano = piano;

  if(note<31)
    nstrings = 1;
  else if(note<41)
    //XXX
    nstrings = 2;
  else 
    nstrings = 3;

  //nstrings = 1;
  nstringsi = 1.0/(float)nstrings;
  
  for(int k=0;k<nstrings;k++) {
    stringT[k] = new dwgs(Fs);
    stringHT[k] = new dwgs(Fs);
  }
  hammer = new Hammer(Fs);

  longHP.create(0.15,0.5);

  outputdelay.setDelay(Fs);
  energy = 0.0;
  maxEnergy = 0.0;
  bActive = false;
}

bool PianoNote :: isDone()
{
	return (energy<1e-6*maxEnergy);
}

void PianoNote :: triggerOn(float velocity, float *tune) 
{
  //fprintf(stderr,"note = %d velocity =  %g\n",note,velocity);
  Value *v = piano->vals;
  float f0 = 27.5; //A0 note 21
  float L = .14 + 1.4/(1+exp(-3.4+1.4*log(f/f0)));
  L = 0.075 + 2.0/(1.0+exp(-4.8+2.1*log(f/f0)));
  L *= v[pStringLength];
  //L = .115;
  float p = 2.0+1.0*log(f/f0)/log(4192/f0);
  p *= v[pHammerCompliance];


  //float m = .06 * (1.0 - 0.9*pow((float)log(f/f0)/log(4192/f0),(float)0.1));
  //m=.018 * (1.0 - .7*pow(log(f/f0)/log(4192/f0),0.8));
  float m = .02 * pow(0.5 - log(f/4192),0.3);
  m = .013 - .005 *log(f/f0)/log(4192/f0);
  m *= v[pHammerMass];
  //m = .0092;

  //float K,p;
  float K = 40/pow((float).7e-3,(float)p);
  float alpha = 0.1e-4*log(f/f0)/log(4192/f0);
  //K *= 1.;

  int N = note - 20;
  float eps =  0.9894 + 8.8e-5 * N;
  //p = 3.7 + .015 * N;
  //p = 3.7;
  //K = 15.5e3 / pow(1e-3,p) * exp(.059*N) * (1 - eps);
  //K = 183.0 / pow(1e-3,p) * exp(.045*N);
  //float tau = 1e-6 * (2.72 - 0.02 * N + 9e-5 * N * N);
  //alpha = tau / (1 - eps);
  //alpha = 1e-6 * (148 + 1.83 * N - 5.5e-2 * N*N + 8.5e-4 * N*N*N);
  K *= v[pHammerSpringConstant];
  alpha *= v[pHammerHysteresis];


  float r = .0025*pow((float)(0.9+.6*log(f/f0)),(float)-1.4);
  //r = .000475;
  r *= v[pStringRadius];
  float S = PI * r * r;
  mu = S*v[pStringDensity];
  T = (2*L*f)*(2*L*f)*mu;
  float rcore = (r<.0006)?r:.0006;
  float Score = PI * rcore * rcore;
  float mucore = Score * v[pStringDensity];
  Z = sqrt(T*mu);
  float E = v[pYoungsModulus] * 1e9;
  float B = (PI*PI*PI)*E*(rcore*rcore*rcore*rcore)/(4.0*L*L*T);
  //B *= 5;

  //float vLong = sqrt(E*S/mu);
  float vLong = sqrt(E*S/mu);
  vTran = sqrt(T/mu);
  float longFreq1 = vLong / (2 * L);
  //float Zlong = sqrt(E*S*mu);
  float Zlong = sqrt(E*Score*mucore);

  tranBridgeForce = Z;
  longBridgeForce = Zlong*square(vLong/vTran) / L / Fs;
  longTranBridgeForce = 0.5 * Zlong * (vLong / vTran) / vTran;
  
  Z2i = 1.0/(2*Z);
  alphasb = (2*Z)/(Z*nstrings + v[pBridgeImpedance]);

  float Zb = v[pBridgeImpedance];
  float hp = v[pHammerPosition];
  //hp = .0081/L;


  float mBridge = v[pBridgeMass];
  float kBridge = v[pBridgeSpring];
  

  //fprintf(stderr,"f = %g, r = %g mm, L = %g, T = %g, hpos = %g, hm = %g, Z = %g, K = %g, B = %g, Zb = %g, alpha = %g, p = %g, vTran=%g, vLong=%g, mBridge=%G, kBridge=%g\n",f,1000*r,L,T,hp,m,Z,K,B,Zb,alpha,p,vTran,vLong,mBridge,kBridge);
 
  float ZbH = v[pBridgeHorizontalImpedance];
  float Zhv = v[pVerticalHorizontalImpedance];
  float khv = kBridge * 0.;

  

  float ES = E*S;
  float gammaL = v[pLongitudinalGamma];
  float gammaL2 = v[pLongitudinalGammaQuadratic];


  int upsampleMin = 1;
  int decimationMax = MaxDecimation;
	for(int k=0;k<nstrings;k++) {
    float fk;
    if(tune) {
      fk = f*(1 + (tune[k]-1) * v[pStringDetuning]);
    } else {
      fk = f*(1 + (TUNE[nstrings-1][k]-1) * v[pStringDetuning]);
    }
    upsampleMin = max(upsampleMin, stringT[k]->getMinUpsample(fk, hp));
  }
  upsample = upsampleMin;

	for(int k=0;k<nstrings;k++) {
    float fk;
    if(tune) {
      fk = f*(1 + (tune[k]-1) * v[pStringDetuning]);
    } else {
      fk = f*(1 + (TUNE[nstrings-1][k]-1) * v[pStringDetuning]);
    }
    decimationMax = min(decimationMax, stringT[k]->getMaxDecimation(fk, upsample,v[pDecimation]));
  }
  decimation = decimationMax;

  bridge.create(upsample*Fs,
                mBridge,kBridge,Zb,
                mBridge,kBridge,ZbH,
                Zhv,khv,
                Z*nstrings);

  hammer->set(upsample,m,K,p,Z,alpha);

	for(int k=0;k<nstrings;k++) {
    float fk;
    if(tune) {
      fk = f*(1 + (tune[k]-1) * v[pStringDetuning]);
    } else {
      fk = f*(1 + (TUNE[nstrings-1][k]-1) * v[pStringDetuning]);
    }
    stringT[k]->set(upsample,fk,v[pStringDecay],v[pStringLopass],B,L,longFreq1,gammaL,gammaL2,v[pHammerPosition],Z);
    stringHT[k]->set(upsample,fk,v[pStringDecay],v[pStringLopass],B,L,longFreq1,gammaL,gammaL2,v[pHammerPosition],Z);
	}
	hammer->strike(velocity);
	maxEnergy = 0.0;
	energy = 0.0;
  
	outputdelay.clear();
	bActive = true;
  longDelay = 10;
  t = 0;
  tTran = 0;
  tLong = -longDelay;
  tLongRead = 0;
  tTranRead = 0;
  nSamplesReady = 0;
  decimationNoHammer = decimation;
  decimation = 1;
  nLongNeeded = decimation + 2;
  
  //nLongNeeded = 3;
  downSampleFilter.create(upsample);
  downSampleDelayNeeded = downSampleFilter.getDelay();

  fprintf(stderr,"upsample/decimation  %d %d\n",upsample, decimationNoHammer);
}

void PianoNote :: triggerOff() 
{
  Value *v = piano->vals;
  float gammaLDamped = v[pLongitudinalGammaDamped];
  float gammaL2Damped = v[pLongitudinalGammaQuadraticDamped];
	for(int k=0;k<nstrings;k++) {
		stringT[k]->damper(v[pDampedStringDecay],v[pDampedStringLopass],gammaLDamped,gammaL2Damped);
		stringHT[k]->damper(v[pDampedStringDecay],v[pDampedStringLopass],gammaLDamped,gammaL2Damped);
	}
}

void PianoNote :: deActivate()
{
	bActive = false;
}

bool PianoNote :: isActive() 
{
	return bActive;
}


PianoNote :: ~PianoNote() {
 for(int k=0;k<nstrings;k++) {
    delete stringT[k];
    delete stringHT[k];
  } 
  delete hammer;
}

Piano :: ~Piano() { 
  for(int k=PIANO_MIN_NOTE;k<=PIANO_MAX_NOTE;k++)
	  delete noteArray[k];
  delete input;
  delete soundboard;
}

void Piano :: process(float *out, int samples) 
{
  for(int i=0;i<samples;i++) {
    PianoNote *v = voiceList;
    float output = 0;
    do {
      if(v) output += v->goUpsampled();
    } while(v && (v=v->next) && (v!= voiceList));  
#ifdef FDN_REVERB
    out[i] = vals[pVolume] * soundboard->reverb(output);
#else
    //    if(isnan(output)) abort();
    out[i] =  vals[pVolume] * output;
    input[i] =  vals[pVolume] * output;
#endif
  }


#ifndef FDN_REVERB
  soundboard->fft_conv(input, out, samples);
#endif
}

void Piano :: process(float **in, float **out, int samples, int offset) 
{
	process(out[0]+offset,samples);
	memcpy(out[1]+offset, out[0]+offset, samples * sizeof(float));
}

void Piano :: process(float **inS, float **outS, int sampleFrames) 
{ 
  int delta = 0;
  BlockEvents end;
  end.delta = sampleFrames;
  end.eventStatus = isNull;
  blockEvents[numBlockEvents++] = end;    
 

  PianoNote *v = voiceList;
  PianoNote *remove[NUM_NOTES];
  int k=0;
  int n=0;
  do {
	  if(v && v->isDone()) {
		  v->deActivate();	    
		  remove[k++] = v;
	  }	 
  } while(v && (v=v->next) && (v!=voiceList));  
  for(int j=0;j<k;j++)
    removeVoice(remove[j]);
	
  for(int i=0;i<numBlockEvents;i++) {    
	
    BlockEvents event = blockEvents[i];
    int nextDelta = event.delta;
	
	if(event.byte1>=PIANO_MIN_NOTE && event.byte1<=PIANO_MAX_NOTE) {
		if(event.eventStatus == isNote) {

			PianoNote *v = noteArray[event.byte1];		
			if(event.byte2) {		
				if(!(noteArray[event.byte1]->isActive())) {
					addVoice(v);			
        }
        float velocity = vals[pMaxVelocity]*pow((float)(event.byte2/127.0),(float)2.0);
				v->triggerOn(velocity,NULL);
			} else {	
				if(v) {
          v->triggerOff();
				}
			}
		}
	}
  

  nextDelta = std::min(nextDelta, sampleFrames);
  process(inS, outS, nextDelta - delta, delta);
  
	delta = nextDelta;
  }  
  numBlockEvents = 0; 
}

void Piano :: triggerOn(int note, float velocity, float *tune)
{
  PianoNote *v = noteArray[note];
  addVoice(v);
  v->triggerOn(velocity,tune);
}

#include "audioeffectx.h"

AudioEffect* createEffectInstance (audioMasterCallback audioMaster)
{
  return new Piano(audioMaster,NumParams);
}

void Piano::resume()
{
  VstEffect::resume();
	//wantEvents(1);
  blockSize = updateBlockSize();
  Fs = updateSampleRate();
  init(Fs,blockSize);
}


//-----------------------------------------------------------------------------------------
void Piano::setParameterLiteral (VstInt32 index, float value) {
  Value &p = vals[index];
  p.v = value;
}

void Piano::setParameter (VstInt32 index, float value)
{
  Value &p = vals[index];
  p.f = value;

  switch(index) {
  case pYoungsModulus:
    p.v = 200 * exp(4.0 * (value - 0.5));
    break;
  case pStringDensity:
    p.v = 7850 * exp(4.0 * (value - 0.5));
    break;
  case pHammerMass:
    p.v = exp(4.0 * (value - 0.5));
    break;
  case pStringLength:
    p.v = exp(2.0 * (value - 0.25));
    break;
  case pStringRadius:
    p.v = exp(2.0 * (value - 0.25));
    break;
  case pHammerCompliance:
    p.v = 2.0 * value;
    break;
  case pHammerSpringConstant:
    p.v = (2.0 * value);
    break;
  case pHammerHysteresis:
    p.v = exp(4.0 * (value-0.5));
    break;
  case pBridgeImpedance:
    p.v = 4000.0 * exp(12.0 * (value - 0.5));
    break;
  case pBridgeHorizontalImpedance:
    p.v = 40000.0 * exp(12.0 * (value - 0.5));
    break;
  case pVerticalHorizontalImpedance:
    p.v = 400.0 * exp(12.0 * (value - 0.5));
    break;
  case pHammerPosition:
    p.v = value * 0.5;
    break;
  case pSoundboardSize:
    p.v = value;
#ifdef FDN_REVERB
    soundboard->set(vals[pSoundboardSize],vals[pSoundboardDecay],vals[pSoundboardLopass]);
#endif
    break;
  case pStringDecay:
    p.v = 0.25 * exp(6.0 * (value - 0.25));
    break;
  case pStringLopass:
    p.v = 5.85 * exp(6.0 * (value - 0.5));
    break;
  case pDampedStringDecay:
    p.v = 8.0 * exp(6.0 * (value - 0.5));
    break;
  case pDampedStringLopass:
    p.v = 25.0 * exp(6.0 * (value - 0.5));
    break;
  case pSoundboardDecay:
    p.v = 20.0 * exp(4.0 * (value - 0.5));
#ifdef FDN_REVERB
    soundboard->set(vals[pSoundboardSize],vals[pSoundboardDecay],vals[pSoundboardLopass]);
#endif
    break;
  case pSoundboardLopass:
    p.v = 20.0 * exp(4.0 * (value - 0.5));
#ifdef FDN_REVERB
    soundboard->set(vals[pSoundboardSize],vals[pSoundboardDecay],vals[pSoundboardLopass]);
#endif
    break;
  case pLongitudinalGamma:
    p.v = 1e-2 * exp(10.0 * (value - 0.5));
    break;
  case pLongitudinalGammaQuadratic:
    p.v = 2.0e-2 * value;
    break;
  case pLongitudinalGammaDamped:
    p.v = 5e-2 * exp(10.0 * (value - 0.5));
    break;
  case pLongitudinalGammaQuadraticDamped:
    p.v = 4.0e-2 * value;
    break;
  case pLongitudinalMix:
    p.v = (value==0.0)?0.0:1e0 * exp(16.0 * (value - 0.5));
    break;
  case pVolume:
    p.v = 5e-3 * exp(8.0 * (value - 0.5));
    break;
  case pMaxVelocity:
    p.v = 10 * exp(8.0 * (value - 0.5));
    break;
  case pStringDetuning:
    p.v = 1.0 * exp(10.0 * (value - 0.5));
    break;
  case pBridgeMass:
    p.v = 10.0 * exp(10.0 * (value - 0.5));
    break;
  case pBridgeSpring:
    p.v = 1e5 * exp(20.0 * (value - 0.5));
    break;
  case pDecimation:
    p.v = (value * 200);
    break;
  }

}

//-----------------------------------------------------------------------------------------
float Piano::getParameter (VstInt32 index)
{
  return vals[index].f;
}

//-----------------------------------------------------------------------------------------
void Piano::getParameterName (VstInt32 index, char* label)
{
  vst_strncpy (label, params[index].name, kVstMaxParamStrLen);
}

//-----------------------------------------------------------------------------------------
void Piano::getParameterDisplay (VstInt32 index, char* text)
{
  float2string (vals[index], text, kVstMaxParamStrLen);
}

//-----------------------------------------------------------------------------------------
void Piano::getParameterLabel (VstInt32 index, char* label)
{
  vst_strncpy (label, params[index].label, kVstMaxParamStrLen);
}

//------------------------------------------------------------------------
bool Piano::getEffectName (char* name)
{
  vst_strncpy (name, "Qiano", kVstMaxEffectNameLen);
  return true;
}

//------------------------------------------------------------------------
bool Piano::getProductString (char* text)
{
  vst_strncpy (text, "Qiano", kVstMaxProductStrLen);
  return true;
}

//------------------------------------------------------------------------
bool Piano::getVendorString (char* text)
{
  vst_strncpy (text, "Mune", kVstMaxVendorStrLen);
  return true;
}

//-----------------------------------------------------------------------------------------
VstInt32 Piano::getVendorVersion ()
{ 
  return 1000; 
}

VstInt32 Piano::canDo(char* text)
{
	if (strcmp(text, "receiveVstEvents") == 0)
		return 1;
	if (strcmp(text, "receiveVstMidiEvent") == 0)
		return 1;	
	return -1;
}


int main(int c, char **v) 
{
  int note = atoi(v[1]);
  int t = atof(v[2]);
  int vel = atoi(v[3]);
  int N = atoi(v[4]);
  float magic = atof(v[5]);

  float *in[2];
  float *out[2];

  in[0] = new float[t];
  in[1] = new float[t];
  out[0] = new float[t];
  out[1] = new float[t];

  AudioEffect *effect = createEffectInstance(vst2xPluginHostCallback);
  int blockSize = 512;
  effect->setBlockSize(blockSize);
  effect->resume();
  effect->setParameter(pDecimation,magic);

  
  for(int k=0; k<N; k++) {
    MidiEvent midiEvents[2];
    midiEvents[0].noteOn(note,vel,0,0);
    midiEvents[1].noteOff(note,vel,0,t/2);
    
    processMidiEvents(effect, midiEvents, 2);


    for(int i=0; i<t; i += blockSize) {
      int block = std::min(blockSize,t - i);
      //fprintf(stderr,"block = %d %d\n",block,blockSize);
      float *in2[2];
      float *out2[2];
      in2[0] = in[0] + i;
      in2[1] = in[1] + i;
      out2[0] = out[0] + i;
      out2[1] = out[1] + i;

      processReplacing(effect, in2,out2,block);
    }
  }

  for(int k=0; k<t; k++) {
    //printf("%g\n",out[0][k]);
  }
}




extern "C" void qianoNote(int note, int N, float velocity, float *x, float *tune)
{
  int blockSize = 512;
  int Fs = 44100;
  Piano *effect = (Piano*)createEffectInstance(vst2xPluginHostCallback);
  effect->init(Fs,blockSize);
  effect->triggerOn(note,velocity,tune?(tune+3*(note-1)):NULL);


  float in[512];

  for(int i=0; i<N; i += blockSize) {
    int block = std::min(blockSize,N - i);
    effect->process(x+i,block);
  }
}



/*
Zhh - h --- Z 
     Zhv
Zvv - v --- Z 

  h --- Z
 Zhv
  v --- Z

*/

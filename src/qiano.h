#ifndef PIANO_H
#define PIANO_H

#define PIANO_MIN_NOTE 22
#define PIANO_MAX_NOTE 108

#undef FDN_REVERB
#define FDN_REVERB 1


#include "filter.h"
#include "dwgs.h"
#include "reverb.h"
#include "hammer.h"
#include "types.h"
#include "midi.h"
#include "vsteffect.h"

enum {
  MaxDecimation = 16,
  MaxLongDelay = 16,
  TranBufferSize = MaxDecimation + MaxLongDelay + 1,
  NumParams = 30
};

class Param {
public:
  int tag;
  char name[64];
  char label[64];
};

class Value {
public:
  operator float() const {
    return v;
  }
  float f;
  float v;
};

enum parameters {
  pYoungsModulus = 0,
  pStringDensity,
  pHammerMass,
  pStringLength,
  pStringRadius,
  pHammerCompliance,
  pHammerSpringConstant,
  pHammerHysteresis,
  pBridgeImpedance,
  pBridgeHorizontalImpedance,
  pVerticalHorizontalImpedance,
  pHammerPosition,
  pSoundboardSize,
  pStringDecay,
  pStringLopass,
  pDampedStringDecay,
  pDampedStringLopass,
  pSoundboardDecay,
  pSoundboardLopass,
  pLongitudinalGamma,
  pLongitudinalGammaQuadratic,
  pLongitudinalGammaDamped,
  pLongitudinalGammaQuadraticDamped,
  pLongitudinalMix,
  pVolume,
  pMaxVelocity,
  pStringDetuning,
  pBridgeMass,
  pBridgeSpring,
  pDecimation
};

int getParameterIndex(const char *key);

class Piano;

class PianoNote {
public:
	PianoNote(int note, int Fs, Piano *piano);
	~PianoNote();
	float go();
	float goUpsampled();
	void triggerOn(float velocity, float *tune);
	void triggerOff();
	bool isActive();
	void deActivate();

	PianoNote *next;
	PianoNote *prev;	
	int note;
  Piano *piano;

	bool isDone();
	float maxEnergy;
	float energy;
	Delay<65536> outputdelay;
  static void fillFrequencyTable();
  static double freqTable[NUM_NOTES];    

protected:

  MSD2Filter bridge;


	bool bActive;
	int Fs;
  BiquadHP longHP;
	float alphasb;
	float Z2i;
	float nstringsi;
  float tranBridgeForce;
  float longBridgeForce;
  float longTranBridgeForce;
  float Z;
  float Zhv;
  float vTran;


  float longForces[MaxDecimation + 1];
  float longForcesK[MaxDecimation + 1];

  float tranForces[MaxLongDelay + MaxDecimation + 1];
  float tranForcesH[MaxLongDelay + MaxDecimation + 1];
  float longTranForces[MaxLongDelay + MaxDecimation + 1];

  
  float T;
  float mu;
	float f;
  int nstrings;

	dwgs *stringT[3];	
	dwgs *stringHT[3];	
  Hammer *hammer;

  int t;
  int tTran;
  int tLong;
  int tTranRead;
  int tLongRead;
  int nLongNeeded;
  int nSamplesReady;
  int decimation;
  int decimationNoHammer;
  int longDelay;
  int upsample;
  int downSampleDelayNeeded;
  DownSampleFIR downSampleFilter;
};

class Piano : public VstEffect
{
public:
	Piano(audioMasterCallback audioMaster, int parameters);
	~Piano();

	void addVoice(PianoNote *v);	
	void removeVoice(PianoNote *v);
	void process(float *out, int samples);	
	void process(float **in, float **out, int frameSamples, int offset);
	virtual void process(float **in, float **out, int frameSamples);
  void init(float Fs, int blockSize);
  void triggerOn(int note, float velocity, float *tune);

	//vst crap
  virtual void resume();
	virtual void setParameterLiteral (VstInt32 index, float value);
	virtual void setParameter (VstInt32 index, float value);
	virtual float getParameter (VstInt32 index);
	virtual void getParameterLabel (VstInt32 index, char* label);
	virtual void getParameterDisplay (VstInt32 index, char* text);
	virtual void getParameterName (VstInt32 index, char* text);

	virtual bool getEffectName (char* name);
	virtual bool getVendorString (char* text);
	virtual bool getProductString (char* text);
	virtual VstInt32 getVendorVersion();
	virtual VstInt32 canDo(char* text);

protected:     
  friend class PianoNote;
  Value vals[NumParams];
	PianoNote *voiceList;
	PianoNote *noteArray[NUM_NOTES];
	float Fs;
  //int blockSize;
  float *input;

#ifdef FDN_REVERB
  Reverb *soundboard;
#else
  ConvolveReverb<revSize> *soundboard;
#endif
};

AudioEffect* createEffectInstance (audioMasterCallback audioMaster);
void qianoNote(int note, int N, float velocity, double *x, double *tune);
extern "C" void qianoNote(int note, int N, float velocity, float *x, float *tune);
	
#endif

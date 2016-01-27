#ifndef VSTEFFECT_H
#define VSTEFFECT_H

#include "audioeffectx.h"
#include "midi.h"

class VstEffect : public AudioEffectX 
{
 public:	 
	 VstEffect(audioMasterCallback audioMaster, int parameters);
    ~VstEffect();
	// Processing
	
	virtual void process(float **in, float **out, int frameSamples)=0;
	virtual void processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames);
	//	virtual void processDoubleReplacing (double** inputs, double** outputs, VstInt32 sampleFrames);
	//2.4
	//virtual void resume();

	// Program
	virtual void setProgramName (char* name);
	virtual void getProgramName (char* name);

	// Parameters
	virtual void setParameter (VstInt32 index, float value)=0;
	virtual float getParameter (VstInt32 index)=0;
	virtual void getParameterLabel (VstInt32 index, char* label)=0;
	virtual void getParameterDisplay (VstInt32 index, char* text)=0;
	virtual void getParameterName (VstInt32 index, char* text)=0;

	virtual bool getEffectName (char* name)=0;
	virtual bool getVendorString (char* text)=0;
	virtual bool getProductString (char* text)=0;
	virtual VstInt32 getVendorVersion()=0;
	virtual VstInt32 canDo(char* text)=0;

	// MIDI
	virtual VstInt32 getNumMidiInputChannels ()  { return 1; }
	virtual VstInt32 getNumMidiOutputChannels () { return 0; }

	virtual VstInt32 processEvents (VstEvents* events);
	void reallyProcessEvents(VstEvents *events, BlockEvents *blockEvents, int *numBlockEvents);
	// bool sendVstEventsToHost (VstEvents* events);
float blob;

	

protected:	
	char programName[kVstMaxProgNameLen+1];
	float **replacingInputs;
  
	// MIDI essentials...
	BlockEvents *blockEvents;
	int numBlockEvents;
	int totalEvents ;
};

#endif

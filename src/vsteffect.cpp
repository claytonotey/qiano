#include "vsteffect.h"
#include <stdio.h>
#include <stdlib.h>

VstEffect :: VstEffect(audioMasterCallback audioMaster, int parameters) : AudioEffectX (audioMaster, 1, parameters)
{
  canProcessReplacing ();	// supports replacing output
  //	canDoubleReplacing ();	// supports double precision processing
  vst_strncpy (programName, "Default", kVstMaxProgNameLen);	// default program name
  
  blockEvents = new BlockEvents[EVENTS_QUEUE_MAX];
  numBlockEvents = 0;
  replacingInputs = (float**) calloc(2,sizeof(float*));
  replacingInputs[0] = (float*) calloc(4096,sizeof(float));
  replacingInputs[1] = (float*) calloc(4096,sizeof(float));

}

VstEffect :: ~VstEffect() 
{
  free(blockEvents);
  free(replacingInputs[0]);
  free(replacingInputs[1]);
  free(replacingInputs);
}

void VstEffect :: processReplacing (float** inputs, float** outputs, VstInt32 samples)
{	

  for(int c=0;c<2;c++) {	
	if(!(cEffect.flags & effFlagsIsSynth))
		memcpy(replacingInputs[c], inputs[c], samples * sizeof(float));
    memset(outputs[c],0,samples*sizeof(float));   
  }  
  
 
  process(replacingInputs,outputs,samples);
}

VstInt32 VstEffect :: processEvents (VstEvents* events) {
  reallyProcessEvents(events, blockEvents, &numBlockEvents);
  
  return 1;
}

void VstEffect::setProgramName (char* name)
{
	vst_strncpy (programName, name, kVstMaxProgNameLen);
}

void VstEffect::getProgramName (char* name)
{
	vst_strncpy (name, programName, kVstMaxProgNameLen);
}

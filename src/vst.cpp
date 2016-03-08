#define VST_FORCE_DEPRECATED 0
#include <stdlib.h>
#include <stdio.h>
#include "audioeffectx.h"
#include "aeffectx.h"
#include "vst.h"

extern "C" {

static int _canHostDo(const char* pluginName, const char* canDoString) 
{
  return 0;
}

VstTimeInfo vstTimeInfo;

VstIntPtr VSTCALLBACK vst2xPluginHostCallback(AEffect *effect, VstInt32 opcode, VstInt32 index, VstIntPtr value, void *dataPtr, float opt) 
{
  int result = 0;
  
  switch(opcode) {
    case audioMasterAutomate:
      break;
    case audioMasterVersion:
      // We are a VST 2.4 compatible host
      result = 2400;
      break;
    case audioMasterCurrentId:
      result = effect->uniqueID;
      break;
    case audioMasterIdle:
      break;
    case audioMasterPinConnected: 
      // Deprecated
      break;
     case audioMasterWantMidi:
      // This is called by old VST2.3 plugins to tell us that they are instruments
      break;
    case audioMasterGetTime:
      // These values are always valid
      vstTimeInfo.samplePos = 0;
      vstTimeInfo.sampleRate = 44100;

      // Set flags for transport state
      vstTimeInfo.flags = 0;

      if(value & kVstPpqPosValid) {
        // TODO: This calculation might be wrong
        double quarterNotesPerMinute = 120.0;
        double millisecondsPerBeat = 1000.0 * 60.0 / quarterNotesPerMinute;
        vstTimeInfo.ppqPos = millisecondsPerBeat / 96.0f;
        vstTimeInfo.flags |= 1;
      }
      if(value & kVstTempoValid) {
        vstTimeInfo.tempo = 120.0;
        vstTimeInfo.flags |= 1;
      }
      if(value & kVstBarsValid) {
        vstTimeInfo.flags |= 0;
      }
      if(value & kVstCyclePosValid) {
        vstTimeInfo.flags |= 0;
      }
      if(value & kVstTimeSigValid) {
        vstTimeInfo.timeSigNumerator = 4.0;
        vstTimeInfo.timeSigDenominator = 4.0;
        vstTimeInfo.flags |= 1;
      }
      if(value & kVstSmpteValid) {
        vstTimeInfo.flags |= 0;
      }
      if(value & kVstClockValid) {
        vstTimeInfo.flags |= 0;
      }

      dataPtr = &vstTimeInfo;
      break;
    case audioMasterProcessEvents:
      break;
    case audioMasterSetTime: // Deprecated
      break;
    case audioMasterTempoAt: // Deprecated
      break;
    case audioMasterGetNumAutomatableParameters: // Deprecated
      break;
    case audioMasterGetParameterQuantization: // Deprecated
      break;
    case audioMasterIOChanged:
      break;
    case audioMasterNeedIdle: // Deprecated
      break;
    case audioMasterSizeWindow:
      break;
    case audioMasterGetSampleRate:
      result = 44100;
      break;
    case audioMasterGetBlockSize:
      result = 512;
      break;
    case audioMasterGetInputLatency:
      result = 0;
      break;
    case audioMasterGetOutputLatency:
      result = 0;
      break;
    case audioMasterGetPreviousPlug: // Deprecated
      break;
    case audioMasterGetNextPlug: // Deprecated
      break;
    case audioMasterWillReplaceOrAccumulate: // Deprecated
      break;
    case audioMasterGetCurrentProcessLevel:
      result = kVstProcessLevelUnknown;
      break;
    case audioMasterGetAutomationState:
      result = kVstAutomationUnsupported;
      break;
    case audioMasterOfflineStart:
      break;
    case audioMasterOfflineRead:
      break;
    case audioMasterOfflineWrite:
      break;
    case audioMasterOfflineGetCurrentPass:
      break;
    case audioMasterOfflineGetCurrentMetaPass:
      break;
    case audioMasterSetOutputSampleRate: 
      // Deprecated
      break;
    case audioMasterGetOutputSpeakerArrangement: // Deprecated
      break;
    case audioMasterGetVendorString:
      strncpy((char *)dataPtr, "Mune", kVstMaxVendorStrLen);
      result = 1;
      break;
    case audioMasterGetProductString:
      strncpy((char *)dataPtr, "Host", kVstMaxProductStrLen);
      result = 1;
      break;
    case audioMasterGetVendorVersion:
      result = 1234;
      break;
    case audioMasterVendorSpecific:
      break;
    case audioMasterCanDo:
      result = _canHostDo(NULL, (char *)dataPtr);
      break;
    case audioMasterSetIcon: // Deprecated
      break;
    case audioMasterGetLanguage:
      result = kVstLangEnglish;
      break;
    case audioMasterOpenWindow: // Deprecated
      break;
    case audioMasterCloseWindow: // Deprecated
      break;
    case audioMasterGetDirectory:
      break;
    case audioMasterUpdateDisplay:
      break;
    case audioMasterBeginEdit:
      break;
    case audioMasterEndEdit:
      break;
    case audioMasterOpenFileSelector:
      break;
    case audioMasterCloseFileSelector:
      break;
    case audioMasterEditFile: // Deprecated
      break;
    case audioMasterGetChunkFile: // Deprecated
      break;
    case audioMasterGetInputSpeakerArrangement: // Deprecated
      break;
    default:
      break;
  }
  
  return result;
}

}

static float getSampleRate()
{
  return 44100.0f;
}

static int getBlockSize()
{
  return 512;
}

static void resumePlugin(AudioEffect *effect) 
{
  effect->dispatcher(effMainsChanged, 0, 1, NULL, 0.0f);
}

static void suspendPlugin(AudioEffect *effect)
{
  effect->dispatcher(effMainsChanged, 0, 0, NULL, 0.0f);
}

static void initPlugin(AudioEffect *effect)
{
  effect->dispatcher(effOpen, 0, 0, NULL, 0.0f);
  effect->dispatcher(effSetSampleRate, 0, 0, NULL, (float)getSampleRate());
  effect->dispatcher(effSetBlockSize, 0, getBlockSize(), NULL, 0.0f);
  resumePlugin(effect);
}

static void fillVstMidiEvent(const MidiEvent *midiEvent, VstMidiEvent* vstMidiEvent) 
{
  vstMidiEvent->type = kVstMidiType;
  vstMidiEvent->byteSize = sizeof(VstMidiEvent);
  vstMidiEvent->deltaFrames = midiEvent->deltaFrames;
  vstMidiEvent->midiData[0] = midiEvent->status;
  vstMidiEvent->midiData[1] = midiEvent->data1;
  vstMidiEvent->midiData[2] = midiEvent->data2;
  vstMidiEvent->flags = 0;
  vstMidiEvent->reserved1 = 0;
  vstMidiEvent->reserved2 = 0;
}

void processMidiEvents(AudioEffect *effect, MidiEvent *midiEvent, int numEvents) 
{
  VstEvents vstEvents;
  vstEvents.numEvents = numEvents;
  vstEvents.events[0] = (VstEvent*)malloc(sizeof(VstMidiEvent) * vstEvents.numEvents);
  
  for(int i=0; i<numEvents; i++) {
    VstMidiEvent* vstMidiEvent = (VstMidiEvent*)&(vstEvents.events[0][i]);
    fillVstMidiEvent(&midiEvent[i], vstMidiEvent);
  }

  effect->dispatcher(effProcessEvents, 0, 0, &vstEvents, 0.0f);
  free(vstEvents.events[0]);
}

void processReplacing(AudioEffect *effect, float **in, float **out, VstInt32 samples) 
{
  effect->processReplacing(in, out, samples);
}

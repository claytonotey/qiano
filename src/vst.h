#ifndef VST_H
#define VST_H

extern "C" {
VstIntPtr VSTCALLBACK vst2xPluginHostCallback(AEffect *effect, VstInt32 opcode, VstInt32 index, VstIntPtr value, void *dataPtr, float opt);
}


class MidiEvent {
public:
  void noteOn(int note, int velocity, int channel, int delta) {
    deltaFrames = delta;
    status = 0x90 | channel;
    data1 = note;
    data2 = velocity;
  }
  void noteOff(int note, int velocity, int channel, int delta) {
    deltaFrames = delta;
    status = 0x80 | channel;
    data1 = note;
    data2 = velocity;
  }

  char status;
  char data1;
  char data2;
  char meh;
  unsigned long deltaFrames;
};

void processMidiEvents(AudioEffect *effect, MidiEvent *midiEvent, int numEvents);
void processReplacing(AudioEffect *effect, float **in, float **out, VstInt32 samples);
#endif

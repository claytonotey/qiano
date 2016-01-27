/*---------------------------------------------------------------

   © 2001, Marcberg Soft & Hard GmbH, All Rights Reserved

---------------------------------------------------------------*/

#ifndef __VSTMIDI_H
#define __VSTMIDI_H

#include "audioeffectx.h"

//----------------------------------------------------------------------------- 
// constants & macros

const int PITCHBEND_MAX_RANGE = 36;
const int PITCHBEND_MAX = 16384;
const int PITCHBEND_NULL = 8192;

const int NUM_NOTES = 128;	// 128 midi notes
const double NOTE_UP_SCALAR = 1.059463094359295264561825294946;	// 12th root of 2
const long EVENTS_QUEUE_MAX = 30000;

//----------------------------------------------------------------------------- 
// enums

// these are the MIDI event status types that are handled
enum
{
	isNote,
	isPitchbend,
	isNotesOff,
	isAftertouch,
	isControl,
	isNull
};

//----------------------------------------------------------------------------- 
// types

// this holds MIDI event information
struct BlockEvents {
	int eventStatus;	// the event status MIDI byte
	int byte1;	// the first MIDI data byte
	int byte2;	// the second MIDI data byte
	long delta;	// the delta offset (the sample position in the current block where the event occurs)
	int channel;	// the MIDI channel
};


//----------------------------------------------------------------------------- 
// function prototypes

void reallyProcessEvents(VstEvents *events, BlockEvents *blockEvents, int *numBlockEvents);

#endif

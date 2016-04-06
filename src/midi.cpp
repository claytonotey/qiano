/*---------------------------------------------------------------
       Marc's VST MIDI stuff --- happened February 2001
---------------------------------------------------------------*/

#include "midi.h"
#include "vsteffect.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//-----------------------------------------------------------------------------------------
// this function gets called right before process() or processReplacing() 
// if VstEvents have occured during this processing block

void VstEffect :: reallyProcessEvents(VstEvents *events, BlockEvents *blockEvents, int *numBlockEvents)
{
  VstMidiEvent *midiEvent;
  int status;
  char *midiData;
  long i, j;
  bool sorted;
  BlockEvents tempEvent;


/* note:  This function depends on the base plugin class to zero numBlockEvents 
at the end of each processing block & does not do that itself because it is both 
prossible & allowable for processEvents() to be called more than once per block. */

  for (i = 0; (i < events->numEvents); i++) {
    // check to see if this event is MIDI; if no, then we try the for-loop again
    if ( ((events->events[i])->type) != kVstMidiType )
      continue;


    // cast the incoming event as a VstMidiEvent
    midiEvent = (VstMidiEvent*)events->events[i];


    // address the midiData[4] string from the event to this temp data pointer
    midiData = midiEvent->midiData;

    // save the channel number ...
    blockEvents[*numBlockEvents].channel = midiData[0] & 0x0f;
    // ... & then wipe out the channel (lower 4 bits) for simplicity
    status = midiData[0] & 0xf0;
    
    // looking at notes   (0x9* is Note On status ~ 0x8* is Note Off status)
    if ( (status == 0x90) || (status == 0x80) ) {
      // note-off received, set velocity to zero
      if (status == 0x80)
		blockEvents[*numBlockEvents].byte2 = 0;
      // use the received note velocity value
      else
		blockEvents[*numBlockEvents].byte2 = midiData[2] & 0x7f;
      

      blockEvents[*numBlockEvents].eventStatus = isNote;	// status
      blockEvents[*numBlockEvents].byte1 = midiData[1] & 0x7f;	// note
      blockEvents[*numBlockEvents].delta = midiEvent->deltaFrames;	// timing offset

      (*numBlockEvents)++;
    }

    // looking at pitchbend   (0xe* is Pitchbend status)
    else if (status == 0xe0) {
      blockEvents[*numBlockEvents].eventStatus = isPitchbend;	// status
      blockEvents[*numBlockEvents].byte1 = midiData[1] & 0x7f;	// LSB
      blockEvents[*numBlockEvents].byte2 = midiData[2] & 0x7f;	// MSB
      blockEvents[*numBlockEvents].delta = midiEvent->deltaFrames;	// timing offset
      
      (*numBlockEvents)++;
    }
    else if (status == 0xa0) {
      blockEvents[*numBlockEvents].eventStatus = isAftertouch;	// status
      blockEvents[*numBlockEvents].byte1 = midiData[1] & 0x7f;	// LSB
      blockEvents[*numBlockEvents].byte2 = midiData[2] & 0x7f;	// MSB
      blockEvents[*numBlockEvents].delta = midiEvent->deltaFrames;	// timing offset
     
      (*numBlockEvents)++;
    }

    // all notes off
    else if ( (status == 0xb0) && (midiData[1] == 0x7b) ) {
    blockEvents[*numBlockEvents].eventStatus = isNotesOff;	// status
      blockEvents[*numBlockEvents].delta = midiEvent->deltaFrames;	// timing offset
      
      (*numBlockEvents)++;
	} else if (status == 0xb0) {
      blockEvents[*numBlockEvents].eventStatus = isControl;	// status
      blockEvents[*numBlockEvents].byte1 = midiData[1] & 0x7f;	// LSB
      blockEvents[*numBlockEvents].byte2 = midiData[2] & 0x7f;	// MSB
      blockEvents[*numBlockEvents].delta = midiEvent->deltaFrames;	// timing offset
      
      (*numBlockEvents)++;
	} 


    midiEvent++;
		
    // don't go past the allocated space for the events queue
    if (*numBlockEvents >= EVENTS_QUEUE_MAX)
      *numBlockEvents = EVENTS_QUEUE_MAX;
  }

  // Sort the events in our queue so that their in chronological order.  (bubble sort)
  // The host is supposed to send them in order, but just in case...
  for (i=0; i < (*numBlockEvents-1); i++) {
    // default it to true & change it to false if the next loop finds unsorted items
    sorted = true;
    //
    for (j=0; j < (*numBlockEvents-1-i); j++)
      {
	// swap the neighbors if they're out of order
	if (blockEvents[j+1].delta < blockEvents[j].delta)
	  {
	    tempEvent = blockEvents[j];
	    blockEvents[j] = blockEvents[j+1];
	    blockEvents[j+1] = tempEvent;
	    sorted = false;
	  }
      }
    //
    // no need to go through all (numBlockEvents-1)! iterations
    // if the array is fully sorted already
    if (sorted)   break;
  }
}



/* I wrote this in an email to someone explaining how my MIDI handling works.  
   I figured it was worth throwing in here for anyone else who might look at my source code.
   It's a step by step explanation of what happens.

	In processEvents(), I receive a VstEvent array (which includes a 
counter value for how many items there are in that particular array) & 
then I look at every item.  First I check if it's a MIDI event (like in 
the SDK).  If it is, then I cast it into my VstMidiEvent variable (like in 
the SDK) & then start examining it.  I only want to deal with it if is one 
of 3 types of MIDI events:  a note, a pitchbend message, or a panic 
message.  That is what each of the next three "if"s are all about.  If the 
event is any one of those 3 things, then it gets worked on.
	You are right that the struct array blockEvents[] is my queue.  I fill 
up one item in that array for each interesting event that I receive during 
that processing block.  I store the status in my own way (using an enum 
that's in my header) as either isNote, isPitchbend, or isNotesOff.  This 
is what goes in the blockEvents.eventStatus field.  Then I take MIDI bytes 
1 & 2 & put them into blockEvents.byte1 & .byte2.  If it's a note off 
message, then I put 0 into byte2 (velocity).  By the way, with pitchbend, 
byte 1 is the LSB, but since LSB is, so far as I know, inconsistantly 
implemented by different MIDI devices, my plugin doesn't actually use that 
value in any way.  At the end of each of these "if"s, I update 
numBlockEvents because that is my counter that I look at during process() 
to see how many events I have to deal with during that block.
	& for each event, I store the deltaFrames value.  deltaFrames is the 
number of samples into that processing block that the event occurs.  This 
is what makes sample accurate timing (on MIDI data playback, not with live 
playing, of course) possible.  A VST host will send all of the upcoming 
events for a giving processing block to processEvents() & then the exact 
position within that processing block is given by deltaFrames.  If 
deltaFrames is 0, then the event occurs at the very beginning of the block.  
If deltaFrames is 333, then it occurs 333 samples into that processing 
block.  While it is not required, the SDK claims that hosts generally, as 
a matter of convention, send VstEvents to processEvents() in chronological 
order.  My plugin assumes that they are received in order.  There is no 
sorting done in my plugin.
	Now on to process().  Basically, I divide my process up into 
sub-chunks according to when events occur (which means according to the 
deltaFrames values).  I have two important variables here:   eventcount & 
currentBlockPosition.  The eventcount keeps track of how many of the 
events for that block I have addressed.  I initialize it to -1 so that 
first it will do processing up until the first event & then it will start 
counting events at 0 with the first event.  This is because there most 
likely will be audio to process before any events occur during that block 
(unless the block began with all notes off).  currentBlockPosition stores 
the sample start position of the current sub-chunk.  Basically it is the 
deltaFrames values of the latest event that I am working on.  It obviously 
starts out at 0.
	Next I start a "do" loop that cycles for every event.  First it 
evaluates the duration of the current processing sub-chunk, so it first 
checks to see if there are any more upcoming events in the queue.  If not, 
then the current chunk processes to the end of the processing block, if 
yes, then the current sub-chunk position is subtracted from the upcoming 
event's deltaFrames value.  I move up inputs & outputs accordingly, etc.
	Next comes a "for" loop that goes through my noteTable[] struct array 
& looks for any active notes (i.e. notes with a non-zero velocity value).  
All that I do during this loop is check the velocity of each note & then 
process the audio for that note if it has a non-zero velocity.
	After that "for" loop, I increment eventcount, leave the events loop 
if events are done, update currentBlockPosition, & then call heedEvents().  
heedEvents() is important.  heedEvents() is where I take in the effects of 
the next MIDI event.  Basically I tell heedEvents() which event number I 
am looking at & then it updates any vital stuff so that, when going 
through the next processing sub-chunk, all necessary changes have been 
made to make the next batch of processing take into account the impact of 
the latest event.  So heedEvents() is pretty much just a switch statement 
checking out which type of event is being analyzed & then implementing it 
in the appropriate way.  It would take quite a while to fully explain what 
happens in there & why, so I'll leave it up to you to determine whichever 
parts of it you want to figure out.  I don't think that it is very 
complicated (just many steps), but I could be wrong & if you can't figure 
out any parts of it, let me know & I'll explain.
	There is one more really important thing in process().  At the of the 
whole function, I reset numBlockEvents (the global events counter) to 0.  
This is extremely important because processEvents() only gets called if 
new events are received during a processing block.  If not, then 
processEvents() does not get called, numBlockEvents does not get zeroed 
at the beginning of processEvents(), & process() will process the same 
MIDI events over & over & over for every processing block until a new MIDI 
event is received.  This fact about processEvents() is not explained in 
the SDK & I spent FOREVER with a malfunctioning plugin until I figured 
this out.
*/

#!/bin/bash

ROOT=/Developer/SDKs/MacOSX10.4u.sdk

VSTROOT=../../vst2.4
export VSTSDK_DIR="$VSTROOT/sdk"
export VSTPLUG_DIR="$VSTROOT/plugin"
export VSTGUI_DIR="$VSTROOT/vstgui"
export VST_CFLAGS="-I$VSTSDK_DIR -I$VSTPLUG_DIR -I$VSTGUI_DIR"
export VST_LIBS="-framework Carbon -framework QuickTime -framework System -framework ApplicationServices"

#export CXXFLAGS="-arch i386"
#./configure --disable-universal_binary --enable-vst --enable-debug && make clean && make
./configure --enable-vst --enable-debug && make clean && make

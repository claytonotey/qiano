AUTOMAKE_OPTIONS = subdir-objects

CXX = clang++
AM_CXXFLAGS = @QIANO_CFLAGS@ @VST_CFLAGS@ -I../include -flax-vector-conversions -march=native -framework Accelerate
#-fpascal-strings
libpauknincludedir = $(includedir)
libpaukninclude_HEADERS = qiano.h

lib_LTLIBRARIES = libqiano.la
libqiano_la_SOURCES = dwgs.cpp filter.cpp hammer.cpp reverb.cpp utils.cpp

qianovst_LDADD = libqiano.la @VST_LIBS@
qianovst_LDFLAGS = $(os_ldflags) -bundle

qiano_LDADD = libqiano.la @VST_LIBS@
qiano_LDFLAGS = $(os_ldflags)

filtertest_LDADD = libqiano.la
filtertest_LDFLAGS = $(os_ldflags)

shared_SOURCES = vst.cpp midi.cpp vsteffect.cpp @VSTSDK_DIR@/audioeffect.cpp @VSTSDK_DIR@/audioeffect.h @VSTSDK_DIR@/audioeffectx.cpp @VSTSDK_DIR@/audioeffectx.h @VSTSDK_DIR@/aeffeditor.h @VSTPLUG_DIR@/aeffect.h @VSTPLUG_DIR@/aeffectx.h @VSTSDK_DIR@/vstplugmain.cpp FFTConvolver.cpp Utilities.cpp TwoStageFFTConvolver.cpp AudioFFT.cpp

qiano_SOURCES = qiano.cpp $(shared_SOURCES)

qianovst_SOURCES = qiano.cpp $(shared_SOURCES)

filtertest_SOURCES = filtertest.cpp AudioFFT.cpp

bin_PROGRAMS = qiano qianovst filtertest

clean: cleanapproot clean-am

cleanapproot:
	rm -rf $(top_builddir)/approot

all: $(top_builddir)/approot/Library/Audio/Plug-Ins/VST/qiano.vst $(top_builddir)/mac/vst/Info.plist $(top_builddir)/mac/vst/PkgInfo

$(top_builddir)/approot/Library/Audio/Plug-Ins/VST/qiano.vst: qianovst
	mkdir -p $(top_builddir)/approot/Library/Audio/Plug-Ins/VST/qiano.vst/Contents/MacOS
	cp -f qianovst $(top_builddir)/approot/Library/Audio/Plug-Ins/VST/qiano.vst/Contents/MacOS/qiano
	cp -f $(top_builddir)/mac/vst/Info.plist $(top_builddir)/approot/Library/Audio/Plug-Ins/VST/qiano.vst/Contents
	cp -f $(top_builddir)/mac/vst/PkgInfo $(top_builddir)/approot/Library/Audio/Plug-Ins/VST/qiano.vst/Contents
	mkdir -p $(top_builddir)/approot/Library/Audio/Plug-Ins/VST/qiano.vst/Contents/Resources
#	sudo chown -R root:admin $(top_builddir)/approot/Library

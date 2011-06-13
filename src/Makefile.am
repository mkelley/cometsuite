# We'll need to know the location of the installed FITS composition files
AM_CPPFLAGS = -D_DATADIR='"$(pkgdatadir)"'

bin_PROGRAMS = rundynamics integratelist integratexyz syn2ds9 xyz2fits	\
	xyz2psd xyzinfo
noinst_PROGRAMS = accutest jettest rbvtest vtest randtest

accutest_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.h Vej.cc accutest.cc longlat.h mskerr.h	\
	paramSet.h paramSet.cc particle.h particle.cc planet_lookup.c	\
	rundynamics.h state.h StringConv.h StringConv.cc xyzstream.h	\
	xyzstream.cc

integratexyz_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Integrator.h		\
	Integrator.cc Physical.h Physical.cc Vector.h Vector.cc Vej.h	\
	Vej.cc everhart.c getopt.h getopt.c getopt1.c integratexyz.cc	\
	logFile.h logFile.cc longlat.h mskerr.h paramSet.h paramSet.cc	\
	particle.h particle.cc pfunctions.h pfunctions.cc		\
	planet_lookup.c ra15.h ra15.cc rundynamics.h state.h		\
	StringConv.h StringConv.cc xyzstream.h xyzstream.cc

integratelist_SOURCES = Composition.h CoordTrans.h CoordTrans.cc	\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Integrator.h		\
	Integrator.cc Physical.h Physical.cc Vector.h Vector.cc Vej.h	\
	Vej.cc everhart.c getopt.h getopt.c getopt1.c integratelist.cc	\
	logFile.h logFile.cc longlat.h mskerr.h paramSet.h paramSet.cc	\
	particle.h particle.cc pfunctions.h pfunctions.cc		\
	planet_lookup.c ra15.h ra15.cc rundynamics.h state.h		\
	StringConv.h StringConv.cc xyzstream.h xyzstream.cc

jettest_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.h Vej.cc longlat.h jettest.cc		\
	paramSet.h paramSet.cc particle.h particle.cc pfunctions.h	\
	pfunctions.cc rundynamics.h state.h StringConv.h StringConv.cc

randtest_SOURCES = Distribution.h Distribution.cc randtest.cc

rbvtest_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.cc longlat.h paramSet.h paramSet.cc	\
	particle.h particle.cc pfunctions.h pfunctions.cc rbvtest.cc	\
	rundynamics.h state.h StringConv.h StringConv.cc

rundynamics_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Integrator.h		\
	Integrator.cc Physical.h Physical.cc Vector.h Vector.cc Vej.h	\
	Vej.cc everhart.c getopt.h getopt.c getopt1.c logFile.h		\
	logFile.cc longlat.h mskerr.h paramSet.h paramSet.cc		\
	particle.h particle.cc pfunctions.h pfunctions.cc		\
	planet_lookup.c ra15.h ra15.cc rundynamics.h rundynamics.cc	\
	state.h StringConv.h StringConv.cc xyzstream.h xyzstream.cc

syn2ds9_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.h Vej.cc getopt.h getopt.c getopt1.c	\
	longlat.h mskerr.h paramSet.h paramSet.cc particle.h		\
	particle.cc planet_lookup.c projection.h projection.cc		\
	rundynamics.h state.h syn2ds9.cc StringConv.h StringConv.cc	\
	xyzstream.h xyzstream.cc

vtest_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.h Vej.cc longlat.h paramSet.h		\
	paramSet.cc particle.h particle.cc pfunctions.h pfunctions.cc	\
	rundynamics.h state.h StringConv.h StringConv.cc vtest.cc

xyz2fits_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.h Vej.cc everhart.c getopt.h getopt.c	\
	getopt1.c longlat.h mskerr.h paramSet.h paramSet.cc particle.h	\
	particle.cc particleWeights.h particleWeights.cc pfunctions.h	\
	pfunctions.cc planet_lookup.c projection.h projection.cc	\
	round.h rundynamics.h state.h StringConv.h StringConv.cc	\
	xyz2fits.cc xyzImage.h xyzImage.cc xyzProject.h xyzProject.cc	\
	xyzstream.h xyzstream.cc

xyz2psd_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.h Vej.cc getopt.h getopt.c getopt1.c	\
	longlat.h mskerr.h paramSet.h paramSet.cc particle.h		\
	particle.cc particleWeights.h particleWeights.cc pfunctions.h	\
	pfunctions.cc planet_lookup.c projection.h projection.cc	\
	rundynamics.h state.h StringConv.h StringConv.cc xyz2psd.cc	\
	xyzHist.h xyzHist.cc xyzProject.h xyzProject.cc xyzstream.h	\
	xyzstream.cc

xyzinfo_SOURCES = Composition.h CoordTrans.h CoordTrans.cc		\
	Distribution.h Distribution.cc Dynamical.h Dynamical.cc		\
	GrainTemperature.h GrainTemperature.cc Physical.h Physical.cc	\
	Vector.h Vector.cc Vej.h Vej.cc longlat.h mskerr.h paramSet.h	\
	paramSet.cc particle.h particle.cc pfunctions.h pfunctions.cc	\
	planet_lookup.c rundynamics.h state.h StringConv.h		\
	StringConv.cc xyzinfo.cc xyzstream.h xyzstream.cc

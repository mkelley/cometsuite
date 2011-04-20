/***************************************************************************

  Observes an xyzfile and returns the grain distribution.

  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <valarray>
#include <getopt.h>
#include "xyzHist.h"
#include "rundynamics.h"
#include "StringConv.h"

#define SUBPROJECT "xyz2psd"

using namespace std;

struct runtimePar { // runtime parameters
  bool verbose;
};

bool parseCommandLine(int, char**, runtimePar&, xyzHist&);
void usage();

/******************************************************************************/
int main(int argc, char *argv[])
{
  runtimePar runtime = { false };
  xyzHist bin;

  // default is no thermal weighting and output to stdout
  bin.wavelength(0);
  bin.outfileName("/dev/stdout");

  if (parseCommandLine(argc, argv, runtime, bin)) {
    return EXIT_FAILURE;
  }

  bin.createPSD();
  bin.writePSD();

  return EXIT_SUCCESS;
}

/******************************************************************************/
/** Checks the command line for proper input and loads a parameter
    file if possible. */
bool parseCommandLine(int argc, char** argv, runtimePar& runtime,
		      xyzHist& bin) {
  bool files = false, help = false;

  if (argc > 1) {
    int c;
    int digit_optind = 0;
    option longOptions[] = {
      {"afhro",      required_argument, 0, 'a'},
      {"aper",       required_argument, 0, 0},
      {"ageinvert",  no_argument,       0, 0},
      {"agerange",   required_argument, 0, 0},
      {"betainvert", no_argument,       0, 0},
      {"betarange",  required_argument, 0, 0},
      {"bin",        required_argument, 0, 0},
      {"binsize",    required_argument, 0, 'b'},
      {"fscale",     required_argument, 0, 0},
      {"help",       no_argument,       0, 'h'},
      {"latinvert",  no_argument,       0, 0},
      {"latrange",   required_argument, 0, 0},
      {"max",        required_argument, 0, 'm'},
      {"observer",   required_argument, 0, 0},
      {"platescale", required_argument, 0, 'p'},
      {"psd",        required_argument, 0, 0},
      {"rhlimit",    required_argument, 0, 0},
      {"scattering", required_argument, 0, 0},
      {"size",       required_argument, 0, 's'},
      {"suninvert",  no_argument,       0, 0},
      {"sunrange",   required_argument, 0, 0},
      {"thermal",    required_argument, 0, 't'},
      {"verbose",    no_argument,       0, 'v'},
      {0, 0, 0, 0}
    };

    while (1) {
      int optionIndex = 0;
      int this_option_optind = optind ? optind : 1;

      c = getopt_long(argc, argv, "a:b:d:hm:o:p:s:t:v", longOptions,
		      &optionIndex);
      if (c == -1) break;

      switch(c) {
      case 0:
	if (longOptions[optionIndex].name == "aper") {
	  bin.aperture(StringConv(optarg).toScalar<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "ageinvert") {
	  bin.ageInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "agerange") {
	  bin.ageRange(StringConv(optarg).toValarray<float>() *= 86400);
	  break;
	}
	if (longOptions[optionIndex].name == "betainvert") {
	  bin.betaInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "betarange") {
	  bin.betaRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "bin") {
	  bin.bin(string(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "fscale") {
	  bin.fileScales(StringConv(optarg).toVector<double>());
	  break;
	}
	if (longOptions[optionIndex].name == "latinvert") {
	  bin.latInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "latrange") {
	  bin.latRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "npole") {
	  bin.npole(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "observer") {
	  bin.observerName(string(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "psd") {
	  bin.nuclearPsd(string(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "rhlimit") {
	  bin.rhLimit(atof(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "scattering") {
	  float w = StringConv(optarg).toScalar<float>();
	  if (w <= 0) {
	    bin.scatteringMode(false);
	  } else {
	    bin.wavelength(w);
	    bin.scatteringMode(true);
	  }
	  break;
	}
	if (longOptions[optionIndex].name == "suninvert") {
	  bin.sunInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "sunrange") {
	  bin.sunRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	break;

      case 'a': bin.afrhoSlope(StringConv(optarg).toScalar<float>()); break;
      case 'b': bin.binSize(StringConv(optarg).toScalar<float>()); break;
      case 'd': bin.graindensity(atof(optarg)); break;
      case 'h': help = true; break;
      case 'm': bin.max(StringConv(optarg).toScalar<long>()); break;
      case 'o': bin.outfileName(string(optarg)); break;
      case 't':
	{
	  float w = StringConv(optarg).toScalar<float>();
	  if (w <= 0) {
	    bin.thermalMode(false);
	  } else {
	    bin.wavelength(w);
	    bin.thermalMode(true);
	  }
	}
	break;
      case 'v':
	runtime.verbose = true;
	bin.verbose(true);
	break;
      case '?':
      default: break;
      }
    }

    while (optind < argc) {
      // add an input file
      bin.xyzfileNames(string(argv[optind]));
      files = true;
      optind++;
    }

    if (!files || help) {
      usage();
      return true;
    }

  } else {
    usage();
    return true;
  }

  return false;
}

/******************************************************************************/
/** Prints the program usage and info. */
void usage() {
  cerr << PACKAGE_STRING << "\n";
  cerr << "Usage: " << SUBPROJECT << " [OPTION]... xyzfile1.xyz [xyzfile2.xyz...]\n\
\n\
Mandatory arguments to long options are mandatory for short options too.\n\
    -a, --afrho=VAL         weight the final size distribution by the dust\n\
                            production, afrho = r^afrhoSlope (default: -2)\n\
    --aper=VAL              set the diameter of an aperture to VAL size, set to\n\
                            -1 to include all particles, regardless of distance\n\
                            to the nucleus (units: arcsec, default: 3)\n\
    --ageinvert             use all particles except RANGE[0] < age < RANGE[1]\n\
                            (default: disabled)\n\
    --agerange=RANGE        only use particles with ages inside RANGE:\n\
                            RANGE[0] <= age <= RANGE[1] (units: days, default:\n\
                            use all)\n\
    --betainvert            use all particles except RANGE[0] < beta < RANGE[1]\n\
                            (default: disabled)\n\
    --betarange=RANGE       only use particles with beta values inside RANGE:\n\
                            betaLowerLimit <= beta <= betaUpperLimit (default:\n\
                            use all)\n\
    --bin=STRING            set to size (default), beta, or age to return each\n\
                            respective particle distribution\n\
    -b, --binsize=VAL       set the bin size to VAL, or -1 (default) for no\n\
                            binning; if the particle radius is not defined in\n\
                            the input file, then radius = 0.285 / beta is\n\
                            assumed (units: microns for size, unitless for\n\
                            beta, and days for age)\n\
    -d, --density=VAL       use this grain density when transforming between beta\n\
                            and radius; ignored if radius or grain density is\n\
                            defined in the input file (units: g/cm^3, default: 1)\n\
    --fscale=s1,s2,...      a list of values with which to scale each input file\n\
                            s1 corresponds to the first file; if there are more\n\
                            files than scales, the last scale will be repeated\n\
                            (default: 1)\n\
    -h, --help              display this help\n\
    --latinvert             plot all particles except RANGE[0] < lat < RANGE[1]\n\
                            (default: disabled)\n\
    --latrange=RANGE        only plot particles with latitude values inside RANGE:\n\
                            latLowerLimit <= beta <= latUpperLimit (default:\n\
                            plot all)\n\
    -m, --max=VAL           use the first VAL particles (default: use all)\n\
    --npole=VAL             direction of the north pole in ecliptic coordinates\n\
    -o FILENAME             output to file FILENAME (default: stdout)\n\
    --observer=STRING       the location of the observer, Earth (default),\n\
                            Spitzer, or a position vector in units of AU\n\
    --psd=STRING            set the nuclear particle size distribution\n\
                            function, current options are ism, a^x (where x\n\
                            is some value), hanner a0 M N, or none; all psds\n\
                            assume the simulation psd is dn/dlog(a) ~ 1\n\
                            (default: none)\n\
    --rhlimit=VAL           remove grains produced at rh > rhlimit; set to -1 for\n\
                            no limit (default: -1)\n\
    --scattering=VAL        weight the final PSD by each particle's scattered\n\
                            light at wavelength VAL (units: microns,default: 0)\n\
    -t, --thermal=VAL       weight the particle distributions by each the\n\
                            thermal emission at wavelength VAL, set to zero for\n\
                            no thermal weighting (units: microns, default: 0)\n\
    -v, --verbose           output more info than necessary\n\
\n\
(c) 2005-2010 Michael S. Kelley\n";
  return;
}

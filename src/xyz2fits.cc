/***************************************************************************

  Converts xyzfiles to FITS files.

  Copyright (C) 2005-2010,2012 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <CCfits/CCfits>
#include "xyzImage.h"
#include "rundynamics.h"
#include "StringConv.h"

#define SUBPROJECT "xyz2fits"

using namespace std;

struct runtimePar { // runtime parameters
  bool verbose;
};

bool parseCommandLine(int, char**, runtimePar&, xyzImage&);
void usage();
void printParameters(xyzImage&);
void writeParameters(xyzImage&);

/******************************************************************************/
int main(int argc, char *argv[])
{
  runtimePar runtime = { false };

  xyzImage image;

  if (parseCommandLine(argc, argv, runtime, image)) {
    return EXIT_FAILURE;
  }

  printParameters(image);

  try {
    image.createImages();
    image.writeImages();
    writeParameters(image);
  } catch (xyzImageFlags flag) {
    if (flag == cantCreateFits) {
      cerr << "Error creating fits image.\n";
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

/******************************************************************************/
/** Checks the command line for proper input and loads a parameter
    file if possible. */
bool parseCommandLine(int argc, char** argv, runtimePar& runtime,
		      xyzImage& image) {
  bool files = false, help = false;

  if (argc > 1) {
    int c;
    int digit_optind = 0;
    option longOptions[] = {
      {"afrho",      required_argument, 0, 'a'},
      {"ageinvert",  no_argument,       0, 0  },
      {"agerange",   required_argument, 0, 0  },
      {"betainvert", no_argument,       0, 0  },
      {"betarange",  required_argument, 0, 0  },
      {"density",    required_argument, 0, 'd'},
      {"ecliptic",   no_argument,       0, 0  },
      {"fscales",    required_argument, 0, 0  },
      {"help",       no_argument,       0, 'h'},
      {"jet",        required_argument, 0, 'j'},
      {"latinvert",  no_argument,       0, 0  },
      {"latrange",   required_argument, 0, 0  },
      {"loninvert",  no_argument,       0, 0  },
      {"lonrange",   required_argument, 0, 0  },
      {"max",        required_argument, 0, 'm'},
      {"npole",      required_argument, 0, 0  },
      {"observer",   required_argument, 0, 0  },
      {"offset",     required_argument, 0, 0  },
      {"period",     required_argument, 0, 0  },
      {"phase",      required_argument, 0, 0  },
      {"platescale", required_argument, 0, 'p'},
      {"porosity",   required_argument, 0, 0  },
      {"psd",        required_argument, 0, 0  },
      {"radinvert",  no_argument,       0, 0  },
      {"radrange",   required_argument, 0, 0  },
      {"rhlimit",    required_argument, 0, 0  },
      {"scattering", required_argument, 0, 0  },
      {"size",       required_argument, 0, 's'},
      {"suninvert",  no_argument,       0, 0  },
      {"sunrange",   required_argument, 0, 0  },
      {"thermal",    required_argument, 0, 't'},
      {"verbose",    no_argument,       0, 'v'},
      {"vlimit",     required_argument, 0, 0  },
      {0, 0, 0, 0}
    };

    while (1) {
      int optionIndex = 0;
      int this_option_optind = optind ? optind : 1;

      c = getopt_long(argc, argv, "a:d:hj:m:o:p:s:t:v", longOptions,
		      &optionIndex);
      if (c == -1) break;

      switch(c) {
      case 0:
	if (longOptions[optionIndex].name == "ageinvert") {
	  image.ageInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "agerange") {
	  image.ageRange(StringConv(optarg).toValarray<float>() *= 86400);
	  break;
	}
	if (longOptions[optionIndex].name == "betainvert") {
	  image.betaInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "betarange") {
	  image.betaRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "ecliptic") {
	  image.ecliptic(true);
	  break;
	}
	if (longOptions[optionIndex].name == "fscales") {
	  image.fileScales(StringConv(optarg).toVector<double>());
	  break;
	}
	if (longOptions[optionIndex].name == "latinvert") {
	  image.latInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "latrange") {
	  image.latRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "loninvert") {
	  image.lonInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "lonrange") {
	  image.lonRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "npole") {
	  image.npole(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "observer") {
	  image.observerName(string(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "offset") {
	  image.offset(StringConv(optarg).toValarray<float>() /= 3600);
	  break;
	}
	if (longOptions[optionIndex].name == "period") {
	  image.rotPeriod(atof(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "phase") {
	  image.rotPhase(atof(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "psd") {
	  image.nuclearPsd(string(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "radinvert") {
	  image.radInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "radrange") {
	  image.radRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "rhlimit") {
	  image.rhLimit(atof(optarg));
	  break;
	}
	if (longOptions[optionIndex].name == "scattering") {
	  float w = atof(optarg);
	  if (w <= 0) {
	    image.scatteringMode(false);
	  } else {
	    image.wavelength(w);
	    image.scatteringMode(true);
	  }
	  break;
	}
	if (longOptions[optionIndex].name == "suninvert") {
	  image.sunInvert(true);
	  break;
	}
	if (longOptions[optionIndex].name == "sunrange") {
	  image.sunRange(StringConv(optarg).toValarray<float>());
	  break;
	}
	if (longOptions[optionIndex].name == "vlimit") {
	  image.vLimit(atof(optarg));
	  break;
	}
	break; // stop processing long options

      case 'a': image.afrhoSlope(atof(optarg)); break;
      case 'd': image.graindensity(atof(optarg)); break;
      case 'h': help = true; break;
      case 'j':
	{
	  valarray<float> jet = StringConv(optarg).toValarray<float>();
	  if (jet.size() != 3) {
	    cerr << "Need 3 parameters for --jet, but got " << jet.size()
		 << ".\n";
	    return true;
	  }
	  image.setJet(true);
	  image.jet(jet[slice(0, 2, 1)]);
	  image.jetHalfAngle(jet[2] / 2.0);
	  image.rotPeriod(jet[3]);
	}
      case 'm': image.max(atol(optarg)); break;
      case 'o': image.outfileName(string(optarg)); break;
      case 'p': image.platescale(StringConv(optarg).toValarray<float>()); break;
      case 's': image.size(StringConv(optarg).toValarray<long>()); break;
      case 't':
	{
	  float w = atof(optarg);
	  if (w <= 0) {
	    image.thermalMode(false);
	  } else {
	    image.wavelength(w);
	    image.thermalMode(true);
	  }
	}
	break;
      case 'v':
	runtime.verbose = true;
	image.verbose(true);
	break;
      case '?':
      default: break;
      }
    }

    while (optind < argc) {
      // add an input file
      image.xyzfileNames(string(argv[optind]));
      files = true;
      optind++;
    }

    if (!files || help) {
      usage();
      return true;
    }

    return false;  // no errors
  } else {
    usage();
    return true;
  }
}

/******************************************************************************/
/** Prints the program usage and info. */
void usage() {
  cerr << PACKAGE_STRING << "\n";
  cerr << "Usage: " << SUBPROJECT << " [OPTION]... xyzfile1.xyz [xyzfile2.xyz...]\n\
\n\
Mandatory arguments to long options are mandatory for short options too.\n\
\n\
  Generic options:\n\
\n\
    -h, --help              display this help\n\
    -v, --verbose           output more info than necessary\n\
\n\
  Output units/format:\n\
\n\
    --ecliptic              output in ecliptic coordinates\n\
    -o FILENAME             output to file FILENAME (default: out.fits)\n\
    --offset=VAL            Coordinate offsets of the comet position in\n\
                            arcseconds, only affects the WCS header keywords\n\
                            (i.e., this does not include the cos(Dec)\n\
                            correction)\n\
    -p, --platescale=VAL    set the platescale to VAL, may be two values\n\
                            for x and y platescales, or one value, z, in\n\
                            which case the platescale will be -z z (units:\n\
                            arcsec/pix, default: -1 1)\n\
    -s, --size=VAL          set the image size to VAL (units: pixels,\n\
                            default: 512 512)\n\
\n\
  Grain weighting:\n\
\n\
    -a, --afrho=VAL         weight the final image by dust production,\n\
                            afrho = r^afrhoSlope (default: -2)\n\
    --fscales=s1,s2,...     a list of values with which to scale each input\n\
                            file s1 corresponds to the first file; if there\n\
                            are more files than scales, the last scale will\n\
                            be repeated (default: 1)\n\
    --psd=STRING            set the nuclear particle size distribution\n\
                            function, current options are ism, a^x (where\n\
                            x is some value), hanner a0 M N, or none; all\n\
                            psds assume the simulation psd is dn/dlog(a) ~ 1;\n\
                            For reference, dn/dbeta = a^x = beta^(-x-2)\n\
                            (default: none)\n\
    --scattering=VAL        weight the final image by each particle's\n\
                            scattered light at wavelength VAL (units:\n\
                            microns, default: 0)\n\
    -t, --thermal=VAL       weight the final image by each particle's thermal\n\
                            emission at wavelength VAL, set to zero for no\n\
                            thermal weighting (units: microns, default: 24)\n\
\n\
  Modify grain/comet/observer parameters:\n\
\n\
    -d, --density=VAL       use this grain density when transforming\n\
                            between beta and radius; ignored if radius or\n\
                            grain density is defined in the input file\n\
                            (units: g/cm^3, default: 1)\n\
    --npole=VAL             direction of the north pole in ecliptic\n\
                            coordinates\n\
    --observer=STRING       the location of the observer, Earth (default),\n\
                            Spitzer, or a position vector in units of AU\n\
    --period=VAL            set the rotation period in hours (default: 0)\n\
    --phase=VAL             set the rotation phase (at the time of \n\
                            observation) in degrees, where phase = 0 is along\n\
                            the Vernal Equinox. (default: 0)\n\
\n\
  Particle limiting options:\n\
\n\
    --ageinvert             plot all particles except\n\
                            RANGE[0] < age < RANGE[1] (default: disabled)\n\
    --agerange=RANGE        only plot particles with ages inside RANGE:\n\
                            RANGE[0] <= age <= RANGE[1] (units: days,\n\
                            default: plot all)\n\
    --betainvert            plot all particles except\n\
                            RANGE[0] < beta < RANGE[1] (default: disabled)\n\
    --betarange=RANGE       only plot particles with beta values inside\n\
                            RANGE: betaLowerLimit <= beta <= betaUpperLimit\n\
                            (default: plot all)\n\
    --jet=lon,lat,th,per    specify jet location (longitude, latitude),\n\
                            opening angle, and period. (units: deg, hr)\n\
    --latinvert             plot all particles except\n\
                            RANGE[0] < lat < RANGE[1] (default: disabled)\n\
    --latrange=RANGE        only plot particles with latitude values inside\n\
                            RANGE: latLowerLimit <= lat <= latUpperLimit\n\
                            (default: plot all)\n\
    --loninvert             plot all particles except\n\
                            RANGE[0] < lon < RANGE[1] (default: disabled)\n\
    --lonrange=RANGE        only plot particles with longitude values inside\n\
                            RANGE: lonLowerLimit <= lon <= lonUpperLimit\n\
                            where longitudes are defined from 0 to 360 deg\n\
                            (default: plot all)\n\
    -m, --max=VAL           plot the first VAL particles (default: plot\n\
                            all)\n\
    --radinvert             plot all particles except\n\
                            RANGE[0] < radius < RANGE[1] (default: disabled)\n\
    --radrange=RANGE        only plot particles with radii inside\n\
                            RANGE: RANGE[0] <= radius <= RANGE[1]\n\
                            (default: plot all)\n\
    --rhlimit=VAL           remove grains produced at rh > rhlimit; set to\n\
                            -1 for no limit (default: -1)\n\
    --suninvert             plot all particles except\n\
                            RANGE[0] < z_sun < RANGE[1] where z_sun is the\n\
                            Sun-zenith angle (default: disabled)\n\
    --sunrange=RANGE        only plot particles with latitude values inside\n\
                            RANGE: RANGE[0] <= z_sun <= RANGE[1], where\n\
                            z_sun is the Sun-zenith angle (default: plot all)\n\
    --vlimit=VAL            limit the velocities to VAL*sqrt(beta/r_h)\n\
\n"
       << "(c) 2005-2010,2012 Michael S. Kelley\n";
  return;
}

/** Write parameters to the output file or screen. */
void printParameters(xyzImage& image) {
  cout << "\nafrhoSlope\t" << image.afrhoSlope() << "\n";
  cout << "ageInvert\t" << image.ageInvert() << "\n";
  if (image.ageRange().size()) {
    cout << "ageRange\t" << image.ageRange()[0] / 86400 << " "
	 << image.ageRange()[1] / 86400 << "\n";
  } else {
    cout << "ageRange\t" << "none" << "\n";
  }
  cout << "betaInvert\t" << image.betaInvert() << "\n";
  if (image.betaRange().size()) {
    cout << "betaRange\t" << image.betaRange()[0] << " " <<
      image.betaRange()[1] << "\n";
  } else {
    cout << "betaRange\t" << "none" << "\n";
  }
  cout << "radInvert\t" << image.radInvert() << "\n";
  if (image.radRange().size()) {
    cout << "radRange\t" << image.radRange()[0] << " " <<
      image.radRange()[1] << "\n";
  } else {
    cout << "radRange\t" << "none" << "\n";
  }
  cout << "density\t\t" << image.graindensity() << "\n";
  cout << "ecliptic\t" << image.ecliptic() << "\n";
  cout << "fscales\t\t";
  for (int i=0; i<image.fileScales().size(); i++)
    cout << image.fileScales()[i] << " ";
  cout << "\n";
  cout << "latInvert\t" << image.latInvert() << "\n";
  if (image.latRange().size()) {
    cout << "latRange\t" << image.latRange()[0] << " " <<
      image.latRange()[1] << "\n";
  } else {
    cout << "latRange\t" << "none" << "\n";
  }
  cout << "lonInvert\t" << image.lonInvert() << "\n";
  if (image.lonRange().size()) {
    cout << "lonRange\t" << image.lonRange()[0] << " " <<
      image.lonRange()[1] << "\n";
  } else {
    cout << "lonRange\t" << "none" << "\n";
  }
  cout << "max\t\t" << image.max() << "\n";
  cout << "output file\t" << image.outfileName() << "\n";
  cout << "observer\t" << image.observerName() << "\n";
  cout << "WCS offset\t" << image.offset().lambda * 3600 << " " <<
    image.offset().beta * 3600 << "\n";
  cout << "period\t\t" << image.rotPeriod() << "\n";
  cout << "phase\t\t" << image.rotPhase() << "\n";
  cout << "platescale\t" << image.platescale()[0] << " " <<
    image.platescale()[1] << "\n";
  cout << "PSD\t\t" << image.nuclearPsd() << "\n";
  cout << "image size\t" << image.size()[0] << " " << image.size()[1] << "\n";
  cout << "scatteringMode\t" << image.scatteringMode() << "\n";
  cout << "thermalMode\t" << image.thermalMode() << "\n";
  cout << "wavelength\t" << image.wavelength() << "\n";
  cout << "north pole\t" << image.npole()[0] << " " << image.npole()[1] << "\n";
  cout << "rhLimit\t\t" << image.rhLimit() << "\n";
  cout << "sunInvert\t" << image.sunInvert() << "\n";
  if (image.sunRange().size()) {
    cout << "sunRange\t" << image.sunRange()[0] << " " <<
      image.sunRange()[1] << "\n";
  } else {
    cout << "sunRange\t" << "none" << "\n";
  }
  cout << "vLimit\t\t" << image.vLimit() << "\n";
  cout << endl;
}

/** Write parameters to the output file or screen. */
void writeParameters(xyzImage& image) {
  image.addFitsComment("--------------------------------------------");
  image.addFitsComment(string("rundynamics runtime parameters"));
  image.addFitsComment("--------------------------------------------");
  image.addFitsKeyword("AFRHOSLP", image.afrhoSlope(), "Q_d \\propto r^afrhoSlope");
  image.addFitsKeyword("AGEINVRT", image.ageInvert(), "Invert the logic on agerange");
  if (image.ageRange().size()) {
    image.addFitsKeyword("AGERNG1", image.ageRange()[0] / 86400, "Age lower limit (days)");
    image.addFitsKeyword("AGERNG2", image.ageRange()[1] / 86400, "Age upper limit (days)");
  } else {
    image.addFitsKeyword("AGERNG1", "none", "Age lower limit");
    image.addFitsKeyword("AGERNG2", "none", "Age upper limit");
  }
  image.addFitsKeyword("BETINVRT", image.betaInvert(), "Invert the logic on betarange");
  if (image.betaRange().size()) {
    image.addFitsKeyword("BETRNG1", image.betaRange()[0], "Beta lower limit");
    image.addFitsKeyword("BETRNG2", image.betaRange()[1], "Beta upper limit");
  } else {
    image.addFitsKeyword("BETRNG1", "none", "Beta lower limit");
    image.addFitsKeyword("BETRNG2", "none", "Beta upper limit");
  }
  image.addFitsKeyword("DENSITY", image.graindensity(), "Grain density (g/cm^3)");
  image.addFitsKeyword("ECLIPTIC", image.ecliptic(), "Output in ecliptic coordinates");
  image.addFitsKeyword("LATINVRT", image.latInvert(), "Invert the logic on latrange");
  if (image.latRange().size()) {
    image.addFitsKeyword("LATRNG1", image.latRange()[0], "Latitude lower limit");
    image.addFitsKeyword("LATRNG2", image.latRange()[1], "Latitude upper limit");
  } else {
    image.addFitsKeyword("LATRNG1", "none", "Latitude lower limit");
    image.addFitsKeyword("LATRNG2", "none", "Latitude upper limit");
  }
  image.addFitsKeyword("MAX", image.max(), "Maximum number of particles to output");
  image.addFitsKeyword("OUTFILE", image.outfileName(), "Output file name");
  image.addFitsKeyword("OBSERVER", image.observerName(), "Observer name or coordinates (AU)");
  image.addFitsKeyword("OFFSET1", image.offset().lambda * 3600, "WCS axis 1 offset (arcsec)");
  image.addFitsKeyword("OFFSET2", image.offset().beta * 3600, "WCS axis 2 offset (arcsec)");
  image.addFitsKeyword("PERIOD", image.rotPeriod(), "Nucleus rotation period (hr)");
  image.addFitsKeyword("PHASE", image.rotPhase(), "Nucleus phase at time of observation (deg)");
  image.addFitsKeyword("PLTSCL1", image.platescale()[0], "Image axis 1 platescale (arcsec/pixel)");
  image.addFitsKeyword("PLTSCL2", image.platescale()[1], "Image axis 2 platescale (arcsec/pixel)");
  image.addFitsKeyword("PSD", image.nuclearPsd(), "Ejected particle size distribution");
  image.addFitsKeyword("SIZE1", image.size()[0], "Image axis 1 length");
  image.addFitsKeyword("SIZE2", image.size()[1], "Image axis 2 length");
  image.addFitsKeyword("SUNINVRT", image.sunInvert(), "Invert the logic on sunrange");
  if (image.sunRange().size()) {
    image.addFitsKeyword("SUNRNG1", image.sunRange()[0], "Sun-zenith angle lower limit");
    image.addFitsKeyword("SUNRNG2", image.sunRange()[1], "Sun-zenith angle upper limit");
  } else {
    image.addFitsKeyword("SUNRNG1", "none", "Sun-zenith angle lower limit");
    image.addFitsKeyword("SUNRNG2", "none", "Sun-zenith angle upper limit");
  }
  image.addFitsKeyword("RHLIMIT", image.rhLimit(), "Q_d = 0 for r_h > rhlimit");
  image.addFitsKeyword("SCATTER", image.scatteringMode(), "Light scattering on/off");
  image.addFitsKeyword("THERMAL", image.thermalMode(), "Thermal emisison on/off");
  image.addFitsKeyword("WAVELEN", image.wavelength(), "Simulated wavelength");
  image.addFitsKeyword("VLIMIT", image.vLimit(), "Vej = [0, v0*sqrt(beta/r_h)]");
  image.addFitsKeyword("NPOLE1", image.npole()[0], "North pole lambda");
  image.addFitsKeyword("NPOLE2", image.npole()[1], "North pole beta");
}

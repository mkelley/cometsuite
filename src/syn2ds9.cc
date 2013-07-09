/***************************************************************************

  Creates DS9 regions from rundynamics syndyne output.

  Copyright (C) 2005-2008,2010,2012 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <valarray>
#include <getopt.h>
#include "paramSet.h"
#include "particle.h"
#include "xyzstream.h"
#include "projection.h"
#include "longlat.h"
#include "StringConv.h"

#define SUBPROJECT "syn2ds9"

using namespace std;

int parseCommandLine(int, char**, vector<string>& ,string&, string&,
		     valarray<float>&, bool&, bool&, bool&);
long writeRegions(xyzstream&, ofstream&, paramSet, long, string,
		  valarray<float>, bool, bool);
void usage();

/******************************************************************************/
int main(int argc, char *argv[])
{
  string outfile, observerName = "Earth";
  valarray<float> offset(0.0, 2);
  vector<string> infiles;
  bool verbose = false, color = true, ecliptic = false;
  switch (parseCommandLine(argc, argv, infiles, outfile, observerName,
			   offset, verbose, color, ecliptic)) {
  case CL_NOERROR:
    break;
  case CL_HELP:
    return EXIT_SUCCESS;
  default:
    return EXIT_FAILURE;
  }

  xyzstream xyzfile;
  ofstream of;

  of.open(outfile.c_str());
  if (of.fail()) {
    cerr << "Error opening " << outfile << endl;
    return EXIT_FAILURE;
  }

  long n = 0, expectedTotal = 0, max;
  for (int i=0; i<infiles.size(); i++) {
    xyzfile.xyzopen(infiles[i].c_str(), xyzstream::READ);
    if (xyzfile.fail()) {
      cerr << "Error opening " << infiles[i] << endl;
      return EXIT_FAILURE;
    }

    // get the parameters
    paramSet parameters;
    string header = xyzfile.readHeader();
    stringstream str;
    if (verbose) cerr << header << endl;
    str << header;
    parameters.loadParameters(str);
    xyzfile.initData(parameters);

    long max;
    if (parameters.isSyndynes()) {
      max = static_cast<long>(parameters.beta().size() * parameters.steps());
    } else {
      cerr << "Files created with Make Comet do not contain syndynes!\n\n";
      xyzfile.close();
      continue;
    }

    expectedTotal += max;
    n += writeRegions(xyzfile, of, parameters, max, observerName, offset,
		      color, ecliptic);

    if (verbose) {
      cerr.precision(12);
      cerr << "Successfully read " << n << " out of " << max <<
	" particles\n";
    }

    xyzfile.close();
  }

  if (verbose) cerr << "DS9 regions written to " << outfile << "\n\n";
  of.close();

  if (n < expectedTotal) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

/******************************************************************************/
/** Checks the command line for proper input and loads a parameter
    file if possible. */
int parseCommandLine(int argc, char** argv, vector<string>& infiles,
		     string& outfile, string& observerName,
		     valarray<float>& offset, bool& verbose, bool& color,
		     bool& ecliptic) {
  bool files = false, help = false;

  // default output file is stdout
  outfile = "/dev/stdout";

  if (argc > 1) {
    int c;
    int digit_optind = 0;
    option longOptions[] = {
      {"color", 1, 0, 'c'},
      {"ecliptic", 0, 0, 0},
      {"help", 0, 0, 'h'},
      {"observer", 1, 0, 0},
      {"offset", 1, 0, 0},
      {"verbose", 1, 0, 'v'},
      {0, 0, 0, 0}
    };

    while (1) {
      int optionIndex = 0;
      int this_option_optind = optind ? optind : 1;

      c = getopt_long(argc, argv, "c:ho:v", longOptions,
		      &optionIndex);
      if (c == -1) break;

      switch(c) {
      case 0:
	if (longOptions[optionIndex].name == "ecliptic") {
	  ecliptic = true;
	  break;
	}
	if (longOptions[optionIndex].name == "observer") {
	  observerName = optarg;
	  break;
	}
	if (longOptions[optionIndex].name == "offset") {
	  offset = StringConv(optarg).toValarray<float>();
	  offset /= 3600;
	  break;
	}
	break; // stop processing long options

      case 'c':	color = StringConv(optarg).toBool(); break;
      case 'h': help = true; break;
      case 'o': outfile = optarg; break;
      case 'v': verbose = true;	break;
      case '?':
      default: break;
      }
    }

    while (optind < argc) {
      // add an input file
      infiles.push_back(string(argv[optind]));
      files = true;
      optind++;
    }

    if (!files || help) {
      usage();
      return CL_HELP;
    }

    return CL_NOERROR;
  } else {
    usage();
    return CL_HELP;
  }
}

/******************************************************************************/
/** Converts the syndynes to a regions file. */
long writeRegions(xyzstream& xyzfile, ofstream& of, paramSet parameters,
		  long max, string observerName, valarray<float> offset,
		  bool color, bool ecliptic) {
  // some region header info and initial parameters
  of << "# Region file format: DS9 version 4.0\n";
  parameters.writeParameters(of, "# rundynamics ", false);
  of << "global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n";
  if (ecliptic) {
    of << "ecliptic\n";
  } else {
    of << "fk5\n";
  }

  // set up the observer and add in the offset
  projection observer(observerName, parameters.obsDate(), ecliptic);
  observer.offset(offset);

  // set up the available colors
  vector<string> colorList;
  colorList = vector<string>(7);
  colorList[0] = "red";
  colorList[1] = "green";
  colorList[2] = "blue";
  colorList[3] = "cyan";
  colorList[4] = "magenta";
  colorList[5] = "yellow";
  colorList[6] = "white";

  // try and read all particles, connect the dots along the way
  long n;
  long s = 0;
  static int currentColor = 0;
  of.precision(9);
  particle p;
  longlat last;
  for (n=0; n<max; n++) {
    xyzfile.readParticle(p);
    longlat radec = observer.observe(p.fstate().r);

    if (s == 0) {
      last = radec;
      s++;
      //      of << "# beta " << p.beta() << "\n";
      //      of << "# color " << colorList[currentColor] << "\n";
    } else if (s >= parameters.steps()-1) {
      s = 0;
      currentColor++;
      currentColor = currentColor % 7;
    } else {
      of << "line(" << last.lambda << "," << last.beta << "," <<
	radec.lambda << "," << radec.beta << ") # line=0 0";

      if (color) of << "  color=" << colorList[currentColor];

      of << "\n";

      last = radec;
      s++;
    }
    
    if (p.error) break;
  }

  return n;
}

/******************************************************************************/
/** Prints the program usage and info. */
void usage() {
  cerr << "Usage: " << SUBPROJECT << " [OPTION]... xyzfile1.xyz [xyzfile2.xyz...]\n\
\n\
Mandatory arguments to long options are mandatory for short options too.\n\
    -c, --color=BOOL        make each syndyne a differnt color (default: on)\n\
    -h, --help              display this help\n\
    -o FILENAME             output to file FILENAME (default: stdout)\n\
    --observer=STRING       the name of the observer, Earth (default) or\n\
                            Spitzer\n\
    --offset=VAL            Coordinate offsets of the comet position in\n\
                            arcseconds (i.e., this does not include the\n\
                            cos(Dec) correction)\n\
    -v, --verbose           output more info than necessary\n\
\n\
BOOL may be one of {true, yes, on, 1, false, no, off, 0}.\n\
\n\
(c) 2005-2008,2010,2012 Michael S. Kelley\n";
  return;
}

/***************************************************************************

  The integratelist main program.

  Copyright (C) 2011,2012 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <getopt.h>
#include "rundynamics.h"
#include "logFile.h"
#include "paramSet.h"
#include "Integrator.h"
#include "ra15.h"
#include "state.h"
#include "xyzstream.h"

#define COMMENT string("# ")

#define SUBPROJECT "integratelist"

using namespace std;

int parseCommandLine(int, char**, vector<string>&, paramSet&, string&, double&);
void usage();
void timeAndStatus(long, clock_t&, string, bool, logFile&);

int main(int argc, char *argv[])
{
  paramSet temp, parameters, outParameters;
  particle p;
  state istate;
  clock_t start = clock();
  vector<string> xyzfileNames;
  string outfile;
  double jd = -1;
  double dt_obs = 0;  // seconds
  long max = 0;
  long n = 0;
  int status;

  switch (parseCommandLine(argc, argv, xyzfileNames, parameters, outfile, jd)) {
  case CL_NOERROR:
    break;
  case CL_HELP:
    return EXIT_SUCCESS;
  case CL_BADINPUT:
    cerr << "Error opening the first xyz file.\n";    
  default:
    return EXIT_FAILURE;
  }

  // Outputs
  xyzstream inxyz;
  logFile log(outfile);
  log.add(COMMENT + PACKAGE_STRING + "\n");

  // Scan the input files to determine how many grains we will read in
  for (int i=0; i<xyzfileNames.size(); i++) {
    inxyz.xyzopen(xyzfileNames[i], xyzstream::READ);
    if (inxyz.fail()) {
      cerr << "Error opening file: " << xyzfileNames[i] << "\n";
      return EXIT_FAILURE;
    }
    temp = inxyz.readParameters();
    if (temp.isSyndynes()) {
      max += static_cast<long>(temp.beta().size() * temp.steps());
    } else {
      max += temp.nParticles();
    }
    inxyz.close();
  }

  // For now, RADAU15 is the only integrator allowed
  Integrator *integrator = new ra15();
  for (int i=0; i<xyzfileNames.size(); i++) {
    log.add(COMMENT + "using input file: " + xyzfileNames[i] + "\n");
    log.add(COMMENT + "number radius beta age minStep nFC\n");

    cout << "Opening " << xyzfileNames[i] << " for input." << endl;
    inxyz.xyzopen(xyzfileNames[i], xyzstream::READ);
    if (inxyz.fail()) {
      cerr << "Error opening file: " << xyzfileNames[i] << "\n";
      return EXIT_FAILURE;
    }

    if (i == 0) {
      // The output file is only initialized once.  This means that
      // the output parameters set in the first input file will be
      // used for all output
      outParameters = parameters;
      if (jd > 0) {
	outParameters.obsDate(jd);
      } else {
	jd = outParameters.obsDate();
	outParameters.obsDate(jd);
      }
      outParameters.program("Integrate XYZ");
      outParameters.outFile(outfile);
      outParameters.nParticles(max);

      // Make sure the final state is saved
      outParameters.setSavedData("r_f");
      outParameters.setSavedData("v_f");
      outParameters.setSavedData("t_f");

      if (integrator->setup(outParameters)) {
	cerr << "Error setting up the Integrator.\n";
	return EXIT_FAILURE;
      }
    }

    // read in grains until the expected number is reached or an error
    // occurs; update the cumulative sums after every grain
    // get number of particles for syndynes or make comet

    if (parameters.isSyndynes()) {
      max = static_cast<long>(parameters.beta().size() * parameters.steps());
    } else {
      max = parameters.nParticles();
    }

    dt_obs = (jd - parameters.obsDate()) * 86400;
    if (DEBUG) cerr << "Adding " << dt_obs << " seconds to obsDate" << endl;
    for (n=0; n<max; n++) {
      inxyz.readParticle(p);
      if (DEBUG) cerr << n << " read in:\n" << p << endl;

      if (p.error) {
	cout << "Particle read error... exiting.\n\n";
	break;
      }

      if (p.beta() == -99) {
	// flag an error since we don't want to integrate the orbit
	status = true;
      } else {
	istate = p.istate();
	istate.t = -(p.age() + dt_obs);
	p.istate(istate);
	integrator->p(p);
	if (DEBUG) cerr << "initial state:\n" << p << endl;
	status = integrator->calculateOne();
	if (DEBUG) {
	  cerr << "Solution:\n";
	  cerr << "final r: " << integrator->p().fstate().r << endl;
	  cerr << "final v: " << integrator->p().fstate().v << endl;
	  cerr << "final t: " << integrator->p().fstate().t << "\n\n";
	}
      }

      /* Update the paraticle counter, even if a grain was not used
	 (e.g., it fell outside the integration box).  This way we can
	 easily compute the normalization factor when combining two
	 separate simulations.  Only update the log file if a grain
	 remained valid, otherwise the log file could grow very large,
	 very fast when many grains are being thrown away. */

      if ((status == 1) || DEBUG) {
	timeAndStatus(n, start, integrator->status, false, log);
      }
    }
    timeAndStatus(n, start, integrator->status, true, log);
    inxyz.close();
  }

  integrator->xyzfile.timeStamp();

  // clean up
  delete integrator;

  return EXIT_SUCCESS;
}

/** Checks the command line for proper input and loads a parameter
    file if possible. */
int parseCommandLine(int argc, char** argv, vector<string>& xyzfileNames,
		      paramSet& parameters, string& outfile, double& jd) {
  string keyword[1000], value[1000];
  int i = -1, example = 0;
  bool files = false, help = false;
  xyzstream inxyz;

  if (argc > 1) {
    int c;
    int digit_optind = 0;
    option longOptions[] = {
      {"help",            no_argument,       0, 'h'},
      {"jd",              required_argument, 0, 'j'},
      {"output",          required_argument, 0, 'o'},
      {"planets",         required_argument, 0, 0  },
      {"planetlookup",    required_argument, 0, 0  },
      {"xyzfile",         required_argument, 0, 'x'},
      {0, 0, 0, 0}
    };

    while (1) {
      int optionIndex = 0;
      int this_option_optind = optind ? optind : 1;

      c = getopt_long(argc, argv, "hj:o:t:x:", longOptions,
			  &optionIndex);
      if (c == -1) break;

      // on the off chance there are more than 1000 parameters on the
      // command line
      if (i > 999) i = 998;

      switch(c) {
      case 0:
	if (longOptions[optionIndex].name == "planets") {
	  keyword[++i] = "PLANETS";
	  value[i] = optarg;
	  break;
	}

	if (longOptions[optionIndex].name == "planetlookup") {
	  keyword[++i] = "PLANETLOOKUP";
	  value[i] = optarg;
	  break;
	}

	break;

      case 'h':
	help = true;
	break;

      case 'j':
	jd = atof(optarg);
	break;

      case 'o':
      case 'x':
	outfile = optarg;
	break;

      case '?':
      default:
	cerr << "I don't understand the command line parameter: " << c << endl;
	break;
      }
    }

    while (optind < argc) {
      // add an input file
      xyzfileNames.push_back(string(argv[optind]));
      files = true;
      optind++;
    }

    if (!files || help) {
      usage();
      return CL_HELP;
    }

    // Read in the first parameter set
    inxyz.xyzopen(xyzfileNames[0], xyzstream::READ);
    if (inxyz.fail()) return CL_BADINPUT;
    parameters = inxyz.readParameters();
    inxyz.close();

    // Set any parameters specified on the command line
    for (int j=0; j<=i; j++) {
      // nothing should start with a equals sign, this likely means
      // a short parameter was entered '-s=100'
      if (value[j][0] == '=') value[j][0] = ' ';
      if (DEBUG) cerr << keyword[j] << " " << value[j] << "\n";
      parameters.setParameter(keyword[j], value[j]);
    }

    return CL_NOERROR;
  } else {
    usage();
    return CL_HELP;
  }
}

/** Prints the time spent calculating particle positions and updates
    the log file with the last particle\'s status. */
void timeAndStatus(long n, clock_t& start, string status, bool last,
		   logFile& log) {
  // write out the status string
  log.add(n); log.add(status);

  /* Give the user immediate feedback for the first 25 particles,
     then, every 25 particles. */
  if (n <= 25) log.flush();
  if ((n % 25) == 0) log.flush();

  // print time on every 1000th particle or if we are on the last
  // particle
  if (((n % 1000) == 0) || last) {
      clock_t now = clock();
      cout << "Particle " << n << " completed (";
      cout << static_cast<double>(now - start) /
	static_cast<double>(CLOCKS_PER_SEC) / 1000;
      cout << " cpu seconds/particle)\n";
      start = clock();
  }
}

/** Prints the program usage and info. */
void usage() {
  cerr << PACKAGE_STRING << "\n";
  cerr << "Usage: " << SUBPROJECT << " [OPTION]... xyzfile1.xyz [xyzfile2.xyz...]\n\
Mandatory arguments to long options are mandatory for short options too.\n\
    -h, --help              display this help\n\
    -j, --jd=DATE           advance the xyzfile to DATE (overrides file's JD)\n\
    -o                      same as \'--xyzfile\'\n\
    -x, --xyzfile=FILE      set the output file to FILE\n\
\n\
BOOL may be one of {true, yes, on, 1, false, no, off, 0}.\n\
Parameter values may be enclosed in quotes, e.g. --beta=\"1e-3 2e-3 4e-3\" or\n\
  -b \"0.1 0.01 0.001\".\n\
Command line parameters override the parameter file.\n\
\n\
(c) 2011 Michael S. Kelley\n";
}

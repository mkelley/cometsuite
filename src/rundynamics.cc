/***************************************************************************

  The rundynamics program.

  Copyright (C) 2005,2006,2007,2008,2009 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <getopt.h>
#include "rundynamics.h"
#include "logFile.h"
#include "paramSet.h"
#include "Integrator.h"
#include "ra15.h"

#define COMMENT string("# ")

#define SUBPROJECT "rundynamics"

using namespace std;

int parseCommandLine(int, char**, paramSet&, string&);
void usage();
void timeAndStatus(long, clock_t&, string, bool, logFile&);

/** \todo Parallelize?
    \todo Incorporate new integrators.
 */
int main(int argc, char *argv[])
{
  paramSet parameters;
  string infile;
  int r;

  switch (parseCommandLine(argc, argv, parameters, infile)) {
  case NOERROR:
    break;
  case HELP:
    return EXIT_SUCCESS;
  case NOFILE:
    cerr << "Error parsing the command line.\n";
  case BADINPUT:
    cerr << "Error parsing the parameter file.\n";    
  default:
    return EXIT_FAILURE;
  }

  logFile log(parameters.outFile());
  log.add(COMMENT + PACKAGE_STRING + "\n");
  log.add(COMMENT + "using input file: " + infile + "\n");
  log.add(COMMENT + "number radius beta age minStep nFC\n");

  cout << "\n";
  parameters.writeParameters(cout);
  cout << endl;

  // For now, RADAU15 is the only integrator allowed
  Integrator *integrator = new ra15();

  if (integrator->setup(parameters)) {
    cerr << "Error setting up the Integrator.\n";
    return EXIT_FAILURE;
  }

  clock_t start = clock();
  long n = 0;
  int status;
  while ((status = integrator->calculateOne())) {
    /* Update the paraticle counter, even if a grain was not used
       (e.g., it fell outside the integration box).  This way we can
       easily compute the normalization factor when combining two
       separate simulations.  Only update the log file if a grain
       remained valid, otherwise the log file could grow very large,
       very fast when many grains are being thrown away. */
    n++;
    if (DEBUG) cerr << n << endl;

    if ((status == 1) || DEBUG) {
      timeAndStatus(n, start, integrator->status, false, log);
    }
  }
  timeAndStatus(++n, start, integrator->status, true, log);

  integrator->xyzfile.timeStamp();

  // clean up
  delete integrator;

  return EXIT_SUCCESS;
}

/** Checks the command line for proper input and loads a parameter
    file if possible. */
int parseCommandLine(int argc, char** argv, paramSet& parameters,
		      string& infile) {
  string keyword[1000], value[1000];
  int i = -1, example = 0;
  bool file = false, help = false;

  if (argc > 1) {
    int c;
    int digit_optind = 0;
    option longOptions[] = {
      {"beta",            required_argument, 0, 'b'},
      {"box",             required_argument, 0, 0  },
      {"closeapproaches", required_argument, 0, 0  },
      {"comet",           required_argument, 0, 'c'},
      {"example",         no_argument,       0, 0  },
      {"help",            no_argument,       0, 'h'},
      {"jd",              required_argument, 0, 'j'},
      {"kernel",          required_argument, 0, 'k'},
      {"ltt",             required_argument, 0, 0  },
      {"ndays",           required_argument, 0, 0  },
      {"nparticles",      required_argument, 0, 0  },
      {"orbit",           required_argument, 0, 0  },
      {"pfunc",           required_argument, 0, 0  },
      {"planets",         required_argument, 0, 0  },
      {"planetlookup",    required_argument, 0, 0  },
      {"program",         required_argument, 0, 0  },
      {"steps",           required_argument, 0, 's'},
      {"tol",             required_argument, 0, 't'},
      {"xyzfile",         required_argument, 0, 'x'},
      {0, 0, 0, 0}
    };

    while (1) {
      int optionIndex = 0;
      int this_option_optind = optind ? optind : 1;

      c = getopt_long(argc, argv, "b:c:ehj:k:l:o:s:t:x:", longOptions,
			  &optionIndex);
      if (c == -1) break;

      // on the off chance there are more than 1000 parameters on the
      // command line
      if (i > 999) i = 998;

      switch(c) {
      case 0:
	if (longOptions[optionIndex].name == "closeapproaches") {
	  keyword[++i] = "CLOSEAPPROACHES";
	  value[i] = optarg;
	  break;
	}

	if (longOptions[optionIndex].name == "box") {
	  keyword[++i] = "BOX";
	  value[i] = optarg;
	  break;
	}

	if (longOptions[optionIndex].name == "example") {
	  example += 2;
	  break;
	}

	if (longOptions[optionIndex].name == "ltt") {
	  keyword[++i] = "LTT";
	  value[i] = optarg;
	  break;
	}

	if (longOptions[optionIndex].name == "ndays") {
	  keyword[++i] = "NDAYS";
	  value[i] = optarg;
	  break;
	}

	if (longOptions[optionIndex].name == "nparticles") {
	  keyword[++i] = "NPARTICLES";
	  value[i] = optarg;
	  break;
	}

	if (longOptions[optionIndex].name == "orbit") {
	  keyword[++i] = "ORBIT";
	  value[i] = optarg;
	  break;
	}

	if (longOptions[optionIndex].name == "pfunc") {
	  keyword[++i] = "PFUNC";
	  value[i] = optarg;
	  break;
	}

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

	if (longOptions[optionIndex].name == "program") {
	  keyword[++i] = "PROGRAM";
	  value[i] = optarg;
	  break;
	}

	break;

      case 'b':
	keyword[++i] = "BETA";
	value[i] = optarg;
	break;

      case 'c':
	keyword[++i] = "COMET";
	value[i] = optarg;
	break;

      case 'e':
	example++;
	break;

      case 'h':
	help = true;
	break;

      case 'j':
	keyword[++i] = "JD";
	value[i] = optarg;
	break;

      case 'k':
	keyword[++i] = "KERNEL";
	value[i] = optarg;
	break;

      case 's':
	keyword[++i] = "STEPS";
	value[i] = optarg;
	break;

      case 't':
	keyword[++i] = "TOL";
	value[i] = optarg;
	break;

      case 'o':
      case 'x':
	keyword[++i] = "XYZFILE";
	value[i] = optarg;
	break;

      case '?':
      default:
	cerr << "I don't understand the command line parameter: " << c << endl;
	break;
      }
    }

    if (optind < argc) {
      file = true;
      infile = string(argv[optind]);
    }

    if ((!file || help) && !example) usage();

    if (example) {
      if (example > 1)
	parameters.writeParameters(cout, true);
      else
	parameters.writeParameters(cout, false);
      return 1;
    }

    // Stop here if there is no file to load
    if (!file) return NOFILE;

    // Load the parameter file, return 3 on an error
    if (!parameters.loadParameters(infile)) return BADINPUT;

    // Set any parameters specified on the command line
    for (int j=0; j<=i; j++) {
      // nothing should start with a equals sign, this likely means
      // a short parameter was entered '-s=100'
      if (value[j][0] == '=') value[j][0] = ' ';
      if (DEBUG) cerr << keyword[j] << " " << value[j] << "\n";
      parameters.setParameter(keyword[j], value[j]);
    }

    return NOERROR;  // no errors
  } else {
    usage();
    return HELP;
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
  cerr << "Usage: " << SUBPROJECT << " [OPTION]... parameterFilename.par\n\
\n\
Mandatory arguments to long options are mandatory for short options too.\n\
    -b, --beta=VAL          set the beta values to VAL\n\
    --box==VAL              set the integration box to VAL km or set\n\
                            VAL < 0 to disable (default: disabled)\n\
    --closeapproaches=BOOL  switch close approach handling on/off\n\
    -c, --comet=NAME        set the comet to NAME\n\
    -e                      print an example parameter file\n\
    --example               print a commented example parameter file\n\
    -h, --help              display this help\n\
    -j, --jd=DATE           set the Julian date to DATE\n\
    -k, --kernel=FILE       set the kernel to FILE\n\
    --ltt=BOOL              switch light travel time on/off\n\
    --ndays=VAL             set ndays to VAL\n\
    --nparticles=VAL        set nparticles to VAL\n\
    -o                      same as \'--xyzfile\'\n\
    --orbit=VAL             set orbit to VAL\n\
    --pfunc=FUNC            set the particle function (and parameters) to FUNC\n\
    --planets=VAL           set the planets to VAL\n\
    --planetlookup=BOOL     switch the planet look up table on/off\n\
    --program=NAME          set the program to NAME\n\
    -s, --steps=VAL         set steps to VAL\n\
    -t, --tol=VAL           set the tolerance to VAL\n\
    -x, --xyzfile=FILE      set the output file to FILE\n\
\n\
BOOL may be one of {true, yes, on, 1, false, no, off, 0}.\n\
Parameter values may be enclosed in quotes, e.g. --beta=\"1e-3 2e-3 4e-3\" or\n\
  -b \"0.1 0.01 0.001\".\n\
Command line parameters override the parameter file.\n\
\n\
(c) 2005-2010 Michael S. Kelley\n";
}

/***************************************************************************

  Implements rundynamics parameter files.

  Copyright (C) 2005-2010,2012 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <ctime>
#include <unistd.h>
#include "paramSet.h"
#include "rundynamics.h"
#include "getxyz.h"
#include "StringConv.h"

using namespace std;

/** The parameter set class constructor using default parameters. */
paramSet::paramSet() {
  _program = "syndynes";
  _comet = "encke";
  _spkKernel = "encke.bsp";
  //obsDate(2004, 1, 1.0); // note: specifying calendar date doesn't work :)
  _obsDate = 2450643.54170;
  _outFile = "output.xyz";
  _labelFormat = "";
  _pFunc = "";
  _tolerance = 0.01;
  _planets = 511;
  _planetLookUp = false;
  _closeApproaches = true;
  _box = -1;
  _lightTravelTime = false;
  _data = "radius graindensity beta age origin r_i v_ej r_f";  // default for v0.7.3+
  beta(StringConv("1e-3 2e-3 4e-3 6e-3 8e-3 1e-2 1e-1").toVector<double>());
  _nDays = 200;
  _steps = 31;
  _orbit = 1;
  _nParticles = 100000;

  // default for v0.5.0 - v0.7.2
  // _data = "beta age origin r_i v_ej r_f";

  // even older default:
  // _data = "beta radius composition density fractaldim age origin v_ej r_i v_i t_i r_f v_f t_f";

  /* type(name), where type is d for double, i for integer and may
     include [n] if the variable is an array (e.g., d[3]).  Make sure
     each output name is unique such that it doesn't appear in any
     other output name (e.g., cook and cookie) and that there is no
     whitespace in the name. */
  _dataList.push_back(" d(beta)");         _unitList.push_back(" none");
  _dataList.push_back(" d(radius)");       _unitList.push_back(" micron");
  _dataList.push_back(" d(graindensity)"); _unitList.push_back(" g/cm^3");
  _dataList.push_back(" d(age)");          _unitList.push_back(" s");
  _dataList.push_back(" d[2](origin)");    _unitList.push_back(" deg");
  _dataList.push_back(" d[3](v_ej)");      _unitList.push_back(" km/s");
  _dataList.push_back(" d[3](r_i)");       _unitList.push_back(" km");
  _dataList.push_back(" d[3](v_i)");       _unitList.push_back(" km/s");
  _dataList.push_back(" d(t_i)");          _unitList.push_back(" s");
  _dataList.push_back(" d[3](r_f)");       _unitList.push_back(" km");
  _dataList.push_back(" d[3](v_f)");       _unitList.push_back(" km/s");
  _dataList.push_back(" d(t_f)");          _unitList.push_back(" s");
  _dataList.push_back(" c[16](label)");    _unitList.push_back(" none");

  _nData = _dataList.size();
}

/** The parameter set class constructor reading in a parameter file. */
paramSet::paramSet(const string filename) {
  if (!loadParameters(filename))
    cerr << "Error loading parameter file: " << filename;
}

/** The parameter set class constructor copying from another parameter set. */
paramSet::paramSet(const paramSet& p) { *this = p; }

/** The parameter set class constructor reading parameters from a stream. */
paramSet::paramSet(istream& inf) {
  if (!loadParameters(inf))
    cerr << "Error loading parameters from stream\n";
}

/** Returns true if the parameter is a valid data name, false
    otherwise. */
bool paramSet::isValidData(const string name) {
  string _name = getDataName(name);
  for (int i=0;i<_nData;i++) {
    int c = (_dataList[i]).find(_name);
    if (c != string::npos) return true;
  }
  return false;
}

/** Takes a string of the format "type(name)", "type[n](name)", or
    "name", with any amount of whitespace, and returns "name". */
string paramSet::getDataName(const string name) {
  int c;
  string _name = name;

  // remove extra spaces
  while ((c = _name.find(" ")) != string::npos) {
    _name.erase(c);
  }

  if ((c = _name.find("(")) != string::npos) {
    // remove everything up to "("
    _name.erase(0, c+1);
  }

  if ((c = _name.find(")")) != string::npos) {
    // remove ")" and everything after
    _name.erase(c, _name.length());
  }
  return _name;
}

/** Disables all parameters for saving. */
void paramSet::unsetSavedData() {
  _data = "";
  return;
}

/** Returns true if the parameter named by the passed string is
    enabled for saving, false otherwise. */
bool paramSet::isSavedData(const string name) {
  int c;
  string _name = getDataName(name);

  // remove extra spaces
  while ((c = _name.find(" ")) != string::npos) {
    _name.erase(c);
  }

  // make sure it is a legal name
  if (!isValidData(_name)) {
    cerr << name << " is not a valid data name\n";
    return false;
  }

  c = _data.find(_name);
  if (c != string::npos) {
    return true;
  } else {
    return false;
  }
}

/** Returns the description of the passed datum name, or blank if the
    name is invalid. */
string paramSet::getDataDescription(const string name) {
  string _name = getDataName(name);
  for (int i=0;i<_nData;i++) {
    int c = (_dataList[i]).find(_name);
    if (c != string::npos) return _dataList[i];
  }
  return "";
}

/** Enable a parameter for saving. */
void paramSet::setSavedData(const string name) {
  setSavedData(name, true);
}

/** Disable a parameter for saving. */
void paramSet::unsetSavedData(const string name) {
  setSavedData(name, false);
}

/** Enable or disable a parameter for saving. */
void paramSet::setSavedData(const string name, const bool enable) {
  int c;
  string _name = getDataName(name);

  if (enable) {
    if (isSavedData(_name)) {
      return;
    } else {
      _data.append(" ");
      _data.append(_name);
      return;
    }
  } else {
    if (isSavedData(_name)) {
      while ((c = _data.find(_name)) != string::npos) {
	_data.erase(c, _name.length());
      }
    }
  }
  return;
}

/** Enables saving of a list of data values (e.g., those specified by
    the SAVE parameter or the internal DATA parameter). */
void paramSet::setSavedDataList(const string dataList) {
  stringstream data;
  data << dataList;
  if (DEBUG) cerr << "Setting saved data list to " << dataList << " (paramSet)\n";
  while (!data.eof()) {
    string datum;
    data >> datum;
    datum = getDataName(datum);
    if (isValidData(datum)) setSavedData(datum);
  }
  return;
}

/** Return a string listing all the saved data names. */
string paramSet::savedData() {
  return _data;
}

/** Return a string describing the saved data, i.e.,
    type(name) or type[n](name). */
string paramSet::savedDataDescription() {
  string description;
  stringstream data;
  data << _data;
  while (!data.eof()) {
    string datum;
    data >> datum;
    if (isValidData(datum))
      description.append(getDataDescription(datum));
  }
  return description;
}

/** Return a string describing the units of the saved data. */
string paramSet::savedDataUnits() {
  string units;
  stringstream data;
  data << savedDataDescription();
  while (!data.eof()) {
    string datum;
    data >> datum;
    for (int i=0;i<_nData;i++) {
      int c = (_dataList[i]).find(datum);
      if (c != string::npos)
	units.append(_unitList[i]);
    }
  }
  return units;
}

/** Loads a parameter file. Returns true on success. */
bool paramSet::loadParameters(const string filename) {
  ifstream inf(filename.c_str(), ios::in);

  if (!inf.good()) {
    cerr << "Error opening file: " << filename << "\n";
    return false;
  }

  bool ret = loadParameters(inf);
  inf.close();
  return ret;
}

/** Reads parameters from a stream. Returns true on success.

    \todo Refine so that spaces aren't required after the colon.
*/
bool paramSet::loadParameters(istream& inf) {
  while(!inf.eof()) {
    char buf[255];
    string parameter, value;
    stringstream str;
  
    inf.getline(buf, 255);
    str << buf;

    // ignore comments
    str >> parameter;
    if (parameter[0] != '#') {
      // make sure this is a valid parameter, i.e. ends with a colon
      if (parameter[parameter.length()-1] == ':') {
	// remove the :
	parameter.resize(parameter.length()-1);

	// get the parameter value(s), if it exists
	str.getline(buf, 255);
	if (!str.fail()) {
	  value = buf;

	  // remove leading whitespace, esp. important for string values
	  int n = value.find_first_not_of(" ");
	  if (n > 0) value.erase(0, n);

	  // if a value exists, then set the parameter
	  if (value.length() > 0) setParameter(parameter, value);
	}
      }
    }
  }

  char kfile[255];
  strcpy(kfile, spkKernel().c_str());
  findkernel(comet().c_str(), kfile);
  spkKernel(kfile);

  return true;
}

/** Saves a parameter file.
 *  Returns true on success. */
bool paramSet::saveParameters(const string filename) {
  ofstream saveFile(filename.c_str());
  if (!saveFile.is_open()) return false;

  time_t t;
  time(&t);
  saveFile << "# " << filename << "  " << asctime(localtime(&t));
  writeParameters(saveFile);
  saveFile.close();

  cout << "Parameters saved to " << filename << endl;
  return true;
}

/** Writes parameters to a stream without a prefix. */
void paramSet::writeParameters(ostream& os) {
  writeParameters(os, "", false);
}

/** Writes parameters to a stream without a prefix. */
void paramSet::writeParameters(ostream& os, const bool extended) {
  writeParameters(os, "", extended);
}

/** Writes parameters to a stream with a prefix. */
void paramSet::writeParameters(ostream& os,
			       const string pre,
			       const bool extended) {
  int p = os.precision();
  os.precision(16);

  os << pre << "# " << PACKAGE_STRING << "\n";
  if (extended) {
    os << pre << "\n";
    os << pre << "################################################################################\n";
    os << pre << "# The name of the program you wish to run.\n";
  }
  os << pre << "# Valid program names: syndynes, make comet, integratexyz\n";
  os << pre << "PROGRAM: " << program() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "################################################################################\n";
  }
  os << pre << "# Parameters common to all programs.\n";
  if (extended) {
    os << pre << "#\n";
    os << pre << "# Numbers may be entered as an integer, as a floating point value, or in\n";
    os << pre << "# exponential notation.\n";
    os << pre << "#\n";
    os << pre << "# Valid boolean values: true/on/yes/1 and false/off/no/0 are equivalent\n";
    os << pre << "\n";
    os << pre << "# The name of the comet in the SPICE kernel.  Capitalization is\n";
    os << pre << "# unimportant.  Be sure to check your kernel for the appropriate comet\n";
    os << pre << "# name.  For example, C/2001 Q4 (NEAT) is referred by its NAIF ID\n";
    os << pre << "# 1000351.\n";
  }      
  os << pre << "COMET: " << comet() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# The SPICE kernel file name.  By default " << PACKAGE_NAME << " removes all\n";
    os << pre << "# non-alphanumeric characters from the comet name and appends\n";
    os << pre << "# \".bsp\".  For example: COMET: shoemaker-levy 9 --> KERNEL:\n";
    os << pre << "# shoemakerlevy9.bsp\n";
  }
  os << pre << "KERNEL: " << spkKernel() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# The Julian date of the observation.\n";
  }
  os << pre << "JD: " << obsDate() << endl;

  if (extended) {
    os << pre << "# The name of the output file.  The convention is to use a suffix of\n";
    os << pre << "# .xyz\n";
  }
  os << pre << "XYZFILE: " << outFile() << endl;

  if (extended) {
    os << pre << "# The label format string.  Set this if you want each particle in your\n";
    os << pre << "# xyzfile to contain a label string.  Use \"%d\" if you want the particle's\n";
    os << pre << "# number to be saved, making the each partcile's string unique (if you want\n";
    os << pre << "# labels, then you probably want this).  Labels are only saved if \"label\" is\n";
    os << pre << "# added to the SAVE parameter list.  For example: LABEL: fast%d --> creates\n";
    os << pre << "# labels named fast1, fast2, fast3, etc..  Labels may be no longer than 16\n";
    os << pre << "# characters, including the particle number, if used.\n";
  }
  os << pre << "LABEL: " << labelFormat() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# The particle function generator.  Primarily used by Make Comet but\n";
    os << pre << "# may also be used by Syndynes to generate ejection velocities.  The\n";
    os << pre << "# PFUNC string has the format:\n";
    os << pre << "#\n";
    os << pre << "#   PFUNC: pfunc_name parameters; pfunc_name parameter; ...\n";
    os << pre << "#\n";
    os << pre << "# i.e., a list of particle functions and parameters separated by a\n";
    os << pre << "# semi-colon.  Capitalization is not important.  The particle\n";
    os << pre << "# functions are read in left to right order, and the last function\n";
    os << pre << "# will override any previously set variables.  The particle functions\n";
    os << pre << "# may be any of the following modules or templates.  Those marked with\n";
    os << pre << "# ** are not yet implemented.  Quotes should be omitted.  Logarithms\n";
    os << pre << "# are base 10.  Use log10(radius) = -999 to force beta = 0.  Any\n";
    os << pre << "# radius less than or equal to zero will force beta = 0.  Units \n";
    os << pre << "# are AU, km/s, days, and degrees.\n";
    os << pre << "#\n";
    os << pre << "# See setup() in pfunctions.cc or \n";
    os << pre << "# cometsuite/doc/rundynamics/classpfunctions.html for more information.\n";
    os << pre << "#\n";
    os << pre << "# Modules -\n";
    os << pre << "#    RADIUS a\n";
    os << pre << "#    RADIUS min max [steps]\n";
    os << pre << "#    LOGRADIUS log10(a)\n";
    os << pre << "#    LOGRADIUS log10(min) log10(max) [steps]\n";
    os << pre << "#    AGE t\n";
    os << pre << "#    AGE min max [steps]\n";
    os << pre << "#    COMPOSITION name\n";
    os << pre << "#      Available materials (and their shorthand names):\n";
    os << pre << "#        geometric (g) - the default composition\n";
    os << pre << "#        amorphous carbon (ac)\n";
    os << pre << "#        amorphous olivine 50 (ol50)\n";
    os << pre << "#      For syndynes use geometric.  See pfunctions::composition()\n";
    os << pre << "#      more details.\n";
    os << pre << "#    BULKDENSITY rho\n";
    os << pre << "#      Bulk material density for geometric grains\n";
    os << pre << "#    FRACTALDIM D\n";
    os << pre << "#      Fractal dimension for fractally porous grains, must be one of:\n";
    os << pre << "#      3.0 (default), 2.857, 2.727, 2.609, 2.5\n";
    os << pre << "#    RADIUSLAW u1\n";
    os << pre << "#    RHLAW u2\n";
    os << pre << "#    VELOCITY (\"iso\"|\"cos\"|\"temp\") v0\n";
    os << pre << "#    VELOCITY (\"iso\"|\"cos\"|\"temp\") \"limit\" v0\n";
    os << pre << "#    VELOCITY (\"iso\"|\"cos\"|\"temp\") \"normal\" v0 sigma [absolute]\n";
    os << pre << "#    VELOCITY (\"iso\"|\"cos\"|\"temp\") \"range\" min max [steps]\n";
    os << pre << "#    RHLIMIT max\n";
    os << pre << "#    SUNCONE min max\n";
    os << pre << "#    POLE nPoleLambda nPoleBeta\n";
    os << pre << "#    JET longitude latitude angle period  (pole must already be set)\n";
    os << pre << "#\n";
    os << pre << "# Templates -\n";
    os << pre << "#   BASIC_EJECTION_VELOCITY v0\n";
    os << pre << "#   ** ONE_JET v0 [minAge] maxAge log(radiusMin) log(radiusMax) nPoleLambda nPoleBeta jetLat alpha period\n";
    os << pre << "#   ** ONE_SIMPLE_JET v0 [minAge] maxAge log(radiusMin) log(radiusMax) nPoleLambda nPoleBeta jetLat alpha period\n";
    os << pre << "#   SIMPLE_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n";
    os << pre << "#   ISO_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n";
    os << pre << "#   QV_COS_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n";
    os << pre << "#   QV_TEMP_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n";
    os << pre << "#   RTV log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0\n";
    os << pre << "#   RTV_GRID steps log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0\n";
    os << pre << "#   ** RTVL log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0 axisLam axisBet minLat maxLat\n";
    os << pre << "#\n";
    os << pre << "# For example, the following lines are equivalent:\n";
    os << pre << "#   PFUNC: simple_coma 1.35 2400 -1 4\n";
    os << pre << "#   PFUNC: velocity iso 1.35; age 0 2400; logradius 1 4; q_d iso; suncone 0 90\n";
  }
  os << pre << "PFUNC: " << pFunc() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# The fractional tolerance on the 15th order term in the integration\n";
    os << pre << "# routine.  In practice 1e-2 works well.  This is about the limit of\n";
    os << pre << "# double precision.\n";
  }
  os << pre << "TOL: " << tolerance() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# Parameter to specify which planets to include.  Pick planets from\n";
    os << pre << "# the following table and sum together:\n";
    os << pre << "#\n";
    os << pre << "#   Mercury = 2^0 = 1\n";
    os << pre << "#   Venus   = 2^1 = 2\n";
    os << pre << "#   Earth   = 2^2 = 4\n";
    os << pre << "#   Mars    = 2^3 = 8\n";
    os << pre << "#   Jupiter = 2^4 = 16\n";
    os << pre << "#   Saturn  = 2^5 = 32\n";
    os << pre << "#   Uranus  = 2^6 = 64\n";
    os << pre << "#   Neptune = 2^7 = 128\n";
    os << pre << "#   Pluto   = 2^8 = 256\n";
    os << pre << "#\n";
    os << pre << "# In general, most planets should be used (255 or 511).  In tests\n";
    os << pre << "# where the planet look up table is disabled, Pluto may be dropped\n";
    os << pre << "# (255) for a minor increase in performance (5 - 10%) with a\n";
    os << pre << "# negligible change in accuracy.\n";
  }
  os << pre << "PLANETS: " << planets() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# Switch to enable the planet look up table.  In general, leave this\n";
    os << pre << "# off for Syndynes and small Make Comet simulations for improved\n";
    os << pre << "# accuracy\n";
  }
  os << pre << "PLANETLOOKUP: " << planetLookUp() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# Switch to bypass the planet look up table when a comet approaches a\n";
    os << pre << "# planet.  In general, leave this on.  This is ignored when\n";
    os << pre << "# PLANETLOOKUP is off.\n";
  }
  os << pre << "CLOSEAPPROACHES: " << closeApproaches() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# Restrict the particle-comet distance to be less than or equal to\n";
    os << pre << "# this value, measured in km.  Set to less than zero to disable\n";
    os << pre << "# (default).  Particles with particle-comet distances greater than BOX\n";
    os << pre << "# are deleted from the simulation.  Not intended for use with\n";
    os << pre << "# Syndynes.\n";
  }
  os << pre << "BOX: " << box() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# Enable LTT to take into account the light travel time between an\n";
    os << pre << "# Earth-bound observer and the comet.  Usually this is left off.\n";
  }
  os << pre << "LTT: " << lightTravelTime() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# Specify which variables to save for each particle.  Full available\n";
    os << pre << "# list:\n";
    os << pre << "#   radius graindensity beta age origin v_ej r_i v_i t_i r_f v_f t_f label\n";
    os << pre << "# Unless you know what you are doing, it is best to only add to the\n";
    os << pre << "# default parameter list.\n";
  }
  os << pre << "SAVE: " << savedData() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "################################################################################\n";
  }
  os << pre << "# Syndyne specific section.\n";
  if (extended) {
    os << pre << "\n";
    os << pre << "# A list of beta parameter values.\n";
  }
  os << pre << "BETA:";
  for (int i = 0; i < _beta.size(); i++) {
    os << " " << _beta[i];
  }
  os << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# The maxium particle age.\n";
  }
  os << pre << "NDAYS: " << nDays() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# The number of evenly spaced time steps used to determine the age of\n";
    os << pre << "# each particle.  The age of particle i is: NDAYS / (STEPS - 1) * i,\n";
    os << pre << "# where i ranges from 0 to STEPS-1.\n";
  }
  os << pre << "STEPS: " << steps() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "# Save the comet's position for JD +/- ORBIT days along with the\n";
    os << pre << "# syndynes.  STEPS steps of the the orbit will be saved with the beta\n";
    os << pre << "# value -99.\n";
  }
  os << pre << "ORBIT: " << orbit() << endl;

  if (extended) {
    os << pre << "\n";
    os << pre << "################################################################################\n";
  }
  os << pre << "# Make comet specific section.\n";
  if (extended) {
    os << pre << "\n";
    os << pre << "# The number of particles to generate.\n";
  }
  os << pre << "NPARTICLES: " << _nParticles << endl;
  os.precision(p);
}

/** Returns true if the parameter set is for Syndynes. */
bool paramSet::isSyndynes() {
  string program = _program;
  transform(program.begin(), program.end(), program.begin(), (int(*)(int))toupper);
  if (program.find("SYNDYNES") != string::npos) return true;
  return false;
}

/** Returns true if the parameter set is for Make Comet. */
bool paramSet::isMakeComet() {
  string program = _program;
  transform(program.begin(), program.end(), program.begin(), (int(*)(int))toupper);
  if (program.find("MAKE COMET") != string::npos) return true;
  return false;
}

/** Returns true if the parameter set is for Integrate XYZ. */
bool paramSet::isIntegrateXYZ() {
  string program = _program;
  transform(program.begin(), program.end(), program.begin(), (int(*)(int))toupper);
  if (program.find("INTEGRATE XYZ") != string::npos) return true;
  return false;
}

/** Return the dynamics program name (syndynes, make comet, or
    integrate xyz). */
string paramSet::program() { return _program; }

/** Set the dynamics program name. */
void paramSet::program(const string progName) { _program = progName; }

/** Return the comet name. */
string paramSet::comet() { return _comet; }

/** Set the comet name. */
void paramSet::comet(const string name) { _comet = name; }

/** Return the Spice kernel file name. */
string paramSet::spkKernel() { return _spkKernel; }

/** Change the Spice kernel file name. */
void paramSet::spkKernel(const string kernel) { _spkKernel = kernel; }

/** Return the observation date.
 *  The observation date is on the Julian calendar. */
double paramSet::obsDate() { return _obsDate; }

/** Change the observation date to the given Julian date. */
void paramSet::obsDate(const double jd) { _obsDate = jd; }

/** Change the observation date to the given UT date.
 *  Day is given as a decimal.  This function does not yet work. */
void paramSet::obsDate(const int year, const int month, const double day) {
	_obsDate = 2450000.5;
}

/** Change the observation date to that in the passed string.
 *  The string may be a Julian date or of the format "year month day",
    "year, month, day", or "year-month-day". */
void paramSet::obsDate(const string) { _obsDate = 24500000.5; }

/** Return the output file name. */
string paramSet::outFile() { return _outFile; }

/** Change the output file name. */
void paramSet::outFile(const string filename) { _outFile = filename; }

/** Return the particle function name. */
string paramSet::pFunc() { return _pFunc; }

/** Change the particle function.  Force the entire string to
    lower-case. */
void paramSet::pFunc(const string fname) {
  _pFunc = fname;
  transform(_pFunc.begin(), _pFunc.end(), _pFunc.begin(), (int(*)(int))tolower);
}

/** Return the tolerance of the integrator. */
double paramSet::tolerance() { return _tolerance; }

/** Change the tolerance of the integrator. */
void paramSet::tolerance(const double t) { _tolerance = t; }

/** Return the planets being used. */
int paramSet::planets() { return _planets; }

/** Change the planets being used. */
void paramSet::planets(const int p) { _planets = p; }

/** Return the planet look up flag. */
bool paramSet::planetLookUp() { return _planetLookUp; }

/** Change the planet look up flag. */
void paramSet::planetLookUp(const bool pl) { _planetLookUp = pl; }

/** Return the close approaches flag. */
bool paramSet::closeApproaches() { return _closeApproaches; }

/** Change the close approaches flag. */
void paramSet::closeApproaches(const bool ca) { _closeApproaches = ca; }

/** Return the integration box size. */
double paramSet::box() { return _box; }

/** Change the integration box size. */
void paramSet::box(const double b) { _box = b; }

/** Return the light travel time flag. */
bool paramSet::lightTravelTime() { return _lightTravelTime; }

/** Change the light travel time flag. */
void paramSet::lightTravelTime(const bool ltt) { _lightTravelTime = ltt; }

/** Returns a vector of syndyne beta values. */
vector<double> paramSet::beta() { return _beta; }

/** Set the array of beta values for syndynes. */
void paramSet::beta(const vector<double> b) {
  _beta = b;
}

/** Return the maximum number of days to integrate. */
int paramSet::nDays() { return _nDays; }

/** Change the maximum number of days to integrate. */
void paramSet::nDays(const int n) { _nDays = n; }

/** Return the number of steps to take. */
int paramSet::steps() { return _steps; }

/** Change the number of steps to take. */
void paramSet::steps(const int n) { _steps = n; }

/** Return the length of the comet orbit in days. */
double paramSet::orbit() { return _orbit; }

/** Change the length of the comet orbit in days. */
void paramSet::orbit(const double n) { _orbit = n; }

/** Return the number of particles in the simulation. */
long paramSet::nParticles() { return _nParticles; }

/** Change the number of particles in the simulation. */
void paramSet::nParticles(const long n) { _nParticles = n; }

 /** Return the particle label format string. */
string paramSet::labelFormat() { return _labelFormat; }

/** Set the particle label format string. */
void paramSet::labelFormat(const string s) { _labelFormat = s; }

/** Set a parameter to a value. */
void paramSet::setParameter(const string parameter, const string value) {
  string str = parameter;
  transform(str.begin(), str.end(), str.begin(), (int(*)(int))toupper);
  int p = str.find_first_not_of(" ");
  if (p == string::npos) return;

  if (str == "PROGRAM") { program(value); return; }
  if (str == "COMET") { comet(value); return; }
  if (str == "JD") { obsDate(StringConv(value).toScalar<double>()); return; }
  if (str == "BETA") { beta(StringConv(value).toVector<double>()); return; }
  if (str == "XYZFILE") { outFile(value); return; }
  if (str == "NPARTICLES") { nParticles(StringConv(value).toScalar<long>()); return; }
  if (str == "NDAYS") { nDays(StringConv(value).toScalar<int>()); return; }
  if (str == "STEPS") { steps(StringConv(value).toScalar<int>()); return; }
  if (str == "PFUNC") { pFunc(value); return; }
  if (str == "ORBIT") { orbit(StringConv(value).toScalar<double>()); return; }
  if (str == "TOL") { tolerance(StringConv(value).toScalar<double>()); return; }
  if (str == "KERNEL") { spkKernel(value.c_str()); return; }	    
  if (str == "PLANETS") { planets(StringConv(value).toScalar<int>()); return; }
  if (str == "PLANETLOOKUP") {
    planetLookUp(StringConv(value).toBool()); return;
  }
  if (str == "CLOSEAPPROACHES") {
    closeApproaches(StringConv(value).toBool()); return;
  }
  if (str == "BOX") { box(StringConv(value).toScalar<double>()); return; }
  if (str == "LTT") { lightTravelTime(StringConv(value).toBool()); return; }
  if (str == "NOLTT") { lightTravelTime(!StringConv(value).toBool()); return; }

  if (str == "LABEL") { labelFormat(value); return; }
  if (str == "SAVE") { unsetSavedData(); setSavedDataList(value); return; }

  // other stuff
  if (str == "DATA") { unsetSavedData(); setSavedDataList(value); return; }
  if (str == "UNITS") return;

  cerr << "I don't understand parameter: " << str << endl;
}

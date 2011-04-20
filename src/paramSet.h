/***************************************************************************
  Copyright (C) 2005-2007,2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__PARAMSET)
#define __PARAMSET 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

using namespace std;

/** The parameter file class.  This holds all parameter file
  information and can read or write parameter files. */
class paramSet {
 public:
  paramSet();
  paramSet(string);
  paramSet(const paramSet&);
  paramSet(istream&);

  // file io methods
  bool loadParameters(const string);
  bool loadParameters(istream&);
  bool saveParameters(const string);
  void writeParameters(ostream&);
  void writeParameters(ostream&, const bool);
  void writeParameters(ostream&, const string, const bool);

  // set parameters via command line or parameter file
  void setParameter(const string, const string);

  // parameter check
  bool isSyndynes();
  bool isMakeComet();
  bool isIntegrateXYZ();

  // Saved parameter methods
  bool isValidData(const string);
  string getDataDescription(const string);
  string getDataName(const string);
  void unsetSavedData();
  bool isSavedData(const string);
  void setSavedData(const string);
  void unsetSavedData(const string);
  void setSavedData(const string, const bool);
  void setSavedDataList(const string);
  string savedData();
  string savedDataDescription();
  string savedDataUnits();

  // Methods common to both syndynes and mccomet
  string program();
  void   program(const string);
  string comet();
  void   comet(const string);
  string spkKernel();
  void   spkKernel(const string);
  double obsDate();
  void   obsDate(const double);
  void   obsDate(const int, const int, double);
  void   obsDate(const string);
  string outFile();
  void   outFile(const string);
  string labelFormat();
  void   labelFormat(const string s);
  string pFunc();
  void   pFunc(const string);
  double tolerance();
  void   tolerance(const double);
  int    planets();
  void   planets(const int);
  bool   planetLookUp();
  void   planetLookUp(const bool);
  bool   closeApproaches();
  void   closeApproaches(const bool);
  double box();
  void   box(const double);
  bool   lightTravelTime();
  void   lightTravelTime(const bool);

  // Syndyne specific methods
  vector<double> beta();
  void   beta(const vector<double>);
  int    nDays();
  void   nDays(const int);
  int    steps();
  void   steps(const int);
  double orbit();
  void   orbit(const double);		

  // Mccomet specific methods
  long   nParticles();
  void   nParticles(const long);

 private:
  // vars
  bool _lightTravelTime, _planetLookUp, _closeApproaches;
  int _planets, _nDays, _steps;
  long _nParticles;
  double _obsDate, _orbit, _tolerance, _box;
  vector<double> _beta;
  string _program, _comet, _outFile, _pFunc, _spkKernel, _labelFormat, _data;
  vector<string> _dataList;
  vector<string> _unitList;
  int _nData;
};

#endif

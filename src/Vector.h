/***************************************************************************
  Copyright (C) 2004,2005,2007,2008 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__MSKVECTOR)
#define __MSKVECTOR 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

using namespace std;

/** A simple vector class for easy manipulation of coordinates and the like. */
class Vector {
 public:
  // constructors
  Vector();
  Vector(const Vector&);
  Vector(double*);
  Vector(double, double, double);

  // overloaded operators
  double& operator[](int);
  double operator[](int) const;
  Vector operator+(Vector);
  Vector operator-(Vector);
  void   operator+=(Vector);
  void   operator-=(Vector);
  Vector operator*(double);  // scalar multiplier
  double operator*(Vector);  // dot product
  Vector operator%(Vector);  // cross product
  friend ostream& operator<<(ostream&, Vector);

  // methods
  double length();
  Vector unit();
  double *dblarr();
  double cosangle(Vector);
  Vector rotate(const Vector, const double);

 private:
  double a[3];
};

#endif

/***************************************************************************

  A 3-element vector representing rectangular coordinates.

  Copyright (C) 2004,2005,2006,2007,2008 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <ostream>
#include <cmath>
#include "Vector.h"

using namespace std;

/** An empty vector constructor. */
Vector::Vector() {
  a[0] = a[1] = a[2] = 0;
}

/** Defines the vector to be a copy of the input vector. */
Vector::Vector(const Vector& b) {
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

/** Defines the vector to be a copy of the input double array. */
Vector::Vector(double *b) {
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

/** Defines the vector to be a copy of the input variables. */
Vector::Vector(double x, double y, double z) {
  a[0] = x;
  a[1] = y;
  a[2] = z;
}

/** Returns the address of the i-th element of the vector. */
double& Vector::operator[](int i) {
  if (i <= 0) return a[0];
  if (i >= 2) return a[2];
  return a[i];
}

/** Returns the i-th element of the constant vector. */
double Vector::operator[](int i) const {
  if (i <= 0) return a[0];
  if (i >= 2) return a[2];
  return a[i];
}

/** Adds two vectors. */
Vector Vector::operator+(Vector b) {
  Vector c = *this;
  c[0] += b[0];
  c[1] += b[1];
  c[2] += b[2];
  return c;
}

/** Subtracts two vectors. */
Vector Vector::operator-(Vector b) {
  Vector c = *this;
  c[0] -= b[0];
  c[1] -= b[1];
  c[2] -= b[2];
  return c;
}

/** Adds a vector to this vector. */
void Vector::operator+=(Vector b) { *this = *this + b; }

/** Subtracts a vector from this vector. */
void Vector::operator-=(Vector b) { *this = *this - b; }

/** Returns the product of this vector and a scalar. */
Vector Vector::operator*(double b) {
  Vector c = *this;
  c[0] *= b;
  c[1] *= b;
  c[2] *= b;
  return c;
}

/** Returns the dot product of two vectors. */
double Vector::operator*(Vector b) {
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

/** Returns the cross product of two vectors. */
Vector Vector::operator%(Vector b) {
  Vector c;
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return c;
}

/** Writes the vector to a stream. */
ostream& operator<<(ostream& os, Vector b) {
  return os << b[0] << " " << b[1] << " " << b[2];
}

/** Returns the length of the vector. */
double Vector::length() { return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]); }

/** Returns the vector's unit vector. */
Vector Vector::unit() {
  double len = length();
  Vector b = *this;
  b[0] /= len;
  b[1] /= len;
  b[2] /= len;
  return b;
}

/** Returns a three element double array with the vector's coordinates. */
double *Vector::dblarr() { return a; }

/** Returns the cos of the angle measured from a (this vector) to
    b. */
double Vector::cosangle(Vector b) {
  return unit() * b.unit();
}

/** Returns the vector rotated about another (Goldstein 2nd ed.,
    p. 165).  th is in radians. */
Vector Vector::rotate(const Vector n, const double th) {
  Vector rp;  // r prime
  Vector r0 = *this;
  Vector nhat = n;
  nhat = nhat.unit();

  rp = r0 * cos(-th) + (nhat * (nhat * r0)) * (1.0 - cos(-th)) + 
    (r0 % nhat) * sin(-th);
  return rp;
}

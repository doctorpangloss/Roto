/* 

Copyright (C) 2004, Aseem Agarwala, roto@agarwala.org

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

*/



#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#ifdef __APPLE__
using namespace std;
#endif
#include <math.h> 
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

// Karen's CGD

#define MAX_STEPS 1000

// Matrix class the solver will accept
class implicitMatrix
{
 public:
  virtual void matVecMult(const double x[], double r[]) const = 0;
  virtual int numVar() const = 0;
  virtual ~implicitMatrix() {}
};

// Matrix class the solver will accept
// Ax=r, so x is data
class implicitMatrixWithTrans : public implicitMatrix
{
 public:
  virtual void matVecMult(const double x[], double r[]) const = 0;
  virtual void matTransVecMult(double x[], double r[]) const = 0;
  virtual ~implicitMatrixWithTrans() = 0;
};



// Solve Ax = b for a symmetric, positive definite matrix A
// A is represented implicitely by the function "matVecMult"
// which performs a matrix vector multiple Av and places result in r
// "n" is the length of the vectors x and b
// "epsilon" is the error tolerance
// "steps", as passed, is the maximum number of steps, or 0 (implying MAX_STEPS)
// Upon completion, "steps" contains the number of iterations taken
double ConjGrad(int n, implicitMatrix *A, double x[], const double b[], 
		double epsilon,	// how low should we go?
		int    *steps);

// Some vector helper functions
void vecAddEqual(const int n, double r[], const double v[]);
void vecDiffEqual(const int n, double r[], const double v[]);
void vecAssign(const int n, double v1[], const double v2[]);
void vecTimesScalar(const int n, double v[], const double s);
double vecDot(const int n, const double v1[], const double v2[]);
double vecSqrLen(const int n, const double v[]);
double vecMax(const int n, const double v[]);
double vecAbsMax(const int n, const double v[]);
void sanityCheck(const int n, const double v[]);
void vecPrint(const int n, const double v[]);

//----------------------------------------------------------------------


inline void vecAddEqual(const int n, double r[], const double v[])
{
  for (int i = 0; i < n; ++i)
    r[i] = r[i] + v[i];
}

inline void vecDiffEqual(const int n, double r[], const double v[])
{
  for (int i = 0; i < n; ++i)
    r[i] = r[i] - v[i];
}

inline void vecAssign(const int n, double v1[], const double v2[])
{
  //for (int i = 0; i < n; ++i)
  //v1[i] = v2[i];
  memcpy(v1,v2,n*sizeof(double));
}

inline void vecTimesScalar(const int n, double v[], const double s)
{
  for (int i = 0; i < n; ++i)
    v[i] *= s;
}

inline double vecDot(int n, const double v1[], const double v2[])
{
  double dot = 0;
  for (int i = 0; i < n; ++i)
    dot += v1[i] * v2[i];
  return dot;
}

inline double vecSqrLen(const int n, const double v[])
{
  return vecDot(n, v, v);
}

inline double vecMax(const int n, const double v[]) {
  double max = DBL_MIN;
  for (int i=0; i<n; ++i)
    if (v[i] > max) max = v[i];
  return max;
}

inline double vecAbsMax(const int n, const double v[]) {
  double max = DBL_MIN, tmp;
  for (int i=0; i<n; ++i)
    if ((tmp=fabs(v[i])) > max) max = tmp;
  return max;
}

inline void sanityCheck(const int n, const double v[]) {
  for (int i=0; i<n; ++i)
    assert(fabs(v[i]) < 100);
}

inline void vecPrint(const int n, const double v[]) {
  for (int i=0; i<n; ++i)
    printf("%.4f ", v[i]);
  printf("\n");
}

#endif

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



#ifndef MYSPARSEMAT_H
#define MYSPARSEMAT_H

#include <stdlib.h>
#include <assert.h>
#include "linearSolver.h"
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#ifdef __APPLE__
#define _A _A
#define _T _T
#endif

// SPEED: could eventually get rid of _A, pack directly

class MySparseMat {

 public:
  MySparseMat() { _A=NULL; _packed=false;}

  ~MySparseMat();

  bool packed() const { return _packed; }
  

  void init(const int nv, const bool T);

  void pack();

  void printExpanded() {
    int i, j, count;
    for (i=0; i<_nv; i++) {
      count = MIN(10, _nv-_ind[i]);
      for (j=0; j<_ind[i]; j++)
	printf("%7.2f ",0.);
      for (j=0; j<count; j++)
	printf("%7.2f ",_pckdA[i*10 + j]);
      for (j=_ind[i]+count; j<_nv; j++)
	printf("%7.2f ",0.);
      printf("\n");
    }
  }
  
  void drop(double* here, int dim) {
    int i, j, count;
    double* out = here;
    for (i=0; i<_nv; i++) {
      count = MIN(10, _nv-_ind[i]);
      for (j=0; j<count; j++)
	*(out+_ind[i]+j) = _pckdA[i*10 + j];
      out += dim;
    }
  }

  void dropT(double* here, int dim) {
    int i, j, count;
    double* out = here;
    for (i=0; i<_nv; i++) {
      count = MIN(10, _nv-_ind[i]);
      for (j=0; j<count; j++)
	*(out+_ind[i]+j) = _pckdAT[i*10 + j];
      out += dim;
    }
  }

  void addDiagonal(const float f) {
    double *A=_A;
    for (int i=0; i<_nv; i++, A+=_nv+1)
      (*A) += f;
  }
  
  void addnbyk(const int offset, const double coeff, const double* adder, const int n, const int k) {
    int i,j, shop, khop;
    double* A = _A + offset*_nv + offset;
    for (j=0, shop=0, khop=0; j<n; ++j, shop+=_nv, khop+=k) 
      for (i=0; i<k; ++i)
	A[shop + i] += coeff*adder[khop+i];
  }

  void addnbyk(const int roffset, const int coffset, const double coeff, const double* adder, const int n, const int k) {
    int i,j, shop, khop;
    double* A = _A + roffset*_nv + coffset;
    for (j=0, shop=0, khop=0; j<n; ++j, shop+=_nv, khop+=k) 
      for (i=0; i<k; ++i)
	A[shop + i] += coeff*adder[khop+i];
  }

  void mult(const double* x, double* r) const;

  void multT(const double* x, double* r) const;


 private:

  double *_A;
  double *_pckdA;
  double  *_pckdAT;
  
  int *_ind;
  int _nv; // number of variables
  bool _T;
  bool _packed;

};


class MyCongMat : public implicitMatrix  {

 public:

  MyCongMat(MySparseMat *A_, MySparseMat *A_l, int numVec, int l2t) :
    _A(A_), _Al(A_l), _numVec(numVec), l2(l2t) {
    for (int t=0; t<_numVec; ++t)
      assert(_A[t].packed());
  }

  //virtual ~MyCongMat() { delete[] _A; delete[] _Al; }

  void matVecMult(const double x[], double r[]) const;
  
  void printExpanded(FILE* fp);

  int numVar() const { return _numVec*l2; }

 private:
  MySparseMat *_A, *_Al;
  int _numVec, l2;
};




#endif

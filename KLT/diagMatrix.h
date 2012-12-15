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



#ifndef DIAGMATRIX_H
#define DIAGMATRIX_H

// SPEED: consider using C-arrays instead of Ubv
#include <vector>
#include <stdio.h>

// square, symmetric, diagonal matrix
// Only upper half stored, so j >= i is required of all acesses


#include <boost/numeric/ublas/vector.hpp>
typedef boost::numeric::ublas::vector<double> Ubv;
typedef boost::numeric::ublas::zero_vector<double> Zero_Ubv;
typedef unsigned int uint;


class DiagMatrix {

 public:

  // maxWidth is number of off-center diagonals (12)
  DiagMatrix(uint nv, uint maxWidth);

  void matVecMult(const double* x, double* r) const;

  void operator+=(const DiagMatrix& o);

  void operator*=(const double cons);

  const double& operator () (const uint i, const uint j) const;
  
  double& operator () (const uint i, const uint j);

  void takeF(const int i, const int j, const double f);
  // unknown ordering of i,j
  void takeUnrdF(const int i, const int j, const double f);

  int dim() const { return _nv; }

  void outputMat(FILE* fp, int offset=0);

  ~DiagMatrix();

  void clear();

  // assumes j>=i, does not check
  // also checks if within maxWidth, so DOES NOT really mean withinBounds of square matrix
  bool withinBounds(const uint i, const uint j) const;

 private:

  // methods
  void allocateDiag(const uint i);

  // data
  unsigned int _nv, _maxWidth;
  Ubv* _diags; // nv of these, though many are empty, others should be reserved full-length
  

};

inline bool DiagMatrix::withinBounds(const uint i, const uint j) const {
  return (i<_nv && j<_nv && j-i<_maxWidth);
}

inline void DiagMatrix::takeF(const int i, const int j, const double f) {
  (*this)(i,j) += f;
}


// unknown ordering of i,j
inline void DiagMatrix::takeUnrdF(const int i, const int j, const double f) {
  if (j>=i)
    takeF(i,j,f);
  else
    takeF(j,i,f);
}



#endif

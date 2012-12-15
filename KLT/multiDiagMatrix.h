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



#ifndef MULTIDIAGMATRIX_H
#define MULTIDIAGMATRIX_H


#include <boost/numeric/ublas/vector_sparse.hpp>
typedef boost::numeric::ublas::compressed_vector<double> Ubcv;


#include "diagMatrix.h"

class MultiDiagMatrix {

 public:

  MultiDiagMatrix(int numBlocks, const std::vector<int>& maxWidths, const std::vector<int>& dims);

  void matVecMult(const double* x, double* r) const;

  void operator+=(const MultiDiagMatrix& o);

  void operator*=(const double cons);

  const double& operator () (const uint i, const uint j) const;
  
  double& operator () (const uint i, const uint j);

  int dim() const { return _nv; }

  void takeF(const int i, const int j, const double f);
  // unknown ordering of i,j
  void takeUnrdF(const int i, const int j, const double f);

  void outputMat(FILE* fp);

  ~MultiDiagMatrix();

  void clear();

  
 private:

  // methods

  int convertToBlock(uint& bi, uint& bj, const uint i, const uint j) const;

  // data
  uint _numBlocks;
  std::vector<int> _maxWidths;
  std::vector<int> _dims;
  uint _nv; // total size of matrix
  std::vector<DiagMatrix*> _blocks;

  uint* _varToBlock;
  std::vector<uint> _blockToStartVar;
  Ubcv* _junk;  // _nv of these
};


inline int MultiDiagMatrix::convertToBlock(uint& bi, uint& bj, const uint i, const uint j) const { 
  assert(i<_nv && j<_nv);
  int blockNum = _varToBlock[i];
  int offset = _blockToStartVar[ blockNum ];
  bi = i - offset;
  bj = j - offset;  
  return blockNum;
}

inline void MultiDiagMatrix::takeUnrdF(const int i, const int j, const double f) {
  if (j>=i)
    takeF(i,j,f);
  else
    takeF(j,i,f);
}


#endif

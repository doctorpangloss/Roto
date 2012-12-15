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



#include "multiDiagMatrix.h"
#include "myassert.h"

MultiDiagMatrix::MultiDiagMatrix(int numBlocks, const std::vector<int>& maxWidths, const std::vector<int>& dims) : 
_numBlocks(numBlocks), _maxWidths(maxWidths), _dims(dims) {
  uint i;
  assert(_maxWidths.size() == _numBlocks && _dims.size() == _numBlocks);
  _blocks.reserve(numBlocks);
  _blockToStartVar.reserve(numBlocks);

  _nv = 0;
  for (i=0; i<_numBlocks; ++i) {
    _nv += _dims[i];
    _blocks[i] = new DiagMatrix(_dims[i], _maxWidths[i]);    
    if (i==0)
      _blockToStartVar.push_back(0);
    else
      _blockToStartVar.push_back( _blockToStartVar[i-1] + _dims[i-1] );
  }

  // create and fill lookup-table
  _varToBlock = new uint[_nv];
  uint index = 0;
  for (i=0; i<numBlocks; ++i)
    for (uint j=0; j<_dims[i]; ++j, ++index)
      _varToBlock[index] = i;
  assert(index==_nv);

  _junk = new Ubcv[_nv];
  for (i=0; i<_nv; ++i)
    _junk[i].resize(_nv-i);


}

void MultiDiagMatrix::matVecMult(const double* x, double* r) const {
  const double *myx = x;
  double *myr = r;
  uint i;
  for (i=0; i<_numBlocks; ++i) { 
    _blocks[i]->matVecMult(myx,myr);
    myx += _dims[i];
    myr += _dims[i];
  }

  
  // multiply junk
  Ubcv::const_iterator j, vi_end;
  for (i=0; i<_nv; ++i) {
    const Ubcv& vi = _junk[i];
    vi_end = vi.end();
    for (j = vi.begin(); j != vi_end; ++j) {
      assert(j.index() > 0);
      assert(j.index() + i < _nv);
      r[i]             += *j * x[j.index() + i];
      r[j.index() + i] += *j * x[i];   // reflection
    }
  }
  
}

void MultiDiagMatrix::operator+=(const MultiDiagMatrix& o) {
  uint i;
  assert(o._numBlocks == _numBlocks && o._nv == _nv);
  for (i=0; i<_numBlocks; ++i)
    *(_blocks[i]) += *(o._blocks[i]);
      
  for (i=0; i<_nv; ++i)
    _junk[i] += o._junk[i];
}

void MultiDiagMatrix::operator*=(const double cons) {
  uint i;
  for (i=0; i<_numBlocks; ++i)
    *(_blocks[i]) *= cons;

  // mult junk
  for (i=0; i<_nv; ++i)
    _junk[i] *= cons;
}

const double& MultiDiagMatrix::operator () (const uint i, const uint j) const {
  assert(j>=i);
  uint bi, bj, blockNum;
  blockNum = convertToBlock(bi,bj,i,j);
  const DiagMatrix& block = *(_blocks[blockNum]);
  if (block.withinBounds(bi,bj)) {  // part of diagonals
    return  block(bi,bj);
  }
  else {  // not part of diagonals
    return _junk[i](j-i);    
  }
}

  
double& MultiDiagMatrix::operator () (const uint i, const uint j) {
  assert(0);
  assert(j>=i);
  uint bi, bj, blockNum;
  blockNum = convertToBlock(bi,bj,i,j);
  DiagMatrix& block = *(_blocks[blockNum]);
  if (block.withinBounds(bi,bj)) {  // part of diagonals
    return  block(bi,bj);
  }
  else {  // not part of diagonals
    //printf("%d %d\n",i,j);
    return _junk[i](j-i);    
  }
}

void MultiDiagMatrix::takeF(const int i, const int j, const double f) {
  
  assert(j>=i);
  uint bi, bj, blockNum;
  blockNum = convertToBlock(bi,bj,i,j);
  DiagMatrix& block = *(_blocks[blockNum]);
  if (block.withinBounds(bi,bj)) {  // part of diagonals
    block(bi,bj) += f;
  }
  else {  // not part of diagonals
    //printf("%d %d\n",i,j);
    _junk[i](j-i) += f;    
  }
}

void MultiDiagMatrix::outputMat(FILE* fp) {
  uint i;
  for (i=0; i<_numBlocks; ++i) 
    _blocks[i]->outputMat(fp, _blockToStartVar[i]);

  
  // output junk
  Ubcv::const_iterator j;
  for (i=0; i<_nv; ++i) {
    const Ubcv& vi = _junk[i];
    for (j = vi.begin(); j != vi.end(); ++j) {
      fprintf(fp,"%d %d %.10e\n", i+1, j.index() + i+1, *j);
      fprintf(fp,"%d %d %.10e\n",  j.index() + i+1, i+1,*j);
    }
  }
  
}

MultiDiagMatrix::~MultiDiagMatrix() {
  for (uint i=0; i<_numBlocks; ++i)
    delete _blocks[i];
  delete[] _varToBlock;
  delete[] _junk;
}

void MultiDiagMatrix::clear() {
  uint i;
  for (i=0; i<_numBlocks; ++i)
    _blocks[i]->clear();

  // Add Junk
  for (i=0; i<_nv; ++i)
    _junk[i].clear();

  
  /*
  _junk[0](300) += 10;
  _junk[0](301) += 10;
  _junk[1](300) += 10;
  _junk[1](301) += 10; */
}



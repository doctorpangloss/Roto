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



#include "diagMatrix.h"
#include "myassert.h"

DiagMatrix::DiagMatrix(uint nv, uint width) : _nv(nv), _maxWidth(width) {
  //assert(_maxWidth <= _nv);
  if (_maxWidth > _nv)
    _maxWidth = _nv;
  _diags = new Ubv[_maxWidth];
}

// SPEED:  get rid of brackets, use iterators, pointer de-references
void DiagMatrix::matVecMult(const double* x, double* r) const {
  memset(r,0,_nv*sizeof(double));
  uint di=0, dj, s;   
  // do main diagonal separately because no reflection
  
  /*
  uint s=0;
  const Ubv& d = _diags[di];
  s = d.size();
  for (dj=0; dj<s; ++dj)  
    r[dj] += d[dj] * x[di+dj];
  
  // other diagonals
  for (di=1; di<_maxWidth; ++di) {
    const Ubv& d = _diags[di];
    
    s = d.size();
    for (dj=0; dj<s; ++dj) { 
      r[dj]    += d[dj] * x[di+dj];
      r[di+dj] += d[dj] * x[dj];    // reflection
    }    
  }
  */

  const Ubv& d = _diags[di];
  Ubv::const_iterator itj, dend;
  const double *xptr, *r_xptr;
  double *rptr, *r_rptr;

  //dend = d.end();
  s = d.size(); dj = 0;
  rptr = r;  xptr = x+di;
  for (itj=d.begin(); dj<s; ++itj, ++rptr, ++xptr, ++dj)  
    *rptr += *itj * *xptr;
  
  // other diagonals
  for (di=1; di<_maxWidth; ++di) {
    const Ubv& d = _diags[di];
    
    //dend = d.end();
    s = d.size(); dj = 0;
    rptr = r;      xptr = x+di;
    r_rptr = r+di; r_xptr = x;

    for (itj=d.begin(); dj<s; 
	 ++itj, ++rptr, ++xptr, ++r_rptr, ++r_xptr, ++dj) {
      *rptr   += *itj * *xptr;
      *r_rptr += *itj * *r_xptr;    // reflection
    }    
  }



}

void DiagMatrix::operator+=(const DiagMatrix& o) {
  for (uint i=0; i<_maxWidth; ++i) {
    
    if (o._diags[i].size() == 0) {}
    
    else if (_diags[i].size() == 0) _diags[i] = o._diags[i];
    
    else {
      assert(_diags[i].size() == o._diags[i].size());
      _diags[i] += o._diags[i];
    }
  } 
}

void DiagMatrix::operator*=(const double cons) {
  for (uint i=0; i<_maxWidth; ++i)
    _diags[i] *= cons;
}

// access
const double& DiagMatrix::operator () (const uint i, const uint j) const {
  assert(i<_nv && j<_nv && j>=i);
  uint di = j-i, dj = i;
  assert(di<_maxWidth);
  const Ubv& d = _diags[di];
  return d[dj];
}
  
// modification
double& DiagMatrix::operator () (const uint i, const uint j) {
  assert(i<_nv && j<_nv && j>=i);
  uint di = j-i, dj = i;
  assert(di<_maxWidth);
  Ubv& d = _diags[di];

  if (d.size() == 0)
    allocateDiag(di);

  return d(dj);
}




DiagMatrix::~DiagMatrix() {
  delete[] _diags;
}

void DiagMatrix::clear() {
  for (uint i=0; i<_maxWidth; ++i)
    _diags[i].clear();
}

void DiagMatrix::outputMat(FILE* fp, int offset) {
  for (uint i=0; i<_maxWidth; ++i) {
    const Ubv& v = _diags[i];
    for (uint j=0; j<v.size(); ++j) {
      if (v[j] != 0) {
	fprintf(fp,"%d %d %.10e\n",offset + j+1, offset + i+j+1,v[j]);
	if (i!=0) fprintf(fp,"%d %d %.10e\n", offset + i+j+1, offset + j+1,v[j]);
      }
    }
  }
}


//------- PRIVATE  -----------------------------------

void DiagMatrix::allocateDiag(const uint i) {
  assert(i<_maxWidth);
  _diags[i].resize(_nv-i);
}

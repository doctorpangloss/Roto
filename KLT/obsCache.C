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



#include <assert.h>
#include <string.h>
#include "obsCache.h"


ObsCache::ObsCache(SplineKeeper* sk) : _sk(sk) {
  _hits = _totalObs = 0;
  memset(_vars0,0,4*sizeof(int));
  memset(_vars1,0,4*sizeof(int));
  memset(_mat,0,256*sizeof(double));
  memset(_g,0,16*sizeof(double));
  _obsWaiting = 0;
}

void ObsCache::newObs() {
  memset(_J,0,16*sizeof(double));
}

/*
void ObsCache::obsFinished() {
  //++_totalObs;
  //if (_ok)
  //++_hits;
  }*/


bool ObsCache::checkVars(const int v[4], const int i) {
  //return false; 
  if (i==0) {
    if (memcmp(v, _vars0, 4*sizeof(int)) == 0) {   // same
      return true;
    }
    else {
      writeBack();
      memcpy(_vars0,v,4*sizeof(int)); 
    }
      
  }
  else if (i==1) {
    if (memcmp(v, _vars1, 4*sizeof(int)) == 0) {   // same
      return true;
    }
    else {
      writeBack();
      memcpy(_vars1,v,4*sizeof(int));
    }
  }
  else
    assert(0);

  return false;
}

void ObsCache::takeJac(const double J[8], const int ind) {
  int i;
  const double* j = J;

  if (ind==0) {
    for (i=0; i<8; ++i, ++j)
      _J[i] += *j;
  }
  else {
    for (i=8; i<16; ++i, ++j)
      _J[i] += *j;
  }
}

/*
void ObsCache::takeJac(double J1[8], double J2[8],  double r) {
  memcpy(_J,J1,8*sizeof(double));
  memcpy(_J+8,J2,8*sizeof(double));

  int i,j;

  // SPEED: better indicing
  for (i=0; i<16; ++i) {
    _g += _J[i]*r;
    for (j=i; j<16; ++j)
      _mat[i*16 + j] += _J[i]*_J[j];
  }
  
}
*/

void ObsCache::takeJac(const TrackSample& ds, const int ind, const double xcoeff, const double ycoeff) {
  if (ind==0) {
    if (xcoeff != 0) {
      _J[0] += ds._a * xcoeff;
      _J[2] += ds._b * xcoeff;
      _J[4] += ds._c * xcoeff;
      _J[6] += ds._d * xcoeff;
    }
    if (ycoeff != 0) {
      _J[1] += ds._a * ycoeff;
      _J[3] += ds._b * ycoeff;
      _J[5] += ds._c * ycoeff;
      _J[7] += ds._d * ycoeff;
    }
  }
  else if (ind==1) {
    if (xcoeff != 0) {
      _J[8]  += ds._a * xcoeff;
      _J[10] += ds._b * xcoeff;
      _J[14] += ds._d * xcoeff;
      _J[12] += ds._c * xcoeff;
    }
    if (ycoeff != 0) {
      _J[11] += ds._b * ycoeff;
      _J[9]  += ds._a * ycoeff; 
      _J[13] += ds._c * ycoeff;
      _J[15] += ds._d * ycoeff;
    }
  }
}

// SPEED: keyframe or edge terms do not require full matrix, only half
void ObsCache::outer_prod(const double r) { 
  ++_hits;
  ++_totalObs;
  ++_obsWaiting;
  int i,j;
  
  // SPEED: better indicing
  double* rowptr = _mat;
  for (i=0; i<16; ++i, rowptr+=16) {
    _g[i] += _J[i] * r;
    for (j=i; j<16; ++j)
      rowptr[j] += _J[i]*_J[j];
  }
  
}

void ObsCache::outer_prod(const double r, const double c) {
    ++_hits;
  ++_totalObs;
  ++_obsWaiting;
  int i,j;
  
  // SPEED: better indicing
  double* rowptr = _mat;
  for (i=0; i<16; ++i, rowptr+=16) {
    _g[i] += _J[i] * r * c;
    for (j=i; j<16; ++j)
      rowptr[j] += c * _J[i]*_J[j];
  }

}


void ObsCache::writeBack() { 
  if (_obsWaiting == 0) return;
  
  // SPEED: better indicing
  int vars[8];
  memcpy(vars,  _vars0, 4*sizeof(int));
  memcpy(vars+4,_vars1, 4*sizeof(int));

  int i,j;
  double* rowptr = _mat;
  for (i=0; i<16; ++i, rowptr+=16) {
    assert(!isnan(_g[i]));
    _sk->takeg(vars[i/2]+i%2, _g[i]);
    for (j=i; j<16; ++j) {
      if (rowptr[j] != 0)
	_sk->takeF(vars[i/2]+i%2, vars[j/2]+j%2, rowptr[j]); // SPEED
    }
  }

  memset(_mat,0,256*sizeof(double));
  memset(_g,0,16*sizeof(double));
  _obsWaiting = 0;
}

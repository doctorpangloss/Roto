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



#ifndef HB_OneCurve_H
#define HB_OneCurve_H

using namespace std;

#include <stdio.h>
#include <vector>
#include "../jl_vectors.h"
typedef unsigned short ushort;


#define pow2(n) ( 1 << (n))

class HB_OneCurve {

 public:

  // 2^n+1 yields n+1 levels (0->n), where top level has 4 values (corners)
  // so for 17, pass in 5

  // /2 space w,h
  // full-2 fixer, offset
  // x(w) is space, y(h) is time
  HB_OneCurve(int w, int h, vector<int>& fixer, int offset);
  HB_OneCurve() {}
  void init(int w, int h, vector<int>& fixer, int offset);
  ~HB_OneCurve() { delete[] _qu;}


  void sweepUp(double* r) const; // multiply by S'
  
  void sweepDown(double* r) const; // multipy by S
  //void sweepDown2(double* r) const; // multipy by S

  int totalVar() const { return _totalVar; } // full 2 space

  double* extractMatrix(int which) const; // 0 for S, 1 for S'
  
  void printMatrix(double* m, FILE* fp) const;

  void zeroFixed(double* r) const; 

  int space() const { return _w;}
  int time() const { return _h;}

 private:

  void assertFixed(const double* r) const; 
  int die() const{ assert(0); return 1;}
  
  //int _xoffset, _yoffset;
  //int _xlevels, _ylevels, _maxLevels; // # levels, causing dimension 2^(l-1) + 1
  int _w, _h, _maxdim,_totalVar;
  ushort *_qu;
  mutable double _dummy[2];
  vector<int> _fixedLocs;
};



/*
class HB_OneCurve {

 public:

  HB_OneCurve() {}

  void init(int numPts, int numFrames) { // /2 space input
    _numP = numPts; 
    _numF = numFrames;
    _totalVar = _numP*_numF*2;
    _block = new HB_EvenBlock(_numP, _numF);
    printf("%d points, %d frames\n",_numP,_numF);
  }
  ~HB_OneCurve() { delete _block; }

  int totalVar() const { return _totalVar; } // full 2 space
  
  void sweepUp(double* r) const { // multiply by S'
    _block->sweepUp(r);
  }
  
  void sweepDown(double* r) const { // multipy by S
    _block->sweepDown(r);
  }

  double* extractMatrix(int which) const; // 0 for S, 1 for S'

  void printMatrix(double* m, FILE* fp) const;

 private:

  HB_EvenBlock* _block;
  //vector<HB_EvenBlock> _blocks;
  int _numP, _numF, _totalVar; 
};
*/


#endif

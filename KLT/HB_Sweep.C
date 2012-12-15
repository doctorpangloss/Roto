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



#include "HB_Sweep.h"



HB_Sweep::HB_Sweep(int nCurves, int* nPoints, int frames, vector<int>* fixers) {
  _nCurves = nCurves;
  _curves = new HB_OneCurve[_nCurves];
  int offset=0;
  printf("Initing sweepers\n");
  for (int i=0; i<_nCurves; ++i) {
    _curves[i].init(nPoints[i], frames, fixers[i], offset);
    offset+= _curves[i].totalVar();
  }
  printf("Done Initing sweepers\n");
}

void HB_Sweep::sweepUp(double* r) const {
  for (int i=0; i<_nCurves; ++i) {
    if (_curves[i].space()>1 && _curves[i].time()>1)  // skip 1-dim cases
      _curves[i].sweepUp(r);
    r += _curves[i].totalVar();
  }
}

void HB_Sweep::sweepDown(double* r) const {
  for (int i=0; i<_nCurves; ++i) {
    if (_curves[i].space()>1 && _curves[i].time()>1)  // skip 1-dim cases
      _curves[i].sweepDown(r);
    r += _curves[i].totalVar();
  }
}

void HB_Sweep::zeroFixed(double* r) const {
  for (int i=0; i<_nCurves; ++i) {
    _curves[i].zeroFixed(r);
    r += _curves[i].totalVar();
  }
}

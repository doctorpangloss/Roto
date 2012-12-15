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



#ifndef HB_Sweep_H
#define HB_Sweep_H

#ifdef __APPLE__
using namespace std;
#endif
#include <vector>
//#include "multiTrackData.h"
#include "HB_OneCurve.h"

class HB_Sweep {

 public:

  // nPoints is /2 one frame, fixers full-2 variables
  HB_Sweep(int nCurves, int* nPoints, int frames, vector<int>* fixers);
  ~HB_Sweep() { delete[] _curves; }

  void sweepUp(double* r) const; // multiply by S'

  void sweepDown(double* r) const; // multipy by S

  void zeroFixed(double* r) const;

 private:
  
  HB_OneCurve* _curves;
  int _nCurves;
};





#endif

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



#ifndef DRAWCONTCORR_H
#define DRAWCONTCORR_H

#include <set>
#include "rotoPath.h"

class RotoPath;

// maps samples of a drawpath to spline-based track curves

class DrawContCorr {

 public:

  DrawContCorr(int *corrs, RotoPath** corrPtrs, int n);

  DrawContCorr();

  ~DrawContCorr();

  void fillForward(DrawContCorr* o) const;

  const set<RotoPath*>& getRotoSet() const;

  void vacateRotoCorr();

  bool nextRotosExist() const; // checks if rotoset has valid _nextC's

  float getSampleRotoT(const int i);
  //float getRotoT(const float t);
  RotoPath* getSampleRoto(const int i);

 private:

  int _n; // number of discrete samples
  float *_t;  // other(b) T value
  RotoPath** _ptrs; // same length as _samples,_domains, pointer to relevant RotoPath for each sample
  set<RotoPath*> _rotoset;  // set of rotoPaths corresponded to
};



#endif

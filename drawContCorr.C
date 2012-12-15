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



#include "drawContCorr.h"


DrawContCorr::DrawContCorr(int *corrs, RotoPath** corrPtrs, int n) : _n(n) {
  _t = new float[_n];
  _ptrs = new RotoPath*[_n];
  memcpy(_ptrs, corrPtrs, _n*sizeof(RotoPath*));
  for (int i=0; i<_n; ++i) {
    assert(_ptrs[i]);
    _t[i] = _ptrs[i]->getSampleT(corrs[i]);
    _rotoset.insert(_ptrs[i]);
  }
  
}

DrawContCorr::DrawContCorr() {
  _t = NULL;
  _ptrs = NULL;
  _n = 0;
}

DrawContCorr::~DrawContCorr() {
  if (_t) delete[] _t;
  if (_ptrs) delete[] _ptrs;
}

void DrawContCorr::fillForward(DrawContCorr* o) const {
  assert(_t && _ptrs && _n>0);

  o->_n = _n;
  o->_t = new float[_n];
  o->_ptrs = new RotoPath*[_n];
  for (int i=0; i<_n; ++i) {
    o->_ptrs[i] = _ptrs[i]->nextC();
    assert(o->_ptrs[i]);
    o->_t[i] = _ptrs[i]->getNextT(_t[i]);
  }
}

const set<RotoPath*>& DrawContCorr::getRotoSet() const {
  return _rotoset;
}

void DrawContCorr::vacateRotoCorr() {
  assert(_t && _ptrs && _n>0);

  delete[] _ptrs; _ptrs = NULL;
  delete[] _t; _t = NULL;
  _n = 0;
}

bool DrawContCorr::nextRotosExist() const { // checks if rotoset has valid _nextC's
  if (!(_t && _ptrs && _n>0))
    return false;
  for (int i=0; i<_n; ++i)
    if (!_ptrs[i]->nextC())
      return false;
  return true;
}

float DrawContCorr::getSampleRotoT(const int i) {
  assert(_t && _ptrs && _n>0);
  assert(i>=0 && i<_n);
  return _t[i];
}

RotoPath* DrawContCorr::getSampleRoto(const int i) {
  assert(_t && _ptrs && _n>0);
  assert(i>=0 && i<_n);
  return _ptrs[i];
}

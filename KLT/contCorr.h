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



#ifndef CONTCORR_H 
#define CONTCORR_H

#include <qdatastream.h>
#include "bezSpline.h"
#include "myassert.h"

// Encodes continuous correspondence between splines a and b

class ContCorr {

 public:

  ContCorr(const int numSegs, const int numSamples);

  // corrs should be same length as a's discrete samples
  ContCorr(BezSpline* a, BezSpline* b, int* corrs);

  ContCorr(QDataStream* fp);

  void saveqt(QDataStream* fp) const;
  
  // returns t value of later curve corresponding to this curve's sample i
  float getCorrSample(int i) const {
    assert(!_identity);
    assert(i>=0 && i<_n);
    return _samples[i];
  }

  float getT(const float t) const; // gets corresponding T


  int getNumSamples() const { return _n; }
    
  bool identity() const { return _identity; }
  

  ~ContCorr() {
    /*if (_samples) {
      delete[] _samples; delete[] _domain;
      }*/
    //if (_corrs) delete[] _corrs;
  }

  void splitSamples(const float t);
  
  /*
  const int* getDiscreteCorrs() const { return _corrs; }
  int* getDiscreteCorrs() { return _corrs; }
  int getNumDiscreteCorrs() const { return _numSamples; }
  void setNumDiscreteCorrs(int i) { 
    if (_corrs) delete[] _corrs;
    _numSamples = i; _corrs = new int[_numSamples]; }
  */

 private:
  ContCorr(const ContCorr& o);
  ContCorr();
  int binSearchT(const float t, bool* exact) const;
  

  std::vector<float> _samples, _domain;  // other(b) T value, our(a) T value, respectively
  int _n;
  float _aTop, _bTop; // maximum t value for curves a & b; for debugging, error checking
  bool _identity;

  //int* _corrs;  // SPACE not necessary
  int _numSamples;
};




#endif

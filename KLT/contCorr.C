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



#include "contCorr.h"

ContCorr::ContCorr(const int numSegs, const int numSamples) : _numSamples(numSamples) {
  _identity = true; _n = -1;
  //_samples = _domain = NULL;
  _aTop = _bTop = numSegs;
  //_corrs = new int[_numSamples];
  //for (int i=0; i<_numSamples; ++i)
  //_corrs[i] = i;
}

// corrs should be same length as a's discrete samples
ContCorr::ContCorr(BezSpline* a, BezSpline* b, int* corrs) {
  _identity = false;
  _n = a->getDiscreteCount();
  assert(_n>0);
  _aTop = a->numSegs();
  _bTop = b->numSegs();
  assert(corrs[_n-1] < b->getDiscreteCount());

  //_samples = new float[_n];
  //_domain = new float[_n];
  _samples = std::vector<float>(_n);
  _domain = std::vector<float>(_n);
  int j=0, i=0, k, jdiff;
  float up, base;

  _numSamples = _n;
  //_corrs = new int[_numSamples];
  //memcpy(_corrs, corrs, _numSamples*sizeof(int));

  while (i<_n) {

    do {
      ++j;
    } while (j<_n && corrs[j]==corrs[i]);

    jdiff = j-i;
    base = b->getDiscreteT(corrs[i]);
    if (j < _n)
      up = b->getDiscreteT(corrs[j]);
    else
      up = base;
    assert(jdiff>0);

    for (k=0; k<jdiff; ++k, ++i) {
      _samples[i] = base + (float(k)/float(jdiff)) * (up-base);  
      assert(_samples[i]>=0 && _samples[i] <= _bTop);
      _domain[i] = a->getDiscreteT(i);
      assert(_domain[i]>=0 && _domain[i] <= _aTop);
      if (i>0) assert(_domain[i]>_domain[i-1]);
    }
     
  }

    /*
    for (i=0; i<_n; ++i) 
      printf("%d: %d, %.5f %.5f\n",i,corrs[i], _domain[i], _samples[i]);    
      exit(0);*/

}

ContCorr::ContCorr(QDataStream* fp) {
  int i, dummy;

  *fp >> i;
  assert(i==0 || i==1);
  _identity = (bool)i;

  if (!_identity) {
    *fp >> dummy;
    _n = dummy;
    //_samples = new float[dummy];
    //_domain = new float[dummy];
    _samples = std::vector<float>(dummy);
    _domain = std::vector<float>(dummy);
    //float d1, d2;
    for (i=0; i<dummy; ++i) {
      *fp >> _samples[i] >> _domain[i];
    }
  }
  else {
    //_samples = _domain = NULL;
    _n = -1;
  }

  //_corrs = NULL; // suspect
  _numSamples = -1;

  *fp >> _aTop >> _bTop;


}

void ContCorr::saveqt(QDataStream* fp) const {
  int dummy = _n, i;

  i = _identity;
  *fp << i;

  if (!_identity) {
    assert(_samples.size()>0 && _domain.size()>0);
    *fp << dummy;
    for (i=0; i<dummy; ++i) 
      *fp << _samples[i] << _domain[i];
  }

  *fp << _aTop << _bTop;
    
}


// gets corresponding T
float ContCorr::getT(const float t) const {  
  assert(t>=0 && t<=_aTop);
  if (_identity)    return t;

  
  // binary search to appropriate interval, form line, project point onto line  
  bool exact;
  int index = binSearchT(t, &exact);
  
  float res;
  if (exact)
    res = _samples[index];
  else {
    // calculate slope
    float m = (_samples[index+1] - _samples[index]) / (_domain[index+1] - _domain[index]);
    res = _samples[index] + m*(t - _domain[index]);
  }
  assert(res>=0 && res<=_bTop);
  assert(finite(res));

  //printf("in %f, out %f\n",t,res);
  return res;
}

int ContCorr::binSearchT(const float t, bool* exact) const {
  int first = 0, last = _n-1, mid=0;
  assert(t >= _domain[first] && t<=_domain[last]+.1);
  if (t>=_domain[last]) {// not well-handled by binary search
    *exact = true;
    return _n-1;  
  }
  *exact = false;

  while (first < last) {

    mid = (first + last) >> 1; // div 2
    if (t==_domain[mid]) {
      *exact = true; break;
    }
    if (t < _domain[mid])
      last = mid;
    else
      first = mid+1;

  }  

  assert(mid>=0 && mid<_n);  

  if (!*exact) {
    //assert(first==last);
    mid = last-1;
    assert(mid<_n-1 && t>_domain[mid] && t<_domain[mid+1]);
  }

  return mid;
}


void ContCorr::splitSamples(const float t) {
  assert(t>=0 && t<=_bTop);
  int bt = (int)t, i;
  float alpha = t - float(bt);  assert(alpha>0);

  if (_identity) {
    assert(_aTop == _bTop);
    _identity = false;
    ++_bTop;

    assert(_samples.empty() && _domain.empty());
    for (i=0; i<=_bTop; ++i) {
      _samples.push_back(i); 
    }
    for (i=0; i<=bt; ++i) 
      _domain.push_back(i);
    _domain.push_back(t);
    for (i=bt+1; i<=_aTop; ++i) 
      _domain.push_back(i);

    assert(_domain.size() == _samples.size());
    _n = _samples.size();
  }
  else {
    assert(0); //STUB
    
  }
}


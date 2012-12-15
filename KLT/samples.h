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



#ifndef SAMPLES_H
#define SAMPLES_H

#include "jl_vectors.h"

// SPEED: lots of extra data

class DSample {
 public:
  DSample() {}

  //DSample(const Vec2f& loc) { _loc = loc; }
  DSample(const Vec2f& loc, const Vec2f& tan, const Vec2f& curv, const float baset, const float t) : 
    _loc(loc),  _tangent(tan), _curvature(curv), _baset(baset), _t(t) {}
  //DSample(const Vec2f& loc, const Vec2f& tan, const float t, const float bt)  :
  //_loc(loc),  _tangent(tan), _t(t), _bt(bt) {}
  //DSample(const Vec2f& loc, const float t) :
  //_loc(loc),  _t(t) {}
  DSample(const DSample& o) {
    _loc = o._loc; _tangent = o._tangent; _baset = o._baset; _curvature = o._curvature; _t = o._t;
  }

  bool operator< (const DSample& o) const { return (_baset + _t) < (o._baset + o._t); }
  bool operator< (const double ot) const { return (_baset + _t) < ot; }
  bool operator< (const float ot) const { return (_baset + _t) < ot; }

  float totT() const { return _baset + _t; }

  Vec2f _loc;
  Vec2f _tangent;  // not normalized
  Vec2f _curvature; // SPEED: get rid of
  float _baset;  // integral value, t is _baset + _t
  float _t; // this t 0<=t <=1
};


class TrackSample: public DSample {
 public:
  TrackSample() {} 
  TrackSample(const DSample& o);  // calculates _k, _nnormal

  void assertState() { assert(_a>=0 && _a<=1 &&
			      _b>=0 && _b<=1 &&
			      _c>=0 && _c<=1 &&
			      _d>=0 && _d<=1); }

  float _k; // _k * _tangent gives normalized tangent
  Vec2f _nnormal; // normalized normal

  // SPEED: get rid of some of these
  double _a, _b, _c, _d,    // normal abcd
  _ap, _bp, _cp, _dp;     // tangent version
  //_app, _bpp, _cpp, _dpp; // curvature
  double _t;  // 0<=t <=1
};





#endif

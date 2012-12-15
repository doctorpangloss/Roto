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



#ifndef MULTIDATA_H
#define MULTIDATA_H

using namespace std;

#include <vector>
#include <qmutex.h>
#include <qwaitcondition.h>
#include "../jl_vectors.h"

struct FixedCurvePoint {
  int _i; // /2 space, 1 curve
  int _frame; // which frame, varies between 0 & mt->_numFrames-2 inclusive
};

typedef vector<FixedCurvePoint> FixedV;


struct CurveRecord {

  CurveRecord() { memset(_vars, -1, 12*sizeof(int)); _edgeMin = NULL;}

  ~CurveRecord() { if (_edgeMin) delete[] _edgeMin; }

  // the first 'var' in the following three are actually curve points, /2 space
  int varToCurve(const int b)  const {
    if (b==0) return _vars[0][0];
    if (b==_nPoints-1) return _vars[3][0];
    return _vars[1][0];
  }

  int varToCurveVar(const int b)  const {
    if (b==0) return _vars[0][1];
    if (b==_nPoints-1) return _vars[3][1];
    return (_vars[1][1] + b - 1);
  }

  int varToGlobalVar(const int b)  const {
    if (b==0) return _vars[0][2];
    if (b==_nPoints-1) return _vars[3][2];
    return (_vars[1][2] + b - 1);
  }


  // ---------------- DATA
  // /2 space    
  // is _vars[2][] redundant?
  int _vars[4][3]; // 0, 1, n-2, n-1,( which curve, which var, which globalVar)
  int _nPoints;  // one frame
  int _nVars;    // one frame
  int _packedStart;  // where this curve's variables begins in _Z (packed)
  double* _edgeMin; // length _nPoints
  bool _useEdges; 
  FixedV _fixedLocs;
};


// -----------------------------------------------------------------------------

struct MultiTrackData {

  MultiTrackData() { _Z = NULL; _z_mutex = new QMutex; _z_wait = new QWaitCondition(); }

  void transferThreadStuff(MultiTrackData* o);

  void transferEdgeMins(MultiTrackData* o);
  
  ~MultiTrackData() {
    if (_Z)
      delete[] _Z; 
    //delete[] _nCurvePoints; delete[] _curve_map; delete[] _cp_varMap; 
    if (_trackWidths) delete[] _trackWidths;
    if (_crs)
      delete[] _crs;
    if (_z_wait) delete _z_wait;
    if (_z_mutex) delete _z_mutex;
  }
  
  // /2 space
  int cpToCurveNvars(const int c, const int cp) const {
    int newc = _crs[c].varToCurve(cp);
    return _crs[newc]._nVars;
  }

  bool begWithRest(const int c) const {
    return (c == _crs[c]._vars[0][0]);
  }

  bool endWithRest(const int c) const {
    return (c == _crs[c]._vars[3][0]);
  }
  
  // /2 space, knows about keyframe data in _Z
  // In the following two c,j are curve, curve point index
  int packed_index(int c, int j, const int t) const {
    _int_packed_transform(c,j);
    return _int_packed_index(c,j) + t*_crs[c]._nVars;
  }

  // doesn't account for time, so add t*_crs[c]._nVars
  int packed_index(int c, int j) const {
    _int_packed_transform(c,j);
    return _int_packed_index(c,j);
  }

  // doesn't account for time, but sets nvars to time offset
  int packed_index_n(int c, int j, int& nvars) const {
    _int_packed_transform(c,j);
    nvars = _crs[c]._nVars;
    return _int_packed_index(c,j);
  }

  int _int_packed_index(const int c, const int j) const {
    return (_crs[c]._packedStart + j);
  }

  // convert curve, curve point to curve (where var is actually located), packed curve variable (ignoring time)
  void _int_packed_transform(int& c, int& j) const {
    if (j==0) {
      j = _crs[c]._vars[0][1];
      c = _crs[c]._vars[0][0];
    }
    else if (j==_crs[c]._nPoints-1) {
      j = _crs[c]._vars[3][1];
      c = _crs[c]._vars[3][0];
    }
    else
      j = j + _crs[c]._vars[1][1] - 1;
  }

  int size_Z() const { return (_numFrames+1)*_nVars; }
  // remember the following two are /2 space
  int totalNumVar() const { return (_numFrames-1)*_nVars; }
  int totalVars(const int c) const { return _crs[c]._nVars*(_numFrames-1); }


  // ---------------- DATA

  Vec2f* _Z;       // contains keyframe points as well, thus (_numFrames+1)*_nVars long
  int _numFrames;  // there should be _numFrames-1 sets of variables
  int _nCurves;
  int _nVars;      // one frame, all curves
  CurveRecord* _crs;
  Vec2i* _trackWidths; // _nCurves of these, wlh and whh
  QMutex* _z_mutex;
  QWaitCondition* _z_wait;



};


#endif

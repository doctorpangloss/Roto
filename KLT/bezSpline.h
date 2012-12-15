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



#ifndef BEZSPLINE_H
#define BEZSPLINE_H

#include <math.h>
#include <vector>
#include "jl_vectors.h"
#include "talkFitCurve.h"
#include "dynArray.h"
#include "cubicEval.h"
#include "samples.h"

/*
class Coord : public Vec2f {  // control points of bezier segments
 public:
  Coord(const float d1, const float d2) : Vec2f(d1,d2) {}
  Coord *prev, *next;
};
*/

enum SamplingScheme { RESAMPLE, REEVALUATE, RESAMPLE_CONSISTENTLY };
// RESAMPLE: Completely resample evenly, possibly different # samples over time
// REEVALUATE: Calculates new positions given current t values
// RESAMPLE_CONSISTENTLY: Only for initing, should follow with reevaluate.  Samples
// first curve evenly, other curves at same t values


class BezSpline { 

 public:
  
  BezSpline(TalkFitCurve* inCurve);
  //BezSpline(TalkFitCurve* inCurve, float sample_spacing, const int* numSamples);
  BezSpline(const BezSpline* other);
  BezSpline();
  ~BezSpline();

  int numSegs() const { return _numSegs; }
  int numCtrls() const { return _ctrlpts.size(); }
  bool empty() const { return _ctrlpts.empty(); }

  void finishBuilding(); // if you added controls manually

  void setControl(const Vec2f&c, const int i);
  const Vec2f& getEnd(const int which) const;
  bool isCtrlEnd(int i, int* whichEnd) const; 
  bool isCtrlAnyEnd(const int i) const;
  double distanceToEnd2(const Vec2f& pt, const int side) const;
  double distanceToInternalEnds2(const Vec2f& pt, int* whichCtrl) const;
  double distanceToCtrls2(const Vec2f& pt, int* whichCtrl) const;
  void setEnd(const Vec2f& c, const int side);

  double distToSamples2(const Vec2f& loc, DSample& sample) const;
  
  void translate(const Vec2f& delta);

  // should only be used in conjunction with empty spline
  void clearControls();
  void addCtrl(const Vec2f& c) { _ctrlpts.push_back(c); }
  void setNumSegs(const int i);

  const Vec2f* getCtrl(const int i) const {
    assert(i>=0 && i< (int)_ctrlpts.size());
    return &(_ctrlpts[i]);
  }

  const std::vector<Vec2f>& getControls() const { return _ctrlpts; }
  std::vector<Vec2f>& getControls() { return _ctrlpts; }

  void renderCtrlPts() const;
  void renderHandles() const;

  void calcEval2(CubicEval2* ce, const Vec2f pts[4]);
  // Calculates _ce evaluator eval, using 4 vecs starting at _ctrlpts[pts] 
  void calcEval2(const int eval, const int pts); 

  void calcAllEval2();

  void clearSamples();

  void redoSamples(float sample_spacing, const int *numSamples);

  void reevaluateSamples();
  
  void redoSamplesConsistently(const BezSpline* ref);

  void calcOneSample(DSample* ds, float t) const;

  void recalcOneSample(DSample* ds);

  Vec2f getLoc(float t) const;
  Vec2f getTangent(float t) const;
  // ----------- 

  std::vector<DSample>& getSamples() { return _samples; }

  const DSample& getDiscreteSample(const int i) const {
    return _samples[i]; }

  DSample& getDiscreteSample(const int i) {
    return _samples[i]; }

  Vec2f getDiscreteLoc(const int i) const {
    return _samples[i]._loc;
  }

  Vec2f getDiscreteTangent(const int i) const {
    return _samples[i]._tangent;
  }

  float getDiscreteT(const int i) const {
    const DSample& s = _samples[i];
    return s._t + s._baset;
  }

  int getDiscreteCount() const {
    return _samples.size();
  }

  void print() const;
  void printSamples() const;

  void buildEvals();

  void split(float t);

  float findClosestT(const Vec2f& loc);
  
  //------ for track management

  void printZMap() const;
  void printVarMap() const;

  void setMapLength(const int c) { _mapLength = c; }
  void createZMap() { assert(_mapLength>0); _ZMap = new int[_mapLength]; }
  void createVarMap() { assert(_mapLength>0); _varMap = new int[_mapLength]; }
  int ZMap(const int i) const { assert(i>=0 && i<_mapLength && _ZMap); return _ZMap[i]; }
  int varMap(const int i)  const { assert(i>=0 && i<_mapLength && _varMap); return _varMap[i]; }
  int& ZMap(const int i) { assert(i>=0 && i<_mapLength && _ZMap); return _ZMap[i]; }
  int& varMap(const int i)  { assert(i>=0 && i<_mapLength && _varMap); return _varMap[i]; }
  bool hasZMap() const { return _ZMap != NULL; }
  bool hasvarMap() const { return _varMap != NULL; }
  void clearMaps() { delete[] _ZMap; delete[] _varMap; _ZMap = _varMap = NULL; _mapLength = 0;}

  int finalVariable() const; // /2
  int beginningVariable() const; // /2


  //------ 

 private:
  // takes evenSamples and uses them to fill _samples
  // Does not add first sample
  void addEvenSamples(const float s, DynArray<float, 200>* evenSamples);


  // fills evenSamples with uneven sampling
  // Does not add last sample
  float sample(const float a, const float b, const CubicEval2& c, 
	       DynArray<float,200>* evenSamples) const;
  bool flat(const float p, const float q, const float r, const CubicEval2& ce) const;


  void addSample(float t, const CubicEval2& ce, const float base); 

  //DynArray<DSample,100> _samples;
  std::vector<DSample> _samples;
  std::vector<Vec2f> _ctrlpts;
  int _numSegs;
  // Why are these two different?  Keyframes are in Z, not variables
  int* _ZMap;   // used by multisplinedata, maps control points to indices in Z (/2)
  int* _varMap;  // used by multisplinedata, maps control points to variables (thus no keyframes!) (/2)
  int _mapLength; // length of ZMap, varMap
  std::vector<CubicEval2> _ce;



};


#endif

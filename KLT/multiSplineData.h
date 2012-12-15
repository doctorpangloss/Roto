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



#ifndef MULTISPLINEDATA_H
#define MULTISPLINEDATA_H

#include <qmutex.h>
#include <vector>
#include <qwaitcondition.h>
#include "bezSpline.h"
#include "contCorr.h"
//#include "contEdgeMin.h"
#include "jl_vectors.h"
#include "samples.h"

class FixedControl {
 public:
  int _n;  // whivh control
  int _c;  // which curve
  int _t;  // which time slice (varies from 1 to _numFrames-1)
};

typedef std::vector<FixedControl> FixedControlV;
typedef std::vector<Vec2f> ZVec;
typedef std::vector<TrackSample> SampleHolder;

class MultiSplineData {
 
 public:

  MultiSplineData();
  void finishInit();

  // fills splines with control points, affects later discretization
  void takeControls(ZVec* Z);

  //------------------ OLD way of sampling
  /*
  // discretize curve c at time t at interval pixel intervals into TrackSamples such that further calls to 
  // getDiscreteCount, getSample, getCorrSample are affected
  void discretize(int t, int c, int interval);

  // return number of discrete samples from discretize
  int getDiscreteCount() const {
    return _ts.size(); }

  // Get TrackSample n, from discretization
  // Put full-2 variable numbers into vars for x-components
  void getDiscreteSample(int n, TrackSample *res, int vars[4]) const;
  void getDiscreteSample(int n, TrackSample *res) const;

  // Also, get corresponding sample one frame later
  void getCorrDiscreteSample(int n, TrackSample *res, int vars[4]) const;
  void getCorrDiscreteSample(int n, TrackSample *res) const;
  */
  //------------------

  // Discretizes all curves at pixel interval.  If consistent, same number of samples for all time
  void discretizeAll(int interval, SamplingScheme ss, bool onlyInbetweens = false);

  // returns number of samples at time t, curve c
  int getNumSamples(int t, int c) const; 

  void getDiscreteSample(int t, int c, int n, TrackSample* res, int vars[4]) const;
  void getDiscreteSample(int t, int c, int n, TrackSample* res) const;

  // Also, get corresponding sample one frame later
  // creates new sample
  void getCorrDiscreteSample(int t, int c, int n, TrackSample *res, int vars[4]) const;
  void getCorrDiscreteSample(int t, int c, int n, TrackSample *res) const;
  
  // maps to existing sample in frame t+1, returns sample number (/2)
  int getExistingCorrDiscreteSample(int t, int c, int n, TrackSample *res, int vars[4]) const;
  int getExistingCorrDiscreteSample(int t, int c, int n, TrackSample *res) const;

  //------------------


  // Given curve parameter tp, time t, curve c, fill TrackSample res and vars 
  void getSample(float tp, int t, int c, TrackSample* res, int vars[4]) const;
  void getSample(float tp, int t, int c, TrackSample* res) const;

  // return number of Vec2f's in Z 
  int size_Z() const { return _Z.size(); }


  // full-2, actual variables
  int numVars() const { return (_nTotalVars << 1); } // *2

  // Z2 = Z + sol, watching out for keyframe data
  // returns maximum manhattan step
  // Z & Z2 cannot be the same!
  double createTestSol(const ZVec* Z, const double* sol, ZVec* Z2) const;
  
  // outputs vector of control variables, ignoring keyframes
  void outputControlVars(char* name) const;

  // REQUIRES identity correspondence between two keyframes (same # controls)
  void lerpControls();  

  ~MultiSplineData();

  // Given time t, curve c, control point n, return index in Z
  int tcn(const int t, const int c, const int n) const;
  // gives tcn for first or last control (depending on side)
  int end_tc(const int t, const int c, const int side) const;
  int end_tc_var(const int t, const int c, const int side) const;
  // Gives /2 var number rather than Z index
  int tcn_var(const int t, const int c, const int n) const;

  // given time t, curve c, give number of Bezier segments in curve
  int tc_numSegs(const int t, const int c) const;
  int tc_numControls(const int t, const int c) const { return tc_numSegs(t,c)*3+1;}
  // gives full-2 number of variables for one curve, all time
  int c_numVars(const int c) const;
  int c_numVarsOneFrame(const int c) const;

  void transferThreadStuff(MultiSplineData* o);
  
  void transferEdgeMins(MultiSplineData* o);

  const FixedControlV& getFixedLocs() const { return _fixedControls; }

  BezSpline* getSpline(const int c, const int t) { 
    assert(c>=0 && c<_nCurves && t>=0 && t<=_numFrames);
    return _splines[c*(_numFrames+1) + t]; }
  const BezSpline* getSpline(const int c, const int t) const { 
    assert(c>=0 && c<_nCurves && t>=0 && t<=_numFrames);
    return _splines[c*(_numFrames+1) + t]; }

  const SampleHolder* getSampleHolder(const int c, const int t) const { 
    assert(c>=0 && c<_nCurves && t>=0 && t<=_numFrames);
    return _holders[c*(_numFrames+1) + t]; }

  bool useEdges(const int c) const;
  std::vector<double>* refreshEdgeMins(const int c, const int num);

  void addMask(unsigned char* mask) { _masks.push_back(mask); }
  void takeMasks(MultiSplineData* o);
  void setMaskDims(const int w, const int h) { _maskw=w; _maskh=h; }
  bool occluded(const Vec2f& loc, const int t) const;
  

  // -------------- DATA ------------------

  Vec2i* _trackWidths; // _nCurves of these, wlh and whh
  int _numFrames;  // there should be _numFrames-1 sets of variables
  int _nCurves;
  int _nTotalVars;   // full-2, actual variables (/2)
  QMutex* _z_mutex;
  QWaitCondition* _z_wait;


  ZVec _Z;       // contains keyframe controls as well
  std::vector<const ContCorr*> _conts; // Corresponders, _numFrames*_nCurves of them, curve-major
  std::vector<BezSpline*> _splines; // (_numFrames+1)*_nCurves of them, curve-major
  std::vector<int> _numSegs;        // (_numFrames+1)*_nCurves of them, curve-major, number of segments in curve
  std::vector<SampleHolder*> _holders;  // (_numFrames+1)*_nCurves of them, curve-major
  FixedControlV _fixedControls;
  
  std::vector<bool> _useEdges;
  std::vector< std::vector<double> > _edgeMins;
  std::vector<unsigned char*> _masks; // should be _numFrames+1 of these
  int _maskw, _maskh;

  //ContEdgeMin* _edgemins; // nCurves of these
  //int _nVars0; // num controls in first keyframe (/2)
  //BezSpline _discretized;
  //std::vector<TrackSample> _ts;
  //int _dis_t, _dis_c, _dis_numSegs; // Discretization's time, curve, number of segments
  //BezSpline* _dis;
  //int* _ZMap; 
  //int _ctnloc;
  //int _ctNumSegs;





};



#endif 

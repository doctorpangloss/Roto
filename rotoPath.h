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



#ifndef ROTOPATH_H
#define ROTOPATH_H

#ifdef __APPLE__
using namespace std;
#endif
#include <set>
#include <algorithm>
#include <vector>
#include <list>
#include <fstream>
#include <qdatastream.h>
#include "jl_vectors.h"
#include "dynArray.h"
#include "rotoCorresponder.h"
#include "imgSequence.h"
#include "joint.h"
#include "trackGraph.h"
#include "KLT/talkFitCurve.h"
#include "KLT/contCorr.h"
#include "rotoRegion.h"

class DrawPath;
class DrawCurves;
class TrackGraph;
class SplineCtrlIterator;
class SplineSampleIterator;
class RotoRegion;

typedef std::vector<RotoPath*> PathV;
typedef std::pair<RotoPath*, RotoPath*> RotoPathPair;

class RotoPath : public AbstractPath {

 public:
  RotoPath();
  RotoPath(const RotoPath& other);
  ~RotoPath();

  //--- manual painting
  
  void addCtrl(const Vec2f& v);
  void startManualCreation();
  void finishManualCreation();
  void startFinishManualCreation();
  //---
  

  float getNextT(float s) const; // t value of next corresponding location
  float getPrevT(float s) const; // t value of next corresponding location
  float getSampleNextT(int s) const; // t value of sample's next corresponding location
  float getSampleT(int s) const;
  const Vec2f getNextLoc(int s) const;
  Vec2f getLoc(float t) const {
    assert(_bez); return _bez->getLoc(t);
  }

  Vec2f getTangent(float t) const {
    assert(_bez); return _bez->getTangent(t);
  }

  int getNumControls() const { assert(_bez); return _bez->numCtrls(); } // /2 space
  int getNumSegs() const { assert(_bez); return _bez->numSegs(); }
  const Vec2f* getCtrl(const int i) { assert(_bez); return _bez->getCtrl(i); }
  const Vec2f& getCtrlRef(const int i) { assert(_bez); return *(_bez->getCtrl(i)); }
  const BezSpline* getBez() const { return _bez; }
  bool isCtrlEnd(int i, int* whichEnd) const { assert(_bez); return _bez->isCtrlEnd(i,whichEnd); }
  const Vec2f& getEnd(const int which) const { assert(_bez->numCtrls()>0); return _bez->getEnd(which); }
  void setEnd(const Vec2f& c, const int side) { assert(_bez && _bez->numCtrls()>0); _bez->setEnd(c,side); }
  double distanceToSide2(const Vec2f pt, const int side) const {
    assert(_bez); return _bez->distanceToEnd2(pt,side); }
  double distanceToEnds2(const Vec2f pt) {
    assert(_bez); return MIN(_bez->distanceToEnd2(pt,0), _bez->distanceToEnd2(pt,1)); }
  double distanceToInternalEnds2(const Vec2f pt, int* whichCtrl) {
    assert(_bez); return _bez->distanceToInternalEnds2(pt, whichCtrl); }
  double distanceToCtrls2(const Vec2f pt, int* whichCtrl) {
    assert(_bez); return _bez->distanceToCtrls2(pt, whichCtrl); }
  void setControl(const Vec2f& c, const int i) { assert(_bez); _bez->setControl(c,i); } 
  void translate(const Vec2f& delta);

  double distanceToEnd(const Vec2f pt) const { assert(0); return 0;} 
  double distanceToStart(const Vec2f pt) const { assert(0); return 0;}

  double distToSplineSamples2(const Vec2f& loc, DSample& sample) {
    assert(_bez); return _bez->distToSamples2(loc,sample); }

  void takeBez(const BezSpline* b);

  void handleNewBezCtrls();

  float findClosestT(const Vec2f& loc);

  RotoPathPair split(const int c);

  void addVertex(const double x, const double y) {
    Vec2f p(x,y);
    if (count==0 || p != getElement(count-1)) add(p);
  }

  void addVertex(const Vec2f& v) {
    add(v); }

  void startPrevCorrespondence(RotoPath* prevC);
  void finishPrevCorrespondence();
  void clearPrevCorrespondence();

  bool assertState() const;

  //Vec2f getNextCorrPoint(const int i) const; // returns point in _prevFrame that
  // corresponds to OUR point i
  
  //int getNextCorrIndex(const int i) const; // return index in _nextFrame that corresponds to OUR point i
  //int getPrevCorrIndex(const int i) const; // return index in _prevFrame that corresponds to OUR point i

  void addCorrDraw(DrawPath* dp) const {
    _corrDraws.insert(dp); }

  const set<DrawPath*>* getCorrDraws() const { return &_corrDraws; }
  set<DrawPath*>* getCorrDraws() { return &_corrDraws; }

  int getNumCorrDraws() const { return _corrDraws.size(); }
  void emptyCorrDraws() { _corrDraws.clear(); }

  void getUV(const float i, Vec2f& U, Vec2f& V) const; 
  // unit tangent and perp vector of rotoPath, now done using bez

  void render(bool showTrackPoints, bool visEffort=false) const;
  void renderHandles() const { if (_bez) _bez->renderHandles(); };
  void renderHandle(const int i);
  void renderNextCorrLines() const;

  //void renderNextCorrPS(FILE* fp) const;

  void renderSelected() const;

  void renderTransformSelected();
  void transformSelect(int t);
  void clearTransformSelected() { _transformSelected.clear(); }
  void translateSelected(const Vec2f& delta);
  void uniformScaleSelectedAbout(float mag, const Vec2f& center);
  void rotateSelectedAbout(float dtheta, const Vec2f center);

  SplineCtrlIterator transformSelectedCtrlIterator();
  SplineSampleIterator transformSelectedSampleIterator();

  //void renderEmphasizedEndpoint(int which) const;

  void doublePoints();

  void spliceIn(const RotoPath *other, const int i1, const int i2); // inclusive range

  void save(FILE* fp, int version) const;
  void saveqt(QDataStream * fp, int version) const;
  void load(FILE* fp, int version);
  void loadqt(QDataStream* fp, int version);
  void dumpPoints(FILE* fp) const;

  void saveXML(std::ofstream& fp) const;
  void labelSelfXML(int& labelNum);
  void clearXMLLabels() { _xmlLabel = -1;}

  void notifyDrawDelete(DrawPath* dp);

  void calculateBackCorrs(); // calculates bi-directional correspondence between _prevFrame and this

  void handleCorrForw(RotoPath* other);

  void goBezier(const float tolSqr, const int* numSamples);

  bool okEnds() const;

  RotoPath* nextC() { return _nextFrame; }
  const RotoPath* nextC() const { return _nextFrame; }
  void setNextC(RotoPath* next) {
    _nextFrame = next; }
  RotoPath* prevC() { return _prevFrame; }
  const RotoPath* prevC() const { return _prevFrame; }
  void setPrevC(RotoPath* prev) {
    _prevFrame = prev; }
  //bool nextCoors() { return (_nextCoors != NULL); }
  bool nextCorrs() { return (_nextCont != NULL && _nextFrame != NULL); }
  void takeNextCont(RotoPath* o) { _nextCont = o->_nextCont; o->_nextCont = NULL; }

  RotoPath *next, *prev;

  void buildCan(bool first);
  void setCan(int t);

  void calcLerp();

  //void setLastTranslate(const float dx, const float dy) {
  // _lastTranslation.Set(dx,dy);
  //}

  //Vec2f getLastTranslation() const {
  //return _lastTranslation; }

  void actualTranslate(const int n, const Vec2f& t);
  void actualSetLoc(const int n, const Vec2f& l);

  //void setTrans(const double* trans);
  //const double* getTrans() const { return _trans; }
  //const void getTrans(double* here) const { memcpy(here,_trans,6*sizeof(double)); }
  
  void buildINextCorrs(); // build increasing next coors
  void buildIPrevCorrs(); // build increasing prev coors
  //  void generateInverseCoors(const int* inCoors, const int inNum, 
  //		    int *outCoors, const int outNum);

  int lowHeight() const { return _lh;}
  int highHeight() const { return _hh; }
  void setLowHeight(const int h) { _lh=h;}
  void setHighHeight(const int h) { _hh=h;}
  int trackWidth() const { return _trackw; }
  //int trackSubsampling() const { return _trackSubsampling; }
  //void setTrackSubsampling(const int s);
  bool samplingDirty() const { return _samplingDirty; }

  void setFixed(const bool b) { _fixed=b; }
  bool fixed() { return _fixed; }
  void flipFixed();

  void buildTouched(const bool val);
  void setCtrlTouched(const int which, const bool val);
  void setTouchedAll(const bool val);
  void estimateTouched();
  bool touched(const int i) const;
  bool endTouched(const int side) const;
  bool jointTouched(const int side) const;
  Vec2i getStats() const;

  int pickCtrl(const float x, const float y);

  void nudge(const int i, const Vec2f& loc);

  void printInfo();

  void vacateRotoCorr();
  void vacateNextRotoCorr();
  void vacatePrevRotoCorr();

  void buildccomp(TrackGraph* out, const PathV* toTrack);

  bool trackEdges() const { return _trackEdges; }
  void flipTrackEdges() { _trackEdges = !_trackEdges; }
  void setTrackEdges(const bool i) { _trackEdges = i; }

  void insertJoint(const Joint& j, int w);
  void destroyJoints();
  void destroyEndJoints();
  void destroyBegJoints();
  void destroyAllCoincidentJoints(const int s);
  void replaceMeInJoints(RotoPath* newrp, const int side);
  void removeJoint(RotoPath* rp, const int side); // remove joints referencing rp
  const Joints* bJoints() const { return &_bJoints;}
  const Joints* eJoints() const { return &_eJoints;}
  Joints* bJoints() { return &_bJoints;}
  Joints* eJoints() { return &_eJoints;}
  Joints& joints(const int side);
  const Joints& joints(const int side) const;
  void reconcileJoints(const PathV* flexible);
  void reconcileOneJointToMe(int which);
  void reconcileJointToMe(Joints* js, const Vec2f* rploc);
  void addToJointSet(JointSet& all, const int side);

  static void printJoints(const Joints* j);
  void copyBackJoints(const PathV* subtree);
  void copyForwardJoints(const PathV* subtree);

  void fixEndpoint(int which);
  void fixInternal(int which);
  bool endFixed(const int which) const {   assert(which==0 || which==1); return _endsFixed[which]; }
  bool anyEndFixed() const { return _endsFixed[0] || _endsFixed[1]; }
  bool internalFixed(const int which) const;
  bool somehowFixed(const int which) const;
  void unfixCtrl(const int which);

  ContCorr* getNextCont() { return _nextCont;}

  static void copySetToVector(vector<RotoPath*>* v, set<RotoPath*>* s);

  void fillFromBez();

  void splitSegment(float t);

  bool allJointsOk();

  // non-inclusive range, does not modify paths directly (but maybe objects pointed to
  static void buildBackReconcileJoints(PathV& paths, int bFrame, int aFrame);
  static void buildForwardReconcileJoints(PathV& paths, int aFrame, int bFrame);

  void setRegion(RotoRegion *r, const short side) { _regions[side] = r; }
  RotoRegion* getRegion(const short side) { return _regions[side]; }
  void forgetRegion(const RotoRegion* r);

  void setTrimapWidth(const int i) { _trimapWidth = i; }

  //int *_prevCoors, *_nextCoors; // should be nondecreasing

  BezSpline* _bez;

 private:
  //void buildGLTransMatrix() const;
  void mutualInit();
  void h(const int n, const float i, const float j, const int l, const float invc,
	 const Vec2f* P, Vec2f& result) const;
  void hnew(const int n, const float i, const float j, 
	    const Vec2f* P, Vec2f& result) const;
  void hnorm(const int n, const float i, const float j, const int l, const Vec2f* P, Vec2f& result) const;
  int buildBack(TrackGraph* out, const PathV* toTrack, const int side);
  int buildForw(TrackGraph* out, const PathV* toTrack, const int side);
  void reconcileJoint(Joints* js, const int side, const PathV* flexible);
  void copyBackJoint(Joints* js, const PathV* subtree, int side);
  void copyForwardJoint(Joints* js, const PathV* subtree, int side);
  void reconcileFixEndpointJoint(const int side);

  mutable set<DrawPath*> _corrDraws;
  RotoPath *_prevFrame, *_nextFrame;
  bool _fixed;
  
  int _trimapWidth;
  int _lh, _hh; // low and high heights for rotoPath tracking.  _lh is negative, _hh positive
  int _trackw; // tracking width, used symmetrically
  //int _trackSubsampling;
  //Vec2f _lastTranslation;
  //Vec2f* _subLocs;
  //int _nsubLocs;
  bool _trackEdges;
  bool _samplingDirty; // whether or not the points in this array sample _bez's  current control points

  bool _endsFixed[2];  
  std::set<int> _internalsFixed; // fixed controls internal to curve
  std::vector<bool> _touched;

  Joints _bJoints, _eJoints; // begining and end joints
  ContCorr *_nextCont, *_prevCont; 

  int _xmlLabel;

  RotoRegion* _regions[2]; // TODO: write to file!

  std::list<int> _transformSelected; // ranges (segments), maintained sorted
};


class SplineCtrlIterator {
 public:

  SplineCtrlIterator() {} // starts at beginning
  bool end() const;
  SplineCtrlIterator& operator++();
  const Vec2f& ctrl() const;
  Vec2f& ctrl();

  std::vector<Vec2f>* _ctrls;
  const std::list<int>* _selected;
  uint _i;
  std::list<int>::const_iterator _j;
};

class SplineSampleIterator {
 public:
  SplineSampleIterator() {} // starts at beginning
  bool end() const;
  SplineSampleIterator& operator++();
  const DSample& sample() const;
  DSample& sample();

  std::vector<DSample>* _samples;
  const std::list<int>* _selected;
  int _i;
  std::list<int>::const_iterator _j;
};

#endif

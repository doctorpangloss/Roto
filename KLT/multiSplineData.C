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



#include <algorithm>
#include "multiSplineData.h"
#include "myassert.h"

MultiSplineData::MultiSplineData() {
  //_edgemins = NULL;
  _z_mutex = new QMutex; _z_wait = new QWaitCondition(); 
  _maskw = _maskh = 0;
}

void MultiSplineData::finishInit() {
  //_nVars0 = tcn(1,0,0); // number of controls in first keyframe
  //_nTotalVars = _Z.size() - _nVars0 - (_Z.size() - tcn(_numFrames, 0, 0));
  const BezSpline* b = getSpline(_nCurves-1, _numFrames-1);  
  _nTotalVars = b->finalVariable() + 1;
  printf("half 2 number of variables %d\n", _nTotalVars);
  _holders.reserve(_splines.size());

  int t,c;

  for (c=0; c<_nCurves; ++c) {
    for (t=0; t<=_numFrames; ++t) {
      _holders.push_back(new SampleHolder());
    }
  }

  for (c=0; c<_nCurves; ++c) {
    printf("Curve %d: \n",c);
    for (t=0; t<=_numFrames; ++t) {
      printf("   Time %d\n", t);
      const BezSpline* bez = getSpline(c,t);
      printf("   Spline %x\n",bez);
      assert(bez);
      bez->printZMap();
      if (t>0 && t<_numFrames)
      bez->printVarMap();      
    }
  }

  //lerpControls(); // only for manual mode

}


  // fills splines with control points, affects later discretization
void MultiSplineData::takeControls(ZVec* Z) {
  int i, nctrls, numSegs, c, t;
  for (c=0; c<_nCurves; ++c) {
    for (t=0; t<=_numFrames; ++t) {
      numSegs = tc_numSegs(t,c);  assert(numSegs > 0);
      BezSpline* bez = getSpline(c,t);
      bez->clearControls();
      nctrls = 3*numSegs+1;
      bez->setNumSegs(numSegs);
      for (i=0; i<nctrls; ++i)
	bez->addCtrl((*Z)[tcn(t,c,i)]);
      bez->calcAllEval2();      
    }
  }

}

// HERE CHANGE
double MultiSplineData::createTestSol(const ZVec* Z, const double* sol, ZVec* Z2) const { 
  assert(Z != Z2);
  int c, t, numC, varIndex, zIndex, n;
  double maxStep = 0;
  Vec2f v;
  for (c=0; c<_nCurves; ++c) {
    for (t=1; t<_numFrames; ++t) {   // don't do keyframes
      numC = tc_numControls(t,c);
      for (n=0; n<numC; ++n) {
	zIndex = tcn(t,c,n);
	varIndex = 2*tcn_var(t,c,n);
	v.Set(sol[varIndex], sol[varIndex+1]);
	maxStep = MAX(MAX(maxStep, fabs(sol[varIndex])), fabs(sol[varIndex+1]));
	Vec2f_Add((*Z2)[zIndex], v, (*Z)[zIndex]);
      }
	
    }
  }

  assert(finite(maxStep));
  return maxStep;

}

// outputs vector of control variables, ignoring keyframes
void MultiSplineData::outputControlVars(char* name) const {
  FILE* fp = fopen(name,"w");
  int c, t, numC, zIndex, n;
  //double maxStep = 0;
  Vec2f v;
  for (c=0; c<_nCurves; ++c) {
    for (t=1; t<_numFrames; ++t) {   // don't do keyframes
      numC = tc_numControls(t,c);
      for (n=0; n<numC; ++n) {
	zIndex = tcn(t,c,n);
	fprintf(fp, "%.10f\n%.10f\n", _Z[zIndex].x(), _Z[zIndex].y());
      }
    }
  }
  fclose(fp);
}



/*
// discretize curve c at time t at interval pixel intervals into TrackSamples such that further calls to 
// getDiscreteCount, getSample, getCorrSample are affected
void MultiSplineData::discretize(int t, int c, int interval) {

  _dis_t = t; _dis_c = c;
  _dis_numSegs = tc_numSegs(t,c);
  assert(interval > 0);

  // discretized.clear();
//   int i, nctrls;
//   nctrls = 3*_dis_numSegs+1;
//   _discretized.setNumSegs(_dis_numSegs);
//   for (i=0; i<nctrls; ++i)
//     _discretized.addCtrl(Z[tcn(t,c,i)]);

//   _discretized.calcAllEval2();
//   _discretized.redoSamples(1.f,NULL);

//   _ts.clear();
//   for (i=0; i<_discretized.getDiscreteCount(); ++i) {
//     TrackSample newts(_discretized.getDiscreteSample(i));
//     _ts.push_back(newts);
//   }
  

  _dis = getSpline(c,t);
  _dis->redoSamples(float(interval), NULL);
  _ts.clear();
  for (int i=0; i<_dis->getDiscreteCount(); ++i) {
    TrackSample newts(_dis->getDiscreteSample(i));
    _ts.push_back(newts);
  }
  printf("Discretized time %d, curve %d, to %d samples at interval %d\n",t,c,_dis->getDiscreteCount(),interval);
}
*/


/*
// Get TrackSample n, from discretization
void MultiSplineData::getDiscreteSample(int n, TrackSample *res) const { 
  assert(n<_ts.size());
  *res = _ts[n];
}


// Get TrackSample n, from discretization
// Put full-2 variable numbers into vars for x-components
void MultiSplineData::getDiscreteSample(int n, TrackSample *res, int vars[4]) const { 
  assert(n<_ts.size());
  *res = _ts[n];
  int base = (int) res->_baset;

  base *= 3; int k;
  for (int i=0; i<4; ++i) {
    k = tcn_var(_dis_t, _dis_c, base+i); 
    assert(k < _nTotalVars);
    vars[i] = k << 1; // k*2
  }
}

// Also, get corresponding sample one frame later
void MultiSplineData::getCorrDiscreteSample(int n, TrackSample *res, int vars[4]) const {
  assert(_dis_t < _numFrames); // can't do this on last keyframe
  float ot = _conts[_dis_t*_nCurves + _dis_c]->getT(_ts[n]._t + _ts[n]._baset);
  getSample(ot, _dis_t+1, _dis_c, res, vars);
}

void MultiSplineData::getCorrDiscreteSample(int n, TrackSample *res) const {
  assert(_dis_t < _numFrames); // can't do this on last keyframe
  float ot = _conts[_dis_t*_nCurves + _dis_c]->getT(_ts[n]._t + _ts[n]._baset);
  getSample(ot, _dis_t+1, _dis_c, res);
}
*/
//------------------



// Discretizes all curves at pixel interval.  If consistent, same number of samples for all time
void MultiSplineData::discretizeAll(int interval, SamplingScheme ss, bool onlyInbetweens) {
  int i=0, j, c, t;
  for (c=0; c<_nCurves; ++c) {
    for (t=0; t<=_numFrames; ++t, ++i) { // need <=
      if ((t==0 || t==_numFrames) && onlyInbetweens) continue;
      BezSpline* bs = _splines[i];
      SampleHolder* sh = _holders[i];
      sh->clear();


      if (ss==RESAMPLE)
	bs->redoSamples(float(interval), NULL);
      else if (ss==REEVALUATE)
	bs->reevaluateSamples();
      else if (ss==RESAMPLE_CONSISTENTLY) {
	if (t==0/* || t==_numFrames*/) {   // SWITCH
	  bs->redoSamples(float(interval), NULL);
	  //if (t==0) {
	  //bs->print(); 
	  //}
	}
	else 
	  bs->redoSamplesConsistently(getSpline(c,0));
      
      }
      else
	assert(0);

      //printf("Curve %d, time %d,  has %d samples\n",c,t, bs->getDiscreteCount()); 
      sh->reserve(bs->getDiscreteCount());
      //if (t==0 || t==_numFrames) printf("curve %d time %d discretize to %d samples, scheme %d\n",c,t,bs->getDiscreteCount(), ss);
      // SPEED: for reevaluate, recalculating some things that stay consistent
      for (j=0; j<bs->getDiscreteCount(); ++j) { // SPEED: too many layers to samples
	TrackSample newts(bs->getDiscreteSample(j));
	sh->push_back(newts);
      }
    }
  }
}
  
// returns number of samples at time t, curve c
int MultiSplineData::getNumSamples(int t, int c) const {
  return getSampleHolder(c,t)->size();
}

void MultiSplineData::getDiscreteSample(int t, int c, int n, TrackSample* res, int vars[4]) const {
  assert(t>0  && t<_numFrames);  // no variables first/last frame
  assert(c>=0 && c<_nCurves);
  const SampleHolder* sh = getSampleHolder(c,t);
  assert(n>=0 && n<(int)sh->size());
  *res = (*sh)[n];

  int base = (int) res->_baset;
  base *= 3; int k;
  for (int i=0; i<4; ++i) {
    k = tcn_var(t, c, base+i); 
    assert(k < _nTotalVars);
    vars[i] = k << 1; //2*k;
  }
}

void MultiSplineData::getDiscreteSample(int t, int c, int n, TrackSample* res) const {
  assert(t>=0  && t<=_numFrames);
  assert(c>=0 && c<_nCurves);
  const SampleHolder* sh = getSampleHolder(c,t);
  assert(n>=0 && n<(int)sh->size());
  *res = (*sh)[n];
}

void MultiSplineData::getCorrDiscreteSample(int t, int c, int n, TrackSample *res, int vars[4]) const {
  assert(t>=0  && t<_numFrames-1);  // no variables last frame
  assert(c>=0 && c<_nCurves);

  const SampleHolder& sh = *(getSampleHolder(c,t));
  assert(n>=0 && n<(int)sh.size());
  float ot = _conts[c*_numFrames + t]->getT(sh[n]._t + sh[n]._baset);
  getSample(ot, t+1, c, res, vars);   // creates new sample
}

void MultiSplineData::getCorrDiscreteSample(int t, int c, int n, TrackSample *res) const {
  assert(t>=0  && t<_numFrames);
  assert(c>=0 && c<_nCurves);

  const SampleHolder& sh = *(getSampleHolder(c,t));
  assert(n>=0 && n<(int)sh.size());
  float ot = _conts[c*_numFrames + t]->getT(sh[n]._t + sh[n]._baset);
  getSample(ot, t+1, c, res);  // creates new sample
}

int MultiSplineData::getExistingCorrDiscreteSample(int t, int c, int n, TrackSample *res, int vars[4]) const { 
  assert(t>=0  && t<_numFrames-1);  // no variables last frame
  assert(c>=0 && c<_nCurves);
  //if (t != _numFrames-1) {        // should be identity
  getDiscreteSample(t+1,c,n,res,vars); return n;
    //}

}

int MultiSplineData::getExistingCorrDiscreteSample(int t, int c, int n, TrackSample *res) const { 

  assert(t>=0  && t<_numFrames);
  assert(c>=0 && c<_nCurves);
  if (t != _numFrames-1) {        // should be identity
    getDiscreteSample(t+1,c,n,res); return n;
  }

  getCorrDiscreteSample(t,c,n,res); return -1; // G!

  // SPEED: binary search each time, could pre-compute
  const SampleHolder& mysh = *(getSampleHolder(c,t));
  assert(n>=0 && n<(int)mysh.size());
  float myt = mysh[n]._t + mysh[n]._baset;
  float ot = _conts[c*_numFrames + t]->getT(myt);
  const SampleHolder& latersh = *(getSampleHolder(c,t+1));
  SampleHolder::const_iterator found = std::lower_bound(latersh.begin(), latersh.end(), ot);					  
  
  int index;
  if (found == latersh.begin())
    index = 0;
  else if (found==latersh.end())
    index = latersh.size()-1;
  else {
    float ut =  (*found).totT(), lt = (*(found-1)).totT();
    if (fabs( ut - ot ) < fabs( lt - ot ))
      index = found - latersh.begin();
    else
      index = found - 1 - latersh.begin();
  }
  //printf("getExistingCorrDiscreteSample: myt: %f, ot: %f, n: %d, found %d,  choose %d\n",
  // myt, ot, n, found-latersh.begin(), index); 
    
  getDiscreteSample(t+1,c,index,res);
  return index;
}


//------------------

// tp is baset + t
void MultiSplineData::getSample(float tp, int t, int c, TrackSample* res) const {
  const BezSpline* bez = getSpline(c,t);

  DSample ds;
  bez->calcOneSample(&ds, tp);
  *res = ds;
}

// Given curve parameter tp, time t, curve c, fill TrackSample res and vars 
void MultiSplineData::getSample(float tp, int t, int c, TrackSample* res, int vars[4]) const {
  const BezSpline* bez = getSpline(c,t);

  DSample ds;

  bez->calcOneSample(&ds, tp);
  int base = (int) ds._baset;
  *res = ds;

  base *= 3; int k;
  for (int i=0; i<4; ++i) {
    k = tcn_var(t, c, base+i);
    assert(k < _nTotalVars);
    vars[i] = k << 1; //2*k;
  }
}


MultiSplineData::~MultiSplineData() {
  uint i;
  //for (i=0; i<_conts.size(); ++i)
  //delete _conts[i];
  for (i=0; i<_splines.size(); ++i)
    delete _splines[i];  
  for (i=0; i<_holders.size(); ++i) 
    delete _holders[i];
  for (i=0; i<_masks.size(); ++i)
    delete[] (_masks[i]);

  //if (_edgemins) delete[] _edgemins;
  //for (i=0; i<_edgeMins.size(); ++i)
  //if (_edgeMins[i]) delete[] _edgeMins[i];

  if (_trackWidths) delete[] _trackWidths;
  if (_z_wait) delete _z_wait;
  if (_z_mutex) delete _z_mutex;
}

// CAREFUL! Not same as variable number
// Given time t, curve c, control point n, return index in Z  (/2)
int MultiSplineData::tcn(const int t, const int c, const int n) const {  
  // make sure to check validity of numbers
  //assert(t<=_numFrames);
  //assert(c<_nCurves);
  const BezSpline* b = getSpline(c,t);
  assert(b->hasZMap());
  assert(n<=tc_numSegs(t,c)*3);
  return b->ZMap(n);
}

// Gives /2 var number rather than Z index
int MultiSplineData::tcn_var(const int t, const int c, const int n) const {
  assert(t>0 && t<_numFrames);
  //return (tcn(t,c,n) - _nVars0);
  //assert(c<_nCurves);
  assert(n<=tc_numSegs(t,c)*3);
  const BezSpline* b = getSpline(c,t);  
  assert(b->hasvarMap());
  return b->varMap(n);
}

int MultiSplineData::end_tc(const int t, const int c, const int side) const {
  if (side==0)
    return tcn(t,c,0);
  else
    return tcn(t,c,tc_numControls(t,c)-1);
}

int MultiSplineData::end_tc_var(const int t, const int c, const int side) const {
  if (side==0)
    return tcn_var(t,c,0);
  else
    return tcn_var(t,c,tc_numControls(t,c)-1);
}

// given time t, curve c, give number of Bezier segments in curve
int MultiSplineData::tc_numSegs(const int t, const int c) const {  
  // make sure to check validity of numbers
  assert(t<=_numFrames);
  assert(c<_nCurves);
  //return _splines[c*_numFrames + t]->numSegs();
  return _numSegs[c*(_numFrames+1) + t];
}

int MultiSplineData::c_numVars(const int c) const {
  const BezSpline *last = getSpline(c,_numFrames-1), *first = getSpline(c,1);
  return 2*(last->finalVariable() - first->beginningVariable() + 1);
}

int MultiSplineData::c_numVarsOneFrame(const int c) const {
  const BezSpline *first = getSpline(c,1);
  return 2*(first->finalVariable() - first->beginningVariable() + 1);
}


void MultiSplineData::transferThreadStuff(MultiSplineData* o) {
  delete _z_mutex; _z_mutex = o->_z_mutex; o->_z_mutex = NULL;
  delete _z_wait; _z_wait = o->_z_wait;    o->_z_wait = NULL;
}

void MultiSplineData::takeMasks(MultiSplineData* o) {
  _masks = o->_masks;
  _maskw = o->_maskw;
  _maskh = o->_maskh;
  assert(_masks.size() == _numFrames + 1);
  o->_masks.clear();
}

bool MultiSplineData::occluded(const Vec2f& loc, const int t) const {
  int x = (int) loc.x(), y = (int) loc.y();
  if (x <= 0 || x >= _maskw ||
      y <= 0 || y >= _maskh)
    return false;  // out of bounds is NOT occluded
  assert(t>=0 && t<_masks.size());
  if (_masks[t] == NULL)
    return false;

  return ((_masks[t])[y*_maskw + x]);
}

void MultiSplineData::transferEdgeMins(MultiSplineData* o) {
  _edgeMins = o->_edgeMins;
  _useEdges = o->_useEdges;
  /*
  for (int i=0; i<_nCurves; ++i) {
    _edgemins = o->_edgemins;
    o->_edgemins = NULL;
    }*/
}

bool MultiSplineData::useEdges(const int c) const {
  return _useEdges[c];
}

std::vector<double>* MultiSplineData::refreshEdgeMins(const int c, const int num) {
  //if (_edgeMins[c] != NULL)
    //delete[] _edgeMins[c];
    
    //_edgeMins[c] = new double[num];
  _edgeMins[c].clear();
  _edgeMins[c].reserve(num);
  return &(_edgeMins[c]);}


void MultiSplineData::lerpControls() {
  int c, t, numC, n;
  for (c=0; c<_nCurves; ++c) { // iterate over curves
    
    numC = tc_numControls(0,c);
    assert(numC == tc_numControls(_numFrames,c));
    for (n=0; n<numC; ++n) {  // iterate over controls

      Vec2f a = _Z[tcn(0,          c, n)];
      Vec2f b = _Z[tcn(_numFrames, c, n)];

      for (t=1; t<_numFrames; ++t) {  // iterate over non-keyframe time
	assert(numC == tc_numControls(t,c));
	
	float p = float(t) / float (_numFrames);
	Vec2f_LinInterp(_Z[tcn(t, c, n)], a, b, p);

      } // t
    } // n
  } // c
}


TrackSample::TrackSample(const DSample& o) : DSample(o) { 
  _nnormal.Set(-_tangent.y(), _tangent.x());
  _k = 1. / (_nnormal.Normalize());
  //assert(finite(_k) && !isnan(_k));

  double t, mt;
  _t = o._t;
  assert(_t>=0 && _t<=1); 
  t = _t;
  mt = 1.-t;
  
  _a = mt*mt*mt;
  _b = 3.*t*mt*mt;
  _c = 3.*t*t*mt;
  _d = t*t*t;

  
  _ap = -3.*mt*mt;
  _bp = 9.*t*t-12.*t+3.;
  _cp = 3.*t*(-3.*t+2.);
  _dp = 3.*t*t;
  
  //double speed = sqrt(_tangent.x()*_tangent.x() + _tangent.y() + _tangent.y());
  
  /*
  _app = 6.*mt;  //_app /= speed;
  _bpp = 18.*t - 12.;  //_bpp /= speed;
  _cpp = -18.*t + 6.; //_cpp /= speed;
  _dpp = 6.*t;  //_dpp /= speed;
  */

  //assertState(); 
}

/*
void SamplePreComp::compute(const DSample& s) {
  double t, mt;
  _t = s._t;
  assert(_t>=0 && _t<=1); 
  t = _t;
  mt = 1.-t;
  
  _a = mt*mt*mt;
  _b = 3.*t*mt*mt;
  _c = 3.*t*t*mt;
  _d = t*t*t;
  
  _ap = -3.*mt*mt;
  _bp = 9.*t*t-12.*t+3.;
  _cp = 3.*t*(-3.*t+2.);
  _dp = 3.*t*t;
  
  _app = 6.*mt;
  _bpp = 18.*t - 12.;
  _cpp = -18.*t + 6.;
  _dpp = 6.*t;
}
*/

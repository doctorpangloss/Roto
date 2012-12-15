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



#include "bezSpline.h"
#include <GL/gl.h>
#include <float.h>



BezSpline::BezSpline(TalkFitCurve* inCurve) {
  
  _numSegs = inCurve->getSize();
  _ZMap = NULL;
  _varMap = NULL;
  _mapLength = 0;

  inCurve->SetIterationHead();
  int i;
  const Bez* currBez = inCurve->IterateNext();
  assert(currBez);

  Vec2f newC(currBez->_data[0].x(), currBez->_data[0].y());
  _ctrlpts.push_back(newC);

  _ce.reserve(inCurve->getSize());
  // _ce = new CubicEval[inCurve->getSize()*2];

  int j=0;
  do {
    for (i=1; i<4; i++) {
      newC.Set(currBez->_data[i].x(), currBez->_data[i].y());
      _ctrlpts.push_back(newC);
    }
    CubicEval2 newce;
    calcEval2(&newce,currBez->_data);
    _ce.push_back(newce);
    j++;
  }  while ((currBez = inCurve->IterateNext()) != NULL);


    // CORRECTNESS: internal endpoints should average left and right tangents

  assert(_ctrlpts.size()/3 == _numSegs);
  assert(_ctrlpts.size()%3 == 1);
}


/*
BezSpline::BezSpline(TalkFitCurve* inCurve,float sample_spacing, const int * numSamples) {
  
  _numSegs = inCurve->getSize();
  _ZMap = NULL;
  _varMap = NULL;

  inCurve->SetIterationHead();
  int i,j;
  const Bez* currBez = inCurve->IterateNext();
  assert(currBez);

  Vec2f newC(currBez->_data[0].x(), currBez->_data[0].y());
  _ctrlpts.push_back(newC);
  Vec2f headTan(3.*(currBez->_data[1].x() - currBez->_data[0].x()),
		3.*(currBez->_data[1].y() - currBez->_data[0].y()));
  //headTan.Normalize();
  Vec2f headCurv (6.*currBez->_data[0].x() - 12.*currBez->_data[1].x() + 6.*currBez->_data[2].x(),
		  6.*currBez->_data[0].y() - 12.*currBez->_data[1].y() + 6.*currBez->_data[2].y());
  _samples.add(DSample(newC, headTan, headCurv, 0., 0.)); // these are the evenly spaced samples

  //_ce = new CubicEval[inCurve->getSize()*2];
  _ce.reserve(inCurve->getSize());

  // these are the uneven samples -- go figure....
  // These are just t values....
  DynArray<float, 200> evenSamples;        
  float L = 0, n;
  j=0;
  do {
    for (i=1; i<4; i++) {
      newC.Set(currBez->_data[i].x(), currBez->_data[i].y());
      _ctrlpts.push_back(newC);
    }
    CubicEval2 newce;
    calcEval2(&newce,currBez->_data);
    _ce.push_back(newce);
    L += sample(0,1,newce,&evenSamples);
    j++;
  }
  while ((currBez = inCurve->IterateNext()) != NULL);
  evenSamples.add(1);

  if (numSamples==NULL)
    n = ceil(L/sample_spacing);
  else
    n = float(*numSamples-1);
  float s = L / n; 
  assert(finite(s));
  printf("target distance %f\n",s);
  inCurve->SetIterationHead();
  addEvenSamples(s, &evenSamples);

    // CORRECTNESS: internal endpoints should average left and right tangents

  assert(int(_ctrlpts.size())/3 == _numSegs);
  assert(_ctrlpts.size()%3 == 1);
  assert(numSamples==NULL || *numSamples==_samples.getNumElements());
  //printSamples();
}
*/

BezSpline::BezSpline(const BezSpline* other) {
  _samples = other->_samples;
  _ctrlpts = other->_ctrlpts;
  _numSegs = other->_numSegs;
  _ce = other->_ce;
  assert((int)other->_ce.size() == _numSegs);
  _ZMap = NULL;
  _varMap = NULL;
  _mapLength = 0;
}


BezSpline::BezSpline() { 
  _numSegs=0; 
  _ZMap = NULL;
  _varMap = NULL;
  _mapLength = 0;
}

void BezSpline::finishBuilding() {
  int i = _ctrlpts.size();
  int ci = ((i-1)/3) * 3 + 1;
  if (ci!=i)
    _ctrlpts.resize(ci);
  
  _numSegs = (ci-1)/3;
  buildEvals();
}


void BezSpline::setNumSegs(const int i) {
  _numSegs = i;
  _ctrlpts.reserve(3*i  + 1);
  _ce.reserve(i);
  buildEvals();
}

void BezSpline::clearControls() {
  //clearSamples();
  _numSegs = 0;
  _ce.clear();
  _ctrlpts.clear();
  
}

void BezSpline::setControl(const Vec2f& c, const int i) {
  assert(i>=0 && i<_ctrlpts.size());
  _ctrlpts[i] = c;
}


const Vec2f& BezSpline::getEnd(const int which) const {
  assert(_ctrlpts.size() >=2);
  if (which==0)
    return _ctrlpts[0];
  else if (which==1)
    return _ctrlpts[_ctrlpts.size()-1];
  else {
    assert(0);
    return _ctrlpts[0];
  }
}

void BezSpline::setEnd(const Vec2f& c, const int which) {
  assert(_ctrlpts.size() >=2);
  if (which==0)
    _ctrlpts[0] = c;
  else if (which==1)
    _ctrlpts[_ctrlpts.size()-1] = c;
  else 
    assert(0);
}

bool BezSpline::isCtrlEnd(int i, int* whichEnd) const {
  assert(i>=0 && i<_ctrlpts.size());
  if (i==0) {
    *whichEnd = 0;
    return true;
  }
  if (uint(i)==_ctrlpts.size()-1) {
    *whichEnd = 1;
    return true;
  }
  return false;
}

bool BezSpline::isCtrlAnyEnd(const int i) const {
  assert(i>=0 && i<_ctrlpts.size());
  return (i%3==0);
}

double BezSpline::distanceToEnd2(const Vec2f& pt, const int side) const {
  return getEnd(side).distanceTo2(pt);
}

double BezSpline::distanceToInternalEnds2(const Vec2f& pt, int* whichCtrl) const {
  // iterate over controls, skipping first, last, and internal shape controls
  double res = DBL_MAX;
  for (uint i=3; i<_ctrlpts.size()-1; i+=3) {
    float d = _ctrlpts[i].distanceTo2(pt);
    if (d<res) {
      res = d;
      if (whichCtrl)
	*whichCtrl = i;
    }
  }
  return res;
}

double BezSpline::distanceToCtrls2(const Vec2f& pt, int* whichCtrl) const {
  // iterate over controls
  double res = DBL_MAX;
  for (uint i=0; i<_ctrlpts.size(); ++i) {
    float d = _ctrlpts[i].distanceTo2(pt);
    if (d<res) {
      res = d;
      if (whichCtrl)
	*whichCtrl = i;
    }
  }
  return res;
}

double BezSpline::distToSamples2(const Vec2f& loc, DSample& sample) const {
  assert(_samples.size() > 0);
  double res = DBL_MAX;
  int which=-1;
  for (uint i=0; i<_samples.size(); ++i) {
    float d = _samples[i]._loc.distanceTo2(loc);
    if (d<res) {
      res = d;
      which = i;
    }
  }

  sample = _samples[which];
  return res;
}



void BezSpline::redoSamplesConsistently(const BezSpline* ref) {
  assert(ref->_numSegs == _numSegs);
  _samples.clear();
  for (uint i=0; i<ref->_samples.size(); ++i) {
    DSample ds  = ref->_samples[i];
    recalcOneSample(&ds);
    _samples.push_back(ds);
  }
  
  assert(_samples.size() == ref->_samples.size());
}


void BezSpline::reevaluateSamples() {
  DSample* ds;
  assert(_numSegs > 0);
  assert(_samples.size() > 0);
  for (uint i=0; i<_samples.size(); ++i) {
    ds = &(_samples[i]);
    recalcOneSample(ds);
  }
}


void BezSpline::redoSamples(float sample_spacing, const int *numSamples) {
  //printf("Redo samples\n");

  int i;
  if (_numSegs==0) return;
  assert(_numSegs > 0);
  assert(_ctrlpts.size()/3 == _numSegs);
  assert(_ctrlpts.size()%3 == 1);

  Vec2f headTan(3.*(_ctrlpts[1].x() - _ctrlpts[0].x()),
		3.*(_ctrlpts[1].y() - _ctrlpts[0].y()));
  Vec2f headCurv (6.*_ctrlpts[0].x() - 12.*_ctrlpts[1].x() + 6.*_ctrlpts[2].x(),
		  6.*_ctrlpts[0].y() - 12.*_ctrlpts[1].y() + 6.*_ctrlpts[2].y());
  //headTan.Normalize();
  _samples.clear();
  _samples.push_back(DSample(_ctrlpts[0], headTan, headCurv, 0, 0));


  DynArray<float, 200> evenSamples;
  float L = 0, n;
  for (i=0; i<_numSegs; ++i) 
    L += sample(0,1,_ce[i],&evenSamples);
  evenSamples.add(1);

  if (numSamples==NULL)
    n = ceil(L/sample_spacing);
  else
    n = float(*numSamples-1);
  float s = L / n; 
  assert(finite(s));
  //printf("target distance %f\n",s);
  addEvenSamples(s, &evenSamples);

    // CORRECTNESS: internal endpoints should average left and right tangents

  assert(int(_ctrlpts.size())/3 == _numSegs);
  assert(int(_ctrlpts.size())%3 == 1);
  assert(numSamples==NULL || *numSamples==_samples.size());
  //printSamples();
}



void BezSpline::calcEval2(CubicEval2 *ce, const Vec2f pts[4]) {
  ce->_x.init(3.*(pts[1].x() - pts[2].x()) - pts[0].x() + pts[3].x(),
	     3.*(pts[0].x() + pts[2].x()) - 6.*pts[1].x(),
	     3.*(pts[1].x() - pts[0].x()), pts[0].x());
  ce->_y.init(3.*(pts[1].y() - pts[2].y()) - pts[0].y() + pts[3].y(),
	     3.*(pts[0].y() + pts[2].y()) - 6.*pts[1].y(),
	     3.*(pts[1].y() - pts[0].y()), pts[0].y());
}

void BezSpline::calcEval2(const int eval, const int pts) {

  _ce[eval]._x.init(3.*(_ctrlpts[pts+1].x() - _ctrlpts[pts+2].x()) - _ctrlpts[pts].x() + _ctrlpts[pts+3].x(),
	     3.*(_ctrlpts[pts].x() + _ctrlpts[pts+2].x()) - 6.*_ctrlpts[pts+1].x(),
	     3.*(_ctrlpts[pts+1].x() - _ctrlpts[pts].x()), _ctrlpts[pts].x());

  _ce[eval]._y.init(3.*(_ctrlpts[pts+1].y() - _ctrlpts[pts+2].y()) - _ctrlpts[pts].y() + _ctrlpts[pts+3].y(),
	     3.*(_ctrlpts[pts].y() + _ctrlpts[pts+2].y()) - 6.*_ctrlpts[pts+1].y(),
	     3.*(_ctrlpts[pts+1].y() - _ctrlpts[pts].y()), _ctrlpts[pts].y());
}

void BezSpline::calcAllEval2() {
  int i,j;
  assert(_numSegs*3+1 == _ctrlpts.size());
  assert((int)_ce.size() == _numSegs);
  for (i=0, j=0; i<_numSegs; ++i, j+=3) {
    //  printf("%f %f, %f %f, %f %f, %f %f\n",
    //   _ctrlpts[j].x(), _ctrlpts[j].y(),
    //   _ctrlpts[j+1].x(), _ctrlpts[j+1].y(),
    //   _ctrlpts[j+2].x(), _ctrlpts[j+2].y(),
    //   _ctrlpts[j+3].x(), _ctrlpts[j+3].y());
    calcEval2(i, j);
  }
  //printf("\n");
}

void BezSpline::buildEvals() {
  assert(_ce.size()==0);
  for (int i=0; i<_numSegs; i++)
    _ce.push_back(CubicEval2());
}

BezSpline::~BezSpline() {
  if (_ZMap) delete[] _ZMap;
  if (_varMap) delete[] _varMap;
  //if (_ce) delete[] _ce;
}


void BezSpline::renderCtrlPts() const {
  glPointSize(6.0);  // 
  glBegin(GL_POINTS);
  for (uint i=0; i<_ctrlpts.size(); ++i) {
    if (i%3==0)
      glColor3f(1,0,0); // 0,0,1 // G! normally yellow
    else
      glColor3f(0,1,0); // 0,1,0
    const Vec2f& curr = _ctrlpts[i];
    //if (i%3==0) //G!
    glVertex2f(curr.x(), curr.y());
  }
  glEnd();



}

void BezSpline::renderHandles() const {
  glColor3f(0,1,0);
  glLineWidth(1);
  glBegin(GL_LINES);
  for (uint i=0; i<_ctrlpts.size(); ++i) {
    const Vec2f& curr = _ctrlpts[i];
    if (i%3==0 && i!=0)
      glVertex2f(curr.x(), curr.y());
    glVertex2f(curr.x(), curr.y());
  }
  glEnd();
}


void BezSpline::printSamples()  const {
  for (uint i=0; i<_samples.size(); ++i) {
    const DSample& s = _samples[i];
    printf("Sample t %.4f bt %.4f at (%.4f, %.4f), tangent %.4f %.4f, curv %.4f %.4f\n",s._t,s._baset,s._loc.x(),s._loc.y(),s._tangent.x(), s._tangent.y(),
	   s._curvature.x(), s._curvature.y()
	   );
  }
}

void BezSpline::addSample(float t, const CubicEval2& ce, const float base) {
  if (t>1) {
    printf("Bad news %f\n",t);
    t= 1;
  }
  assert(t>=0 && t<=1);
  Vec2f samp(ce._x.f(t), ce._y.f(t));
  Vec2f goodT(ce._x.deriv_1(t), ce._y.deriv_1(t));
  //float speed2  = goodT.x()*goodT.x() + goodT.y()*goodT.y();
  Vec2f curv(ce._x.deriv_2(t), ce._y.deriv_2(t));
  //printf("Curvature %f %f\n",curv.x(), curv.y());
  //printf("Sample t %f bt %f at (%f, %f), tangent %f %f\n",t,base,samp.x(),samp.y(),goodT.x(), goodT.y());
  //goodT.Normalize();
  _samples.push_back(DSample(samp, goodT, curv, base, t));
  //if (_samples.getNumElements()>1) { 
  //DSample v = _samples.getElement(_samples.getNumElements()-2);
  //printf("dist from last %f\n",samp.distanceTo(v._loc));
  //assert(samp.distanceTo(v._loc) > .01);
  //}
}

// from page 301 of Solomon
void BezSpline::addEvenSamples(const float s, DynArray<float, 200>* evenSamples) {

  float st=0, totSegLen=0, t=0, segLen;
  float xl,yl,a,b;
  bool changece=false;
  int k=evenSamples->getNumElements(), i;
  //printf("Starting with %d samples\n",k);
  float base = 0;
  int cenum=0;
  for (i=1; i<k; i++) {
    if (changece) {
      ++cenum;
      t=0;
      ++base;
    }
    b = evenSamples->getElement(i);
    a = evenSamples->getElement(i-1);
    //assert(b<1);
    if (b==0) {
      b=1.;
      changece=true;
    }
    else 
      changece=false;

    xl = _ce[cenum]._x.f(b) - _ce[cenum]._x.f(a); yl = _ce[cenum]._y.f(b) - _ce[cenum]._y.f(a);    
    segLen =sqrt(xl*xl + yl*yl);
    totSegLen += segLen;
    if (s-st <= segLen) {
      t += ((s-st)/segLen) * (b-a);
      addSample(t,_ce[cenum],base);
      while (segLen >s) {
	if (t>=b-.001) {
	  segLen = -1; st=0;
	}
	else {
	  //printf("starting %f\n",t);  
	  assert(t<1);
	  t += (s/segLen) * (b-a);
	  addSample(t,_ce[cenum],base);
	  a += (s/segLen) * (b-a);
	  segLen -= s;
	}
      } // endwhile
      st = segLen - (s-st);
    }
    else 
      st += segLen;
    t = b;
  } // end for
  if (st > .75*s)
    addSample(1,_ce[cenum],base);
    //else
    //printf("didn't add last sample\n");
}



#define TOL .001

bool BezSpline::flat(const float p, const float q, const float r, const CubicEval2& ce) const {
  Vec2f rLoc(ce._x.f(r), ce._y.f(r));
  float xp = ce._x.f(p) - rLoc.x(); float yp = ce._y.f(p) - rLoc.y();
  float xq = ce._x.f(q) - rLoc.x(); float yq = ce._y.f(q) - rLoc.y();
  // collinearity test:
  // 1. cross product should be small (low triangle area between p,q,r)
  // 2. Dot product should be negative (for strange cases)
  // 3. Added 10/6/03: length should be less than 1
  return ((fabs(xp*yq-xq*yp) < TOL) && (xp*xq + yp*yq <= 0) && ( (xp-xq)*(xp-xq) + (yp-yq)*(yp-yq) < 1.  ));
}

float BezSpline::sample(const float a, const float b, const CubicEval2& ce, 
	     DynArray<float,200>* evenSamples) const {
  if (fabs(a-b) < .001) return 0;
  float t = 0.5; //0.45 + 0.1 * (rand()/(float) RAND_MAX);
  t = a + t*(b-a);
  if (flat (a,b,t,ce)) {
    evenSamples->add(a);
    float xl = ce._x.f(b) - ce._x.f(a), yl = ce._y.f(b) - ce._y.f(a);
    return sqrt(xl*xl + yl*yl);
  }
  else
    return (sample(a,t,ce,evenSamples) + sample(t,b,ce,evenSamples));
}

void BezSpline::clearSamples() {
  _samples.clear();
}


Vec2f BezSpline::getLoc(float t) const {
  int seg = (int)floor(t);
  float rem = t - float(seg);
  Vec2f res;
  if (seg==_numSegs) {
    return _ctrlpts[_ctrlpts.size()-1];
  }
  assert(seg<_numSegs);
  
  res.set_x(_ce[seg]._x.f(rem));
  res.set_y(_ce[seg]._y.f(rem));
  return res;
}

Vec2f BezSpline::getTangent(float t) const {
  int seg = (int)floor(t);
  float rem = t - float(seg);
  Vec2f res;
  if (seg==_numSegs) {
    return _ctrlpts[_ctrlpts.size()-1];
  }
  assert(seg<_numSegs);

  res.set_x(_ce[seg]._x.deriv_1(rem));
  res.set_y(_ce[seg]._y.deriv_1(rem));
  return res;
}

void BezSpline::recalcOneSample(DSample* ds) {
  int seg = (int)ds->_baset;
  float rem = ds->_t;
  assert(seg<_numSegs);
  
  ds->_loc.set_x(_ce[seg]._x.f(rem));
  ds->_loc.set_y(_ce[seg]._y.f(rem));
  ds->_tangent.set_x(_ce[seg]._x.deriv_1(rem));
  ds->_tangent.set_y(_ce[seg]._y.deriv_1(rem));
  //float speed2 = ds->_tangent.x()*ds->_tangent.x() + ds->_tangent.y()*ds->_tangent.y();
  ds->_curvature.set_x(_ce[seg]._x.deriv_2(rem));
  ds->_curvature.set_y(_ce[seg]._y.deriv_2(rem));
}

// t is baset + t
void BezSpline::calcOneSample(DSample* ds, float t) const { 
  int seg = (int)floor(t);
  float rem = t - float(seg);
  Vec2f res;
  if (seg==_numSegs) {
    ds->_loc = _ctrlpts[_ctrlpts.size()-1];
    ds->_tangent.Set(_ce[seg-1]._x.deriv_1(1.), _ce[seg-1]._y.deriv_1(1.));
    //float speed2 = ds->_tangent.x()*ds->_tangent.x() + ds->_tangent.y()*ds->_tangent.y();
    ds->_curvature.Set(_ce[seg-1]._x.deriv_2(1.), _ce[seg-1]._y.deriv_2(1.));
    ds->_t = 1.;
    ds->_baset = _numSegs-1;
    return;
  }
  assert(seg<_numSegs);
  
  ds->_loc.set_x(_ce[seg]._x.f(rem));
  ds->_loc.set_y(_ce[seg]._y.f(rem));
  ds->_tangent.set_x(_ce[seg]._x.deriv_1(rem));
  ds->_tangent.set_y(_ce[seg]._y.deriv_1(rem));
  //float speed2 = ds->_tangent.x()*ds->_tangent.x() + ds->_tangent.y()*ds->_tangent.y();
  ds->_curvature.set_x(_ce[seg]._x.deriv_2(rem));
  ds->_curvature.set_y(_ce[seg]._y.deriv_2(rem));
  ds->_t = rem;
  ds->_baset = seg;
}

void BezSpline::print() const {
  uint i,j;
  for (i=0, j=0; i<(uint)_numSegs; ++i, j+=3) 
    printf("%.10f %.10f, %.10f %.10f, %.10f %.10f, %.10f %.10f\n",
	   _ctrlpts[j].x(), _ctrlpts[j].y(),
	   _ctrlpts[j+1].x(), _ctrlpts[j+1].y(),
	   _ctrlpts[j+2].x(), _ctrlpts[j+2].y(),
	   _ctrlpts[j+3].x(), _ctrlpts[j+3].y());

  printf("///////////////\n");


  for (i=0; i<_samples.size(); ++i) {
    const DSample& ds = _samples[i];
    printf("%d, %.5f: ",i, ds._t + ds._baset);
    printf("l:%.3f %.3f, ", ds._loc.x(), ds._loc.y());
    if (i>0) printf("dl: %f\n",ds._loc.distanceTo(_samples[i-1]._loc));
    else printf("\n");
    //printf("t:%.3f %.3f, ", ds._tangent.x(), ds._tangent.y());
    //printf("c:%.3f %.3f\n", ds._curvature.x(), ds._curvature.y());
  }
  
  
}


void BezSpline::printZMap() const {
  printf("          ZMap: ");
  for (int i=0; i<_mapLength; i++)
    printf("%d ",_ZMap[i]);
  printf("\n");
  
}


void BezSpline::printVarMap() const {
  printf("          VarMap: ");
  for (int i=0; i<_mapLength; i++)
    printf("%d ",_varMap[i]);
  printf("\n");
}

int BezSpline::finalVariable() const {
  assert(_varMap);
  return MAX(_varMap[_mapLength-1], _varMap[_mapLength-2]);
}

int BezSpline::beginningVariable() const { 
  assert(_varMap);
  if (_varMap[1]-_varMap[0] == 1)
    return _varMap[0];
  else
    return _varMap[1];
}

void BezSpline::translate(const Vec2f& delta) {
  uint i;
  for (i=0; i<_ctrlpts.size(); ++i)
    _ctrlpts[i] += delta;

  for (i=0; i<_samples.size(); ++i)
    _samples[i]._loc += delta;
}


float BezSpline::findClosestT(const Vec2f& loc) {
  assert(_samples.size() > 0);
  std::vector<DSample>::iterator v, last=_samples.end()-1;
  float dist2 = FLT_MAX, rest=0, t;
  Vec2f X;
  int i=0;
  for (v=_samples.begin(); v!=last; ++v, ++i) {
    Vec2f L((v+1)->_loc, v->_loc); // L = V1 - V0
    Vec2f P(loc, v->_loc); // P = N - V0
    t = L.Dot2(P) / L.Len();
    t = MAX(MIN(1,t),0); // clamp to 0,1
    Vec2f_LinInterp(X, v->_loc, (v+1)->_loc, t);
    P = loc;
    P -= X;
    float mydist2 = P.Len2();
    if (mydist2 < dist2) {
      rest = v->totT() + t*((v+1)->totT() - v->totT());
      dist2 = mydist2;
    }
  }

  assert(rest != FLT_MAX);
  return rest;
}

void BezSpline::split(float t) {
  // algorithm from ~ page 288 of Solomon
  if (t==0 || t==_numSegs) return;
  int seg = (int)floor(t);
  float tt = t - float(seg);
  assert(t>=0 && t<_numSegs);
  assert(tt>=0 && tt<=1);
  
  int i = seg*3;
  Vec2f P01,P12,P012,P23,P123,P0123;
  Vec2f_LinInterp(P01,   _ctrlpts[i+0], _ctrlpts[i+1], tt);
  Vec2f_LinInterp(P12,   _ctrlpts[i+1], _ctrlpts[i+2], tt);
  Vec2f_LinInterp(P012,  P01,           P12,           tt);
  Vec2f_LinInterp(P23,   _ctrlpts[i+2], _ctrlpts[i+3], tt);
  Vec2f_LinInterp(P123,  P12,           P23,           tt);
  Vec2f_LinInterp(P0123, P012,          P123,          tt);

  _ctrlpts.insert(_ctrlpts.begin()+i+2, 3, Vec2f()); // insert 3 control pts 
  _ctrlpts[i+1] = P01;
  _ctrlpts[i+2] = P012;
  _ctrlpts[i+3] = P0123;
  _ctrlpts[i+4] = P123;
  _ctrlpts[i+5] = P23;
  
  ++_numSegs;
  _ce.push_back(CubicEval2());
  assert(_ctrlpts.size() == _numSegs*3+1);
}


  /*
  float Q1,Q2,Q3,Q4,Q5;
  Vec2f Q6,Q7,B,dB,ddB,dddB;
  DSample k;

  float dt = 1.f/float(T_STEPS);
  Q1 = 3.f*dt;
  Q2 = Q1*dt;
  Q3 = dt*dt*dt;
  Q4 = 2.f*Q2;
  Q5 = 6.f*Q3;
  Vec2f_AddScale(Q6,points[0],points[1],-2.f);
  Q6 += points[2];
  Vec2f_Sub(Q7,points[1],points[2]);
  Q7 *= 3.f; Q7 -= points[0]; Q7 += points[3];
  B = points[0];
  Vec2f_Sub(dB,points[1],points[0]);
  dB *= Q1;
  Vec2f_AddScale(dB,dB,Q6,Q2); Vec2f_AddScale(dB,dB,Q7,Q3);
  Vec2f_WeightedSum(ddB,Q6,Q4,Q7,Q5);
  Vec2f_CopyScale(dddB,Q7,Q5);
  //if (drawFirstPoint)
  //_samples.add(DSample(B,dB,0));
 
  
  float t=0;
  bool stop = 0;
  bool halved=0;
  Vec2f tmp1,tmp2;
  while (!stop) {

    while (fabs(dB.x()) > 1.5 || fabs(dB.y()) > 1.5) {  // halve our step size
      tmp2 = dddB;
      dddB *= .125f;
      tmp1 = ddB;
      Vec2f_WeightedSum(ddB,tmp1,.25f,tmp2,-.125f);
      dB *= 0.5f;
      tmp1 *= -.125f; dB += tmp1;
      tmp2 *=.0625f; dB += tmp2;
      
      dt /=2.;
      halved = 1;
    }
    if (!halved)  // don't undo any halving
      while (fabs(dB.x()) < .75 && fabs(dB.y()) < .75) {  // double it
	tmp2 = dddB;
	dddB *= 8.f;
	tmp1 = ddB;
	Vec2f_WeightedSum(ddB,tmp1,4,tmp2,4);
	dB *= 2.f;
	dB += tmp1;
	
	dt *= 2.f;
      }
    
    halved = 0;
    t += dt;
    B += dB;
    dB += ddB;
    ddB += dddB;

    if (t>1) {
	stop = 1; 
	break; // don't add last point, for some reason it always shows up
    } 

    _samples.add(DSample(B, t));
    
  } // end over Bezier scan conversion loop
    
    //-------------    
    */


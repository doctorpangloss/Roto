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



#include <float.h>
#include <GL/gl.h>
#include "rotoPath.h"
#include "GraphicsGems.h"


namespace FitCurves {
  extern "C" {
#include "FitCurves.h"
    
    void DrawBezCurve(int n, BezierCurve curve) {
      AbstractPath::_talkFit.addBez(curve[0].x,curve[0].y,curve[1].x,curve[1].y,
				     curve[2].x,curve[2].y,curve[3].x,curve[3].y);
    }
  }
}


//#define _1root5_ 0.44721359549995792770 
//#define _2root5_ 0.89442719099991585541

RotoPath::RotoPath() : AbstractPath() {
  mutualInit();
  /*_trans[0] = _trans[3] = 1.;
  _trans[1] = _trans[2] = 0.;
  _trans[4] = _trans[5] = 0;*/
}

RotoPath::RotoPath(const RotoPath& other) : AbstractPath(other) {
  mutualInit();
  _lh = other._lh;
  _hh = other._hh;
  _trackw = other._trackw;
  _trackEdges = other._trackEdges;
  if (other._bez)
    _bez = new BezSpline(other._bez);
  else
    _bez = NULL;
  _samplingDirty = other._samplingDirty;
  _trimapWidth = other._trimapWidth;
  //if (other._trackSubsampling > 0)
  //setTrackSubampling(other._trackSubsampling);
  //memcpy(_trans, other._trans, sizeof(double)*6);
  //_lastTranslation = other._lastTranslation;
}

void RotoPath::mutualInit() {
  _prevFrame = _nextFrame = NULL;
  //_nextCoors = _prevCoors = NULL;
  _lh = -4; _hh = 4;
  _trackEdges = true;
  //setTrackSubampling(3);
  _trackw = 7;
  //_glMatrix = NULL;
  _fixed = true;
  _endsFixed[0] = _endsFixed[1] = false;
  _bez = NULL;
  _nextCont = _prevCont = NULL;
  _samplingDirty = true;
  _xmlLabel = -1;
  _regions[0] = _regions[1] = NULL;
  _trimapWidth = 6;// G!
}

RotoPath::~RotoPath() {
  //if (_nextCoors) delete[] _nextCoors; 
  //if (_prevCoors) delete[] _prevCoors;
  if (_nextCont) delete _nextCont;
  if (_prevCont) delete _prevCont;
  if (_bez) {
    //printf("Rotopath deleting %x\n",_bez);
    delete _bez;
  }
}

/*void RotoPath::setTrackSubsampling(const int s) {
  assert(s>0);
  _trackSubsampling = s;
  if (_subLocs != NULL)
    delete[] _subLocs;
  _nsubLocs = ((getNumElements()-1) / s) + 1;
  assert(_nsubLocs > 0);
  _subLocs = new Vec2f[_nsubLocs];
  int i, s_i;
  for (i=0, s_i=0; i<getNumElements(); i+=s, s_i++) {
    _subLocs[s_i] = getElement(i);
    assert(s_i < _nsubLocs);
  }
  
  _trackw = MAX(5,2*s+1);
  //printf("trackw set to %d\n",_trackw);
}
*/


void RotoPath::clearPrevCorrespondence() {
  if (_prevFrame) {
    _prevFrame->_nextFrame = NULL;
    //if (_prevFrame->_nextCoors) {
    //delete _prevFrame->_nextCoors;
    //_prevFrame->_nextCoors = NULL;
    //}
    if (_prevFrame->_nextCont) {
      delete _prevFrame->_nextCont;
      _prevFrame->_nextCont = NULL;
    }    
  }
}

/*
void RotoPath::finishPrevCorrespondence() {
  if (!_prevFrame) return;
  assert(_prevFrame && _prevFrame->_nextFrame == this && !_prevFrame->_nextCoors);
  _prevFrame->_nextCoors = new int[_prevFrame->getNumElements()];

  RotoCorresponder* myCoor = new RotoCorresponder(_prevFrame, this);
  myCoor->calculate();
  myCoor->getSolution(_prevFrame->_nextCoors);
  delete myCoor;

  _prevCoors = new int[getNumElements()];
  RotoCorresponder* myCoorR = new RotoCorresponder(this, _prevFrame);
  myCoorR->calculate();
  myCoorR->getSolution(_prevCoors);
  delete myCoorR;
  }*/


void RotoPath::calculateBackCorrs() {
  assert(_prevFrame);
  RotoCorresponder* myCoor = new RotoCorresponder(_prevFrame, this);
  myCoor->calculate();
  int* discreteCorrs = new int[_prevFrame->getNumElements()];
  myCoor->getSolution(discreteCorrs);
  assert(!_prevFrame->_nextCont);
  _prevFrame->_nextCont = new ContCorr(_prevFrame->_bez, _bez, discreteCorrs);
  delete myCoor;
  delete[] discreteCorrs;

  myCoor = new RotoCorresponder(this, _prevFrame);
  myCoor->calculate();
  discreteCorrs = new int[getNumElements()];
  myCoor->getSolution(discreteCorrs);
  _prevCont = new ContCorr(_bez, _prevFrame->_bez, discreteCorrs);
  delete myCoor;
  delete[] discreteCorrs;
}

void RotoPath::startPrevCorrespondence(RotoPath* prevC) {
  prevC->_nextFrame = this;
  _prevFrame = prevC;
}

void RotoPath::finishPrevCorrespondence() {
  if (!_prevFrame) return;
  assert(_prevFrame && _prevFrame->_nextFrame == this && !_prevFrame->_nextCont);

  assert(!_samplingDirty && !_prevFrame->_samplingDirty);
  assert(_prevFrame->getNumElements() == _prevFrame->_bez->getDiscreteCount());
  assert(getNumElements() == _bez->getDiscreteCount());

  calculateBackCorrs();

  /*
  _prevCoors = new int[getNumElements()];
  RotoCorresponder* myCoorR = new RotoCorresponder(this, _prevFrame);
  myCoorR->calculate();
  myCoorR->getSolution(_prevCoors);
  delete myCoorR;*/
}

void RotoPath::handleCorrForw(RotoPath* other) {
  _nextFrame = other;
  other->_prevFrame = this;
  other->calculateBackCorrs();
  /*other->_prevCoors = new int[other->getNumElements()];
  _nextCoors = new int[getNumElements()];

  RotoCorresponder* myCoor = new RotoCorresponder(other, this);
  myCoor->calculate();
  myCoor->getSolution(other->_prevCoors);
  delete myCoor;


  RotoCorresponder* myCoorR = new RotoCorresponder(this, other);
  myCoorR->calculate();
  myCoorR->getSolution(_nextCoors);
  delete myCoorR;*/

  
}


bool RotoPath::assertState() const {
  /*
    int i;
    bool result = 1;
  
  Vec2f lastLoc = getElement(0), newLoc;
  for (i=1; i<getNumElements(); i++) {
    newLoc = getElement(i);
    if (!(fabs(newLoc.x()-lastLoc.x()) <=1 &&
	  fabs(newLoc.y()-lastLoc.y()) <=1))
      result = 0;
    assert(result);
    lastLoc = newLoc;
  }
  return result;
  */
  return true;
}

// SPEED: this is slow


// CORRECTNESS!: should use _tangent array!!
void RotoPath::getUV(const float t, Vec2f& U, Vec2f& V) const { 

  assert(_bez);
  U = getTangent(t);
  U.Normalize();
  assert(finite(U.x()) && finite(U.y()));
  V.Set(-U.y(), U.x()); // perpendicular to U

  /*
  //  ..assert(i>=0 && i<getNumElements())
  if (i==0)
    Vec2f_Sub(U,getElement(1), getElement(0));
  else if (i==getNumElements()-1) 
    Vec2f_Sub(U,getElement(getNumElements()-1), getElement(getNumElements()-2));
  else {
    Vec2f_Sub(U,getElement(i+1), getElement(i-1));
  }
  U.Normalize();
  V.Set(-U.y(), U.x()); // perpendicular to U
  */
}


/*Vec2f RotoPath::getNextCorrPoint(const int i) const {
  assert (_nextFrame && i>=0 && i<getNumElements() && _nextCont && _bez);
  return (_nextFrame->getElement(_nextCoors[i]));
  
}

int RotoPath::getNextCorrIndex(const int i) const {
  assert(_nextCoors && i>=0 && i<getNumElements());
  return _nextCoors[i];
  
}

int RotoPath::getPrevCorrIndex(const int i) const {
  assert(_prevCoors && i>=0 && i<getNumElements());
  return _prevCoors[i];
  }
*/




void RotoPath::render(bool showTrackPoints, bool visEffort) const {
  int i;

 // render debug grid

  glLineWidth(3.0); //G! 4 // 8 for pictures
  //if (_trimapWidth != -1)
  //glLineWidth(_trimapWidth); // G!
  // render bezier evaluators
  if (_bez && (getNumElements() == 0) ) { // G! || 1
    assert(!glGetError());
    glEnable(GL_MAP1_VERTEX_3);
    int n = _bez->numCtrls(), i3;
    n = ((n-1)/3) * 3 + 1;
    int numSegs = (n-1)/3;
    //printf("manual render segs %d ctrls %d\n", numSegs, n);
    if (numSegs > 0) {
      GLfloat* ctrls = new GLfloat[3*n];
      for (i3=0, i=0; i<n; ++i, i3+=3) {
	const Vec2f* v = _bez->getCtrl(i);
	ctrls[i3] = v->x(); ctrls[i3+1] = v->y(); ctrls[i3+2] = 0;
	//printf("%f %f %f\n", ctrls[i3], ctrls[i3+1], ctrls[i3+2]);
      }
      
      GLfloat* cptr = ctrls;
      for (i=0; i<numSegs; ++i, cptr+=9) {
	glMap1f(GL_MAP1_VERTEX_3, 0., 1., 3, 4, cptr);
	glMapGrid1f(25,0.,1.);
	//glBegin(GL_LINES);
	//printf("doing j's\n");
	glEvalMesh1(GL_LINE,0,25);
	//for (int j=0; j<25; ++j)
	//glEvalPoint1(j);
	//glEnd();
      }
      assert(!glGetError());
      delete[] ctrls;
    }
  }


  // render tracking points
  GLboolean stipple;
  glGetBooleanv(GL_LINE_STIPPLE,&stipple);
  if (!stipple &&  showTrackPoints) {
    glPushAttrib(GL_CURRENT_BIT);
    glColor3f(0,0,1);
    Vec2f loc;
    glBegin(GL_POINTS);
    /*
    int hw = _trackw / 2;
    float finvc = 1./(2.*_trackSubsampling);
    for (int n=0; n<_nsubLocs-1; n++) {
      for (int j = _lh ; j <= _hh ; j++)  
	for (int i = 0 ; i <= hw ; i++) {
	  hnew(n,i,j, finvc, _subLocs, loc);
	  glVertex3f(loc.x(), loc.y(), 0);
	}
    }
    */

    int hw;
    int ul = getNumElements()-1;
    for (int n=0; n<ul; n++) {
      if (n==ul-1)
	hw=2;
      else
	hw=1;
      for (int j = _lh ; j <= _hh ; j++)  
	for (int i = 0 ; i <= hw ; i++) {
	  hnew(n,i,j, getData(), loc);
	  glVertex3f(loc.x(), loc.y(), 0);
	}
    }    
    
    glEnd();
    glPopAttrib();
  }
  

  glLineWidth(3.0); //G!
  glBegin(GL_LINE_STRIP); // G! (lines or points)
  for (i=0; i<getNumElements(); i++) {
    const Vec2f *dp = getPointerToElement(i);
    glVertex2f(dp->x(),dp->y());
  }
  glEnd();
  

  // render bezier control points
  glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT);
  if (_bez) _bez->renderCtrlPts(); 
  glPopAttrib();


  
  
  if (_bez && !_bez->empty() && (!_bJoints.empty() || !_eJoints.empty())) { // G!  (0 && )
    glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT);
    glColor3f(0,1,1); 
    glPointSize(5.); 
    glBegin(GL_POINTS);
    if (!_bJoints.empty()) {
      const Vec2f& dp = getEnd(0);            
      glVertex2f(dp.x(),dp.y());
    }
    if (!_eJoints.empty()) {
      const Vec2f& dp = getEnd(1);
      glVertex2f(dp.x(),dp.y());
    }
    glEnd();
    glPopAttrib();
  }

  //if (_bez && (anyEndFixed() || !_internalsFixed.empty())   && !_fixed) {
  if (visEffort && _bez && !_bez->empty() && int(_touched.size()) == _bez->numCtrls()) {
    std::set<int>::const_iterator c;
    glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT);
    glColor3f(.627,.125,.941);
    glPointSize(7.); 
    glBegin(GL_POINTS);
    // internal fixed
    /*
    for (c=_internalsFixed.begin(); c!= _internalsFixed.end(); ++c) {  
      const Vec2f* v = _bez->getCtrl(*c);
      glVertex2f(v->x(),v->y());
    }
    // endpoint fixed
    for (int i=0; i<2; ++i) {
      if (_endsFixed[i]) {
	//printf("end fixed render\n");
	const Vec2f& v = _bez->getEnd(i);
	glVertex2f(v.x(),v.y());
      }
    }
    */
    for (int i=0; i<_touched.size(); ++i) {
      if (touched(i)) {
	const Vec2f* v = _bez->getCtrl(i);
	glVertex2f(v->x(),v->y());
      }
    }
    glEnd();
    glPopAttrib();
  }


  // render crosshairs
  /*
  glPushAttrib(GL_CURRENT_BIT);
  glColor3f(0,0,1);
  glTranslatef(myCenter.x(), myCenter.y(),0);
  glBegin(GL_LINES);
  glVertex3i(-8,0,0);
  glVertex3i(8,0,0);
  glVertex3i(0,-8,0);
  glVertex3i(0,8,0);
  glEnd();
  glPopAttrib();
  */

  
  

  //glPopMatrix();
}

/*
void RotoPath::renderEmphasizedEndpoint(int which) const {
   glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT);
    glColor3f(1,0,0);
    glPointSize(3.);
    glBegin(GL_POINTS);
    const Vec2f* dp;
    if (which==0)
      dp = getPointerToElement(0);
    else
      dp = getPointerToElement(getNumElements()-1);
    glVertex2f(dp->x(),dp->y());
    glEnd();
    glPopAttrib();
}
*/

void RotoPath::hnew(const int n, const float i, const float j, 
       const Vec2f* P, Vec2f& result) const { 
  Vec2f Tn(P[n+1],P[n]);
  Tn *= .5f; //invc;
  result = P[n];
  result.Inc(i * Tn.x(),  i * Tn.y());
  result.Inc(j * -Tn.y(), j * Tn.x());
}


void RotoPath::h(const int n, const float i, const float j, const int l, const float invc,
       const Vec2f* P, Vec2f& result) const { 
  Vec2f Tn, Nn;
  assert(n<l);
  if (n==l-1) {
    Vec2f_Sub(Tn, P[l-1],P[l-2]);
    Tn *= invc;
  }
  else if (n==0) {
    Vec2f_Sub(Tn, P[1],P[0]);
    Tn *= invc;
  }
  else {
    Vec2f_Sub(Tn, P[n+1],P[n-1]);
    Tn *= invc;
    Tn /= 2.f;
  }
  Nn.Set(-Tn.y(), Tn.x());

  result = P[n];
  result.Inc(i*Tn.x(), i*Tn.y());
  result.Inc(j*Nn.x(), j*Nn.y());
}

void RotoPath::hnorm(const int n, const float i, const float j, const int l,
       const Vec2f* P, Vec2f& result) const { 
  Vec2f Tn, Nn;
  if (n==l-1) {
    Vec2f_Sub(Tn, P[l-1],P[l-2]);
    Tn.Normalize();
  }
  else if (n==0) {
    Vec2f_Sub(Tn, P[1],P[0]);
    Tn.Normalize();
  }
  else {
    Vec2f_Sub(Tn, P[n+1],P[n-1]);
    Tn.Normalize();
  }
  Nn.Set(-Tn.y(), Tn.x());

  result = P[n];
  result.Inc(i*Tn.x(), i*Tn.y());
  result.Inc(j*Nn.x(), j*Nn.y());
}

float RotoPath::getSampleNextT(int s) const {
  assert(_nextCont && !_samplingDirty && _bez);
  if (_nextCont->identity())
    return _bez->getDiscreteT(s);
  else
    return _nextCont->getT(_bez->getDiscreteT(s));
}

float RotoPath::getSampleT(int s) const {
  assert(_bez);
  return _bez->getDiscreteT(s);
}

float RotoPath::getNextT(float s) const {
  assert(_nextCont && !_samplingDirty && _bez && finite(s));
  if (_nextCont->identity())
    return s;
  else
    return _nextCont->getT(s);
}

float RotoPath::getPrevT(float s) const {
  assert(_prevCont && !_samplingDirty && _bez && finite(s));
  if (_prevCont->identity())
    return s;
  else
    return _prevCont->getT(s);
}


const Vec2f RotoPath::getNextLoc(int s) const { 
  assert(_nextFrame); 
  //printf("%f %f\n",_bez->getDiscreteT(s), getSampleNextT(s));
  return _nextFrame->getLoc(getSampleNextT(s)); 
}

void RotoPath::renderNextCorrLines() const {
  if (!(_nextFrame && _nextCont && !_samplingDirty)) return;
  //const Vec2f aloc, bloc;
  glBegin(GL_LINES);
  for (int i=0; i<getNumElements(); i+=1) { 
    const Vec2f& aloc = getElement(i);
    //bloc = getNextCorrPoint(i);
    const Vec2f& bloc = getNextLoc(i);
    glVertex2f(aloc.x(), aloc.y());
    glVertex2f(bloc.x(), bloc.y());
  }
  glEnd();
}

void RotoPath::renderSelected() const {
  //return; // ONE
  /*Vec2f myCenter;
  getCentroid(myCenter);
  glPushMatrix();

  glTranslatef(myCenter.x(), myCenter.y(),0);

  if (!_glMatrix)
    buildGLTransMatrix();
  glMultMatrixd(_glMatrix);  

  glTranslatef(-myCenter.x(), - myCenter.y(),0);
  */
  /*
  Vec2f loc = getElement(0);
  glRectf(loc.x()-2, loc.y()-2, 
	  loc.x()+2, loc.y()+2);
  loc = getElement(getNumElements()-1);
  glRectf(loc.x()-2, loc.y()-2, 
	  loc.x()+2, loc.y()+2);
  */
  //glPopMatrix();
  glColor3f(153./256., 51./256., 153./256.);
  render(false);
  //glLineWidth(1.); //G!
  //renderHandles();
  
}

void RotoPath::renderTransformSelected() {
  //printf("renderTransformSelected start\n");  
  glColor3f(153./256., 51./256., 153./256.);

  /*
  // render bezier evaluators
  if (_bez && getNumElements() == 0) {
    assert(!glGetError());
    glEnable(GL_MAP1_VERTEX_3);
    int n = _bez->numCtrls(), i3, i;
    n = ((n-1)/3) * 3 + 1;
    int numSegs = (n-1)/3;
    //printf("manual render segs %d ctrls %d\n", numSegs, n);
    if (numSegs > 0) {
      GLfloat* ctrls = new GLfloat[3*n];
      for (i3=0, i=0; i<n; ++i, i3+=3) {
	const Vec2f* v = _bez->getCtrl(i);
	ctrls[i3] = v->x(); ctrls[i3+1] = v->y(); ctrls[i3+2] = 0;
	//printf("%f %f %f\n", ctrls[i3], ctrls[i3+1], ctrls[i3+2]);
      }
      
      GLfloat* cptr;
      //for (i=0; i<numSegs; ++i, cptr+=9) {
      for (std::vector<int>::const_iterator c=_transformSelected.begin();
	   c!= _transformSelected.end(); ++c) {
	printf("renderTransformSelected %d\n",*c);
	cptr = ctrls + (*c) * 9;
	glMap1f(GL_MAP1_VERTEX_3, 0., 1., 3, 4, cptr);
	glMapGrid1f(25,0.,1.);
	//glBegin(GL_LINES);
	//printf("doing j's\n");
	glEvalMesh1(GL_LINE,0,25);
	//for (int j=0; j<25; ++j)
	//glEvalPoint1(j);
	//glEnd();
      }
      assert(!glGetError());
      delete[] ctrls;
    }
  }
  */

  /*
  int i,j;
  //for (i=0; i<_transformSelected.size(); ++i)
  //printf("transformSelected %d\n",i);
  assert(_bez && _bez->getDiscreteCount() > 0);
  glBegin(GL_POINTS);
  int lastbaset = -1;
  bool canRender = false;
  //printf("count %d:",_bez->getDiscreteCount());
  for (i=0; i<_bez->getDiscreteCount(); ++i) {
    j = (int) _bez->getDiscreteSample(i)._baset;
    if (j!=lastbaset) {
      lastbaset = j;
      if (find(_transformSelected.begin(), _transformSelected.end(), j) == _transformSelected.end())
	canRender = false;
      else
	canRender = true;
    }
    
    if (canRender) {
      //  printf("%f ",_bez->getDiscreteSample(i).totT());
      const Vec2f dp = _bez->getDiscreteLoc(i);
      glVertex2f(dp.x(),dp.y());
    }
  }
  glEnd();
  //printf("\n");
  */

  glBegin(GL_POINTS);
  SplineSampleIterator c2 = transformSelectedSampleIterator();
  for (; !c2.end(); ++c2) {
    const Vec2f& v = c2.sample()._loc;
    glVertex2f(v.x(), v.y());
  }
  glEnd();
}

void RotoPath::transformSelect(int t) {
  assert(_bez);
  assert(t<_bez->numSegs() && t>=0);
  _transformSelected.push_back(t);
  _transformSelected.sort();
}

void RotoPath::translateSelected(const Vec2f& delta) {
  assert(_bez);
  SplineCtrlIterator c = transformSelectedCtrlIterator();
  for (; !c.end(); ++c)
    c.ctrl() += delta;

  SplineSampleIterator c2 = transformSelectedSampleIterator();
  for (; !c2.end(); ++c2)
    c2.sample()._loc += delta;

}

void RotoPath::rotateSelectedAbout(float dtheta, const Vec2f center) {
  assert(_bez);
  float costheta =  cos(dtheta), sintheta = sin(dtheta);
  Vec2f a(costheta, -sintheta),
        b(sintheta, costheta);

  SplineCtrlIterator c = transformSelectedCtrlIterator();
  for (; !c.end(); ++c) {
    Vec2f p(c.ctrl(), center);
    Vec2f rotp(p.Dot2(a), p.Dot2(b));
    rotp += center;
    c.ctrl() = rotp;
  }

  SplineSampleIterator c2 = transformSelectedSampleIterator();
  for (; !c2.end(); ++c2) {
    Vec2f p(c2.sample()._loc, center);
    Vec2f rotp(p.Dot2(a), p.Dot2(b));
    rotp += center;
    c2.sample()._loc = rotp;
  }
}

void RotoPath::uniformScaleSelectedAbout(float mag, const Vec2f& center) {
  //printf("mag by %f\n",mag);
  
  //uint i;
  //for (i=0; i<_ctrlpts.size(); ++i) {
  SplineCtrlIterator c = transformSelectedCtrlIterator();
  for (; !c.end(); ++c) {
    Vec2f p(c.ctrl(), center);
    p *= mag;
    p += center;
    c.ctrl() = p;
  }
  
  //for (i=0; i<_samples.size(); ++i) {
  SplineSampleIterator c2 = transformSelectedSampleIterator();
  for (; !c2.end(); ++c2) {
    Vec2f p(c2.sample()._loc, center);
    p *= mag;
    p += center;
    c2.sample()._loc = p;
  }
  
  
}


void RotoPath::doublePoints() {
  int newcount = getNumElements()*2-1, i;
  Vec2f *newData = new Vec2f[newcount];
  for (i=0; i<getNumElements(); i++) // copy over old points
    newData[2*i] = getElement(i);
  for (i=1; i<newcount; i+=2) {
    Vec2f_Average(newData[i],newData[i-1],newData[i+1]);
  }
  delete[] data;
  data = newData;
  space = newcount;
  count = newcount;
}

void RotoPath::reconcileFixEndpointJoint(const int side) {
  Joints& jts = joints(side);
  Joints::iterator c;
  for (c=jts.begin(); c!= jts.end(); ++c) {
    c->_rp->_endsFixed[c->_side] = true;
  }
}

void RotoPath::fixEndpoint(int which) {
  assert(which==0 || which==1);
  _endsFixed[which] = true;
  reconcileFixEndpointJoint(which);
}

void RotoPath::flipFixed() { 
  _fixed = !_fixed; 
  if (!_fixed)
    setTouchedAll(false);
}

void RotoPath::fixInternal(int which) {
  assert(_bez && which > 0 && which < _bez->numCtrls()-1);
  _internalsFixed.insert(which);
}

bool RotoPath::internalFixed(const int which) const {
  return (find(_internalsFixed.begin(), _internalsFixed.end(), which) != _internalsFixed.end());
}

bool RotoPath::somehowFixed(const int which) const {
  assert(which>=0 && which < _bez->numCtrls());
  if (which==0)
    return endFixed(0);
  if (which==_bez->numCtrls()-1)
    return endFixed(1);
  return internalFixed(which);
}

void RotoPath::unfixCtrl(const int which) {
  assert(_touched.size() == _bez->numCtrls());
  assert(which>=0 && which < _bez->numCtrls());
  if (which==0) {
    _endsFixed[0] = 0;
    _touched[0] = false;
  }
  else if (which==_bez->numCtrls()-1) {
    _endsFixed[1] = 0;
    _touched[_bez->numCtrls()-1] = false;
  }
  else {
    _internalsFixed.erase(which);
    _touched[which] = false;
  }
}


/*
void RotoPath::renderNextCorrPS(FILE* fp) const {
  Vec2f a,b;
  if (!_nextFrame) return;
  for (int i=0; i<getNumElements(); i++) {
    a = getElement(i);
    b = _nextFrame->getElement(_nextCoors[i]);
    fprintf(fp,"%f %f moveto\n",a.x(),_globalh - a.y());
    fprintf(fp,"%f %f lineto\n",b.x(),_globalh - b.y());
    fprintf(fp,"closepath stroke\n");
  }
  }*/

void RotoPath::dumpPoints(FILE* fp) const {
  int dummy = getNumElements();
  fwrite(&dummy,sizeof(int),1,fp); // write # of points
  fwrite(getData(),sizeof(Vec2f),dummy,fp); // write points
}

void RotoPath::saveXML(std::ofstream& fp) const {
  assert(_bez && _xmlLabel != -1);
  fp << "<Track: {label" << _xmlLabel << "} ";
  int numSegs = _bez->numSegs();
  const Vec2f* v = _bez->getCtrl(0);
  fp << "M " << v->x() << " " << v->y() << " ";
  for (int i=0, index=1; i<numSegs; ++i) {
    fp << "B ";
    for (int j=0; j<3; ++j, ++index) {
      v = _bez->getCtrl(index);
      fp << v->x() << " " << v->y() << " ";
    }
  }
  fp << ">" << std::endl;
}

void RotoPath::labelSelfXML(int& labelNum) {
  if (_xmlLabel == -1) {
    _xmlLabel = labelNum;
    ++labelNum;
  }
  else
    return;

  // label all later corresponding curves with same label
  RotoPath* curr = this->nextC();
  while (curr) {
    assert(curr->_xmlLabel == -1);
    curr->_xmlLabel = _xmlLabel;
    curr = curr->nextC();
  }
}

void RotoPath::saveqt(QDataStream * fp, int version) const {
  int dummy,i;

  dummy = getNumElements();
  *fp << dummy;               // write # of points
  for (i=0; i<dummy; ++i)
    *fp << getData()[i].x() << getData()[i].y();     // write points

  if (_tangents) {
    *fp << dummy;               // write # of points
    for (i=0; i<dummy; ++i)
      *fp << _tangents[i].x() << _tangents[i].y();     // write points
  }
  else {
    dummy=0; *fp << dummy;
  }

  dummy = getNumCorrDraws();
  *fp << dummy;
  for (set<DrawPath*>::const_iterator c = _corrDraws.begin(); c!=_corrDraws.end(); ++c)
    (*c)->writePtr(fp);

  if (_prevFrame)    // rotopath _prevFrame pointer
    _prevFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  if (_nextFrame)    // _nextFrame pointer 
    _nextFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  dummy=0;
  //if (_prevCoors)    // rotopath _prevCoors intArray
  //writeIntArray(fp,_prevCoors, getNumElements());
  //else
  //*fp << dummy;
  //if (_nextCoors)    // _nextCoors intArray
  //writeIntArray(fp,_nextCoors, getNumElements());
  //else
  //*fp << dummy;

  /*
  // for passing to old code 

  if (_prevCont)    //rotopath _prevCoors intArray
    writeIntArray(fp,_prevCont->getDiscreteCorrs(), _prevCont->getNumDiscreteCorrs());
  else
    *fp << dummy;
  if (_nextCont)    // _nextCoors intArray
    writeIntArray(fp,_nextCont->getDiscreteCorrs(), _nextCont->getNumDiscreteCorrs());
  else
    *fp << dummy;
    */
  //------------------

  *fp << _lh << _hh;    // _lh, _hh, _trackw
  if (version==0) {
    *fp << _trackw;
    dummy=0; *fp << dummy;
  }

  //*fp << (Q_INT8)_fixed;
  writeBool(fp, _fixed);


  if (version >=1) {
    //*fp << (Q_INT8)_trackEdges;
    writeBool(fp, _trackEdges);
    // joints
    dummy = _bJoints.size();
    *fp << dummy;
    for (i=0; i<dummy; ++i) {
      _bJoints[i]._rp->writePtr(fp);
      *fp << _bJoints[i]._side;
    }

    dummy = _eJoints.size();
    *fp << dummy;
    for (i=0; i<dummy; ++i) {
      _eJoints[i]._rp->writePtr(fp);
      *fp << _eJoints[i]._side;
    }
  }

  if (version>=2) {
    //*fp << (Q_INT8)_endsFixed[0] << (Q_INT8)_endsFixed[1];
    writeBool(fp, _endsFixed[0]);
    writeBool(fp, _endsFixed[1]);
  }
  
  
  if (version >=3) {  // write spline-rep

    if (!_bez) 
      dummy = -1;
    else
      dummy = _bez->numSegs();
    *fp << dummy;
    for (i=0; i<_bez->numCtrls(); ++i) {
      const Vec2f& v = *(_bez->getCtrl(i));
      *fp << v.x() << v.y();      
    }

    writeBool(fp, _samplingDirty);
  }

  if (version >=4) {
    if (_prevCont) {
      *fp << 1;  _prevCont->saveqt(fp);
    }
    else
	*fp << 0;
    
    if (_nextCont) {
      *fp << 1;  _nextCont->saveqt(fp);
    }
    else
      *fp << 0;
  }

  if (version >= 5) {
    *fp << _internalsFixed.size();
    for (std::set<int>::const_iterator c = _internalsFixed.begin(); 
	 c != _internalsFixed.end(); ++c)
      *fp << *c;
  }
    
  if (version >= 6) {
    *fp << _touched.size();
    for (std::vector<bool>::const_iterator c = _touched.begin();
	 c != _touched.end(); ++c)
      writeBool(fp, *c);
  }

  if (version >= 7) {
    *fp << _trimapWidth;
  }
  
}

void RotoPath::save(FILE* fp, int version) const {
  int dummy;

  dummy = getNumElements();
  fwrite(&dummy,sizeof(int),1,fp); // write # of points
  fwrite(getData(),sizeof(Vec2f),dummy,fp); // wrbite points
  
  if (_tangents) {
    fwrite(&dummy,sizeof(int),1,fp); // write # points
    fwrite(_tangents,sizeof(Vec2f),dummy,fp);
  }
  else {
    dummy = 0;
    fwrite(&dummy,sizeof(int),1,fp); // write null pointer
  }

  dummy = getNumCorrDraws();
  fwrite(&dummy,sizeof(int),1,fp); // write # of corresponding drawPaths    
  for (set<DrawPath*>::const_iterator c = _corrDraws.begin(); c!=_corrDraws.end(); ++c)
    (*c)->writePtr(fp);
  

  if (_prevFrame)    // rotopath _prevFrame pointer
    _prevFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  if (_nextFrame)    // _nextFrame pointer 
    _nextFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  dummy=0;
  //if (_prevCoors)    // rotopath _prevCoors intArray
  //writeIntArray(fp,_prevCoors, getNumElements());
  //else
  //fwrite(&dummy,sizeof(int),1,fp);
  //if (_nextCoors)    // _nextCoors intArray
  //writeIntArray(fp,_nextCoors, getNumElements());
  //else
  //fwrite(&dummy,sizeof(int),1,fp);

  fwrite(&_lh, sizeof(int),1,fp);
  fwrite(&_hh, sizeof(int),1,fp);
  if (version==0) {
    fwrite(&_trackw, sizeof(int),1,fp); 
    dummy=0;
    fwrite(&dummy, sizeof(int),1,fp);  
  }
  fwrite(&_fixed, sizeof(bool),1,fp);  // is this version 0 or 1? I don't know


  if (version>=1) {

    int i;

    fwrite(&_trackEdges, sizeof(bool),1,fp); // wether to track edges
    
    // joints
    dummy = _bJoints.size();
    fwrite(&dummy, sizeof(int),1,fp); // # beginning joints
    for (i=0; i<dummy; ++i) {
      _bJoints[i]._rp->writePtr(fp);
      fwrite(&(_bJoints[i]._side), sizeof(int), 1, fp);
    }

    dummy = _eJoints.size();
    fwrite(&dummy, sizeof(int),1,fp); // # end joints
    for (i=0; i<dummy; ++i) {
      _eJoints[i]._rp->writePtr(fp);
      fwrite(&(_eJoints[i]._side), sizeof(int), 1, fp);
    }
  }

  if (version>=2) {
    fwrite(_endsFixed, sizeof(bool),2,fp); // _endsFixed
  }

}



void RotoPath::loadqt(QDataStream* fp, int version) {
  int dummy,i, dptr[2];
  
  *fp >> dummy;
  addNElements(dummy);
  for (i=0; i<dummy; ++i)
    *fp >> getData()[i].data()[0] >> getData()[i].data()[1];  

  *fp >> dummy;
  if (!dummy)
    _tangents = NULL;
  else {
    assert(dummy==getNumElements());
    _tangents = new Vec2f[getNumElements()]; 
    for (i=0; i<dummy; ++i)
      *fp >> _tangents[i].data()[0] >> _tangents[i].data()[1];     
  }

  *fp >> dummy;
  if (dummy > 0) {
    for (i=0; i<dummy; i++) {
      *fp >> dptr[0] >> dptr[1];
      _corrDraws.insert(_is->resolveDrawWritePtr(dptr));
    }
  }

  *fp >> dptr[0] >> dptr[1]; // _prevFrame
  if (dptr[0] == -1)
    _prevFrame = NULL;
  else
    _prevFrame = _is->resolveRotoWritePtr(dptr);
  
  *fp >> dptr[0] >> dptr[1]; // _nextFrame
  if (dptr[0] == -1)
    _nextFrame = NULL;
  else
    _nextFrame = _is->resolveRotoWritePtr(dptr);

  /**fp >> dummy; // _prevCoors
  if (dummy==0) {
    _prevCoors = NULL;
    assert(_prevFrame==NULL);
  }
  else {
    _prevCoors = new int[dummy];
    readIntArray(fp,_prevCoors,dummy);
  }

  *fp >> dummy; // _nextCoors
  if (dummy==0) {
    _nextCoors = NULL;
    assert(_nextFrame==NULL);
  }
  else {
    _nextCoors = new int[dummy];
    readIntArray(fp,_nextCoors,dummy);
    }*/
    
  *fp >> _lh >> _hh;
  if (version==0) {
    *fp >> _trackw;
    *fp >> dummy;
  }
  //*fp >> (Q_INT8&)_fixed;
  readBool(fp, _fixed);

  if (version>=1) {
    //*fp >> (Q_INT8&)_trackEdges;
    readBool(fp, _trackEdges);
    // joints

    *fp >> dummy;
    _bJoints.reserve(dummy);
    for (i=0; i<dummy; ++i) {
      Joint j;
      *fp >> dptr[0] >> dptr[1]; 
      j._rp = _is->resolveRotoWritePtr(dptr);
      *fp >> j._side;
      _bJoints.push_back(j);
    }
    assert(jointsOk(_bJoints));
    
    *fp >> dummy;
    _eJoints.reserve(dummy);
    for (i=0; i<dummy; ++i) {
      Joint j;
      *fp >> dptr[0] >> dptr[1]; 
      j._rp = _is->resolveRotoWritePtr(dptr);
      *fp >> j._side;
      _eJoints.push_back(j);
    }
    assert(jointsOk(_eJoints));

  }

  if (version>=2) {
    //*fp >> (Q_INT8&)_endsFixed[0] >> (Q_INT8&)_endsFixed[1];
    readBool(fp, _endsFixed[0]);
    readBool(fp, _endsFixed[1]);
  }

  if (version >=3) {  // load spline-rep
    *fp >> dummy;
    if (dummy==-1) _bez = NULL;
    else {
      _bez = new BezSpline();
      _bez->setNumSegs(dummy);
      dummy = dummy*3 + 1;
      for (i=0; i<dummy; ++i) {
	Vec2f v;
	*fp >> v.data()[0] >> v.data()[1];
	_bez->addCtrl(v);
      }
      handleNewBezCtrls();
    }

    readBool(fp, _samplingDirty);
  }

  if (version >= 4) {
    *fp >> dummy;
    if (dummy == 1) _prevCont = new ContCorr(fp);
    else _prevCont = NULL;
    
    *fp >> dummy;
    if (dummy == 1) _nextCont = new ContCorr(fp);
    else _nextCont = NULL;  
  }

  if (version >= 5) {
    *fp >> dummy;
    int d2;
    for (i=0; i<dummy; ++i) {
      *fp >> d2;
      _internalsFixed.insert(d2);
    }
  }

  if (version >= 6) {
    *fp >> dummy;
    bool d3;
    for (i=0; i<dummy; ++i) {
      readBool(fp, d3);
      _touched.push_back(d3);
    }
  }
  else
    estimateTouched();

  if (version >= 7) {
    *fp >> _trimapWidth;
  }
}

void RotoPath::load(FILE* fp, int version)  {
  int dummy,i, dptr[2], res;
  
  res = fread(&dummy,sizeof(int),1,fp); // read # of points
  i = ferror(fp);
  assert(i==0);
  i=feof(fp);
  assert(i==0);
  assert(res == 1);
  addNElements(dummy);
  res = fread(getData(),sizeof(Vec2f),dummy,fp); // read points
  assert(res==dummy);

  res = fread(&dummy,sizeof(int),1,fp); // read _tangent status
  assert(res==1);
  if (!dummy)
    _tangents = NULL;
  else {
    assert(dummy = getNumElements());
    _tangents = new Vec2f[getNumElements()]; 
    res = fread(_tangents,sizeof(Vec2f),dummy,fp);
    assert(res==dummy);
  }

  res = fread(&dummy,sizeof(int),1,fp); // read # of corresponding drawPaths
  assert(res == 1);
  if (dummy > 0) {
    for (i=0; i<dummy; i++) {
      res = fread(dptr,sizeof(int),2,fp);
      assert(res == 2);
      _corrDraws.insert(_is->resolveDrawWritePtr(dptr));
    }
  }

  res = fread(dptr,sizeof(int),2,fp); // _prevFrame pointer
  assert(res == 2);
  if (dptr[0] == -1)
    _prevFrame = NULL;
  else
    _prevFrame = _is->resolveRotoWritePtr(dptr);

  res = fread(dptr,sizeof(int),2,fp); // _nextFrame pointer
  assert(res == 2);
  if (dptr[0] == -1)
    _nextFrame = NULL;
  else
  _nextFrame = _is->resolveRotoWritePtr(dptr);

  // rotopath _prevCoors intArray
  /*
  res = fread(&dummy,sizeof(int),1,fp); 
  assert(res == 1);
  if (dummy==0)
    _prevCoors = NULL;
  else {
    _prevCoors = new int[dummy];
    fread(_prevCoors,sizeof(int),dummy,fp);
  }

  // rotopath _nextCoors intArray
  res = fread(&dummy,sizeof(int),1,fp); 
  assert(res == 1);
  if (dummy==0)
    _nextCoors = NULL;
  else {
    _nextCoors = new int[dummy];
    fread(_nextCoors,sizeof(int),dummy,fp);
  }
  */

  res = fread(&_lh, sizeof(int),1,fp);  assert(res==1);
  res = fread(&_hh, sizeof(int),1,fp);  assert(res==1);
  if (version==0) {
    res = fread(&_trackw, sizeof(int),1,fp);  assert(res==1); 
    res = fread(&dummy, sizeof(int),1,fp);  assert(res==1); 
  }
  res = fread(&_fixed, sizeof(bool),1,fp);  assert(res==1); // is this version 0 or 1


  if (version>=1) {
    
    res = fread(&_trackEdges, sizeof(bool), 1, fp); assert(res==1); 

    // joints
    res = fread(&dummy,sizeof(int),1,fp);  assert(res == 1); assert(dummy>=0 && dummy<100);
    _bJoints.reserve(dummy);
    for (i=0; i<dummy; ++i) {
      Joint j;
      res = fread(dptr,sizeof(int),2,fp); assert(res==2);
      j._rp = _is->resolveRotoWritePtr(dptr);
      res = fread(&(j._side), sizeof(int), 1, fp); assert(res==1);
      _bJoints.push_back(j);
    }

    res = fread(&dummy,sizeof(int),1,fp);  assert(res == 1); assert(dummy>=0 && dummy<100);
    _eJoints.reserve(dummy);
    for (i=0; i<dummy; ++i) {
      Joint j;
      res = fread(dptr,sizeof(int),2,fp); assert(res==2);
      j._rp = _is->resolveRotoWritePtr(dptr);
      res = fread(&(j._side), sizeof(int), 1, fp); assert(res==1);
      _eJoints.push_back(j);
    }

  }

  if (version>=2) {
    res=fread(_endsFixed,sizeof(bool),2,fp);   // _endsFixed
    assert(res==2);
  }
  

  //setTrackSubsampling(1); //_trackSubsampling); 
}



void RotoPath::notifyDrawDelete(DrawPath* dp) {
  _corrDraws.erase(dp);
}

/*
void RotoPath::setTrans(const double* trans) {
  memcpy(_trans, trans, sizeof(double)*6);
  _trans[0] += 1.;
  _trans[3] += 1.;
  buildGLTransMatrix();
}

void RotoPath::buildGLTransMatrix() const {
  if (!_glMatrix)
    _glMatrix = new double[16];
  memset(_glMatrix, 0, 16*sizeof(double));
  _glMatrix[0] =  _trans[0];
  _glMatrix[1] =  _trans[1];
  _glMatrix[4] =  _trans[2];  //// reflection across y axis for openGL
  _glMatrix[5] =  _trans[3];  //// ditto
  _glMatrix[15] = 1;
  _glMatrix[12] = _trans[4];
  _glMatrix[13] = _trans[5];
}

*/


void RotoPath::actualTranslate(const int n, const Vec2f& t) {
  Vec2f l = getElement(n);
  l += t;
  setElement(n,l);
}

void RotoPath::actualSetLoc(const int n, const Vec2f& l) {
  setElement(n,l);
}

/*
void RotoPath::buildINextCoors() {
  _nextCoors = new int[getNumElements()];
  for (int i=0; i<getNumElements(); i++)
    _nextCoors[i] = i;
}

void RotoPath::buildIPrevCoors() {
  _prevCoors = new int[getNumElements()];
  for (int i=0; i<getNumElements(); i++)
    _prevCoors[i] = i;
}
*/

void RotoPath::buildINextCorrs() {
  //assert(!_nextCont);
  if (_nextCont) delete _nextCont;
  _nextCont = new ContCorr(_bez->numSegs(), getNumElements());
}

void RotoPath::buildIPrevCorrs() {
  //assert(!_prevCont);
  if (_prevCont) delete _prevCont;
  _prevCont = new ContCorr(_bez->numSegs(), getNumElements());
}


// TO DO: get rid of overdraw
void RotoPath::spliceIn(const RotoPath *other, const int i1, const int i2) {
  int l;
  /*
  for (l=0; l<getNumElements(); l++)
    printf("%d: (%f %f)\n", l, getElement(l).x(), getElement(l).y());
  printf("\n");
 
   for (l=0; l<other->getNumElements(); l++)
    printf("%d: (%f %f)\n", l, other->getElement(l).x(), other->getElement(l).y());
  printf("\n");

  printf("%d %d, %d %d\n",getNumElements(), other->getNumElements(),i1,i2);
  */
  /*int oldSize = getNumElements();
  int newsize = i1 + other->getNumElements() + oldSize-i2-1, i,j;
  ensureCapacity(newsize);
  count = MAX(oldSize,newsize);

  for (i=newsize-1,j=oldSize-1; j>i2; i--,j--)
    setElement(i,getElement(j));
  for (i=i1,j=0; j<other->getNumElements(); i++,j++)
    setElement(i,other->getElement(j));
  count=newsize;
  */

  int oldCount = getNumElements();
  Vec2f* oldStuff = new Vec2f[oldCount];
  memcpy(oldStuff, getData(), oldCount*sizeof(Vec2f));
  resetElements();
  for (l=0; l<i1; ++l)
    add(oldStuff[l]);
  for (l=0; l<other->getNumElements(); ++l)
    add(other->getElement(l));
  for (l=i2+1; l<oldCount; ++l)
    add(oldStuff[l]);
  
  delete[] oldStuff;
  
  /*if (_prevCoors) {
    delete[] _prevCoors; _prevCoors=NULL; 
  }
  if (_nextCoors) {
    delete[] _nextCoors; _nextCoors=NULL; 
    }
  */
  
  //for (l=0; l<getNumElements(); l++)
  //printf("%d: (%f %f)\n", l, getElement(l).x(), getElement(l).y());
  //printf("\n");
}

void RotoPath::vacatePrevRotoCorr() {
_prevFrame = NULL;
/*  if (_prevCoors) {
    delete[] _prevCoors; _prevCoors=NULL; 
    }*/
  if (_prevCont) {
    delete _prevCont; _prevCont=NULL; 
  }
}

void RotoPath::vacateNextRotoCorr() {
  _nextFrame=NULL;
  /*if (_nextCoors) {
    delete[] _nextCoors; _nextCoors=NULL; 
    }*/
   if (_nextCont) {
    delete _nextCont; _nextCont=NULL; 
  }
}

void RotoPath::vacateRotoCorr() {
  vacatePrevRotoCorr();
  vacateNextRotoCorr();
}

/*
void RotoPath::generateInverseCoors(const int* inCoors, const int inNum, 
				    int *outCoors, const int outNum) {
  int i=1;
  
  for (int j=0; j<outNum; j++) {
    if (j==inCoors[i]) {
      outCoors[j] = i;
      i++;
      assert(i<inNum);
    }
    else
      outCoors[j] = i-1;
  }
  
}
*/


void RotoPath::nudge(const int i, const Vec2f& loc) {
  AbstractPath::nudge(i, loc);
}

void RotoPath::printInfo() {
  printf("numPoints: %d\n",getNumElements());
  printf("Fixed: %d\n",_fixed);
  printf("prev, this, next: %x %x %x\n",_prevFrame, this, _nextFrame);
  //printf("prev coors, next coors: %x %x\n",_prevCoors, _nextCoors);

  //AbstractPath::print();
  if (_bez) {
    //_bez->print(); //G!
   assert(!_bez->hasZMap());
   assert(!_bez->hasvarMap());
  }
  printf("beg joints: ");
  printJoints(&_bJoints);
  printf("end joints: ");
  printJoints(&_eJoints);
}

Joints& RotoPath::joints(const int side) {
  if (side==0)
    return _bJoints;
  else
    return _eJoints;
}

const Joints& RotoPath::joints(const int side) const {
  if (side==0)
    return _bJoints;
  else
    return _eJoints;
}


void RotoPath::insertJoint(const Joint& j, int w) {
  if (w==0) {
    if (find(_bJoints.begin(), _bJoints.end(), j) == _bJoints.end()) // no duplicates
      _bJoints.push_back(j);
  }
  else {
    if (find(_eJoints.begin(), _eJoints.end(), j) == _eJoints.end()) // no duplicates
      _eJoints.push_back(j);
  }
  //assert(_bJoints.size() < 20 && _eJoints.size() < 20); 
  assert(jointsOk(_bJoints));
  assert(jointsOk(_eJoints));
}

void RotoPath::replaceMeInJoints(RotoPath* newrp, const int i) {
  Joints::iterator c, foundJoint;
  Joints& myjoints = joints(i);
  for (c = myjoints.begin();c != myjoints.end(); ++c) {
    Joints& jnt = c->_rp->joints(c->_side);
    foundJoint = find_if(jnt.begin(), jnt.end(),
			 bind2nd(mem_fun_ref(&Joint::rpEqual), this));
    assert(foundJoint != jnt.end());
    foundJoint->_rp = newrp;
  }
  assert(jointsOk(myjoints));
}

void RotoPath::destroyBegJoints() {
  Joints::iterator c;
  assert(jointsOk(_bJoints));
  for (c = _bJoints.begin();c != _bJoints.end(); ++c)
    c->_rp->removeJoint(this, c->_side);
  _bJoints.clear();
}

void RotoPath::destroyEndJoints() {
  Joints::iterator c;
  assert(jointsOk(_eJoints));
  for (c = _eJoints.begin();c != _eJoints.end(); ++c)
    c->_rp->removeJoint(this, c->_side);
  _eJoints.clear();
}

void RotoPath::destroyJoints() {
  destroyBegJoints();
  destroyEndJoints();
}

void RotoPath::destroyAllCoincidentJoints(const int s) {
  // if joints ok, we know about all involved splines, and just clearing
  // everything will result in no dangling pointers
  assert(allJointsOk());
  Joints& jnt = joints(s);
  Joints::iterator c;
  for (c=jnt.begin(); c!=jnt.end(); ++c)
    c->_rp->joints(c->_side).clear();
  jnt.clear();
}

void RotoPath::removeJoint(RotoPath* rp, const int side) {
  assert(jointsOk(_bJoints));
  assert(jointsOk(_eJoints));
  if (side==0) {    
    Joints::iterator c = find_if(_bJoints.begin(), _bJoints.end(),
				 bind2nd(mem_fun_ref(&Joint::rpEqual), rp));
    if(c!=_bJoints.end())
      _bJoints.erase(c); 
  }
  else {
    Joints::iterator c = find_if(_eJoints.begin(), _eJoints.end(),
				 bind2nd(mem_fun_ref(&Joint::rpEqual), rp));
    if(c!=_eJoints.end())
      _eJoints.erase(c); 
  }

  assert(jointsOk(_bJoints));
  assert(jointsOk(_eJoints));
}

void RotoPath::addToJointSet(JointSet& all, const int side) {
  Joint j(this, side);
  all.insert(j);

  Joints& jnt = joints(side);
  Joints::iterator c;
  JointSet::iterator cf;
  for (c=jnt.begin(); c!=jnt.end(); ++c) {
    cf = find_if(all.begin(), all.end(), 
		 bind2nd(mem_fun_ref(&Joint::rpEqual), c->_rp));

    // rotoPath not yet explored, not in set
    if (cf == all.end()) 
      c->_rp->addToJointSet(all,c->_side);
 
  }
}

int RotoPath::buildBack(TrackGraph* out, const PathV* toTrack, const int side) {
  PathV::const_iterator ok = find(toTrack->begin(), toTrack->end(), this);
  if (ok == toTrack->end()) {  // must be a path to track (already interpolated)
    //out->processConstraints(this, side);
    return 1;
  }
  if (out->pathAlreadyThere(this))   // must not already be on output composite path (not traversing loops)
    return 0;

  out->addPath(this); // have to already add to avoid infinite loop

  Joints* backJoints = (side==1) ? &_bJoints : &_eJoints; // which side are we coming in?
  Joints::iterator backs;
  int res;
  for (backs = backJoints->begin(); backs != backJoints->end(); ++backs) {
    res = backs->_rp->buildBack(out, toTrack, backs->_side);  
    if (res)
      out->processConstraint(this, backs->_rp, !side, backs->_side);
  }
  return 0;
}

int RotoPath::buildForw(TrackGraph* out, const PathV* toTrack, const int side) {
  PathV::const_iterator ok = find(toTrack->begin(), toTrack->end(), this);
  if (ok == toTrack->end()) {  // must be a path to track (already interpolated)
    //out->processConstraints(this, side);
    return 1;
  }
  if (out->pathAlreadyThere(this))   // must not already be on output composite path (not traversing loops)
    return 0;

  out->addPath(this); // have to already add to avoid infinite loop

  Joints* fJoints = (side==1) ? &_bJoints : &_eJoints; // which side are we coming in?
  Joints::iterator forws; int res;
  for (forws = fJoints->begin(); forws != fJoints->end(); ++forws) {
    res = forws->_rp->buildBack(out, toTrack, forws->_side);  
    if (res)
      out->processConstraint(this, forws->_rp, !side, forws->_side);
  }
  return 0;
}

void RotoPath::buildccomp(TrackGraph* out, const PathV* toTrack) {
  int res;
  res = buildBack(out, toTrack, 1); assert(res==0);
  //out.reverse();
  Joints::iterator forws;
  for (forws = _eJoints.begin(); forws != _eJoints.end(); ++forws) {
    res = forws->_rp->buildForw(out, toTrack, forws->_side);  
    if (res)
      out->processConstraint(this, forws->_rp, 1, forws->_side);
  }

  //out->printStructure();
}

void RotoPath::printJoints(const Joints* j) {
  Joints::const_iterator c;
    for (c=j->begin(); c!=j->end(); ++c)
      printf("rp %x side %d, ", c->_rp, c->_side);
    printf("\n");
}


void RotoPath::reconcileJoints(const PathV* flexible) {
  reconcileJoint(&_bJoints, 0, flexible);
  reconcileJoint(&_eJoints, 1, flexible);
}

// flexible are curves that can move; don't want to move pre-existing curves
// rploc is location of owning curve (since not in joints)
// 11/07/03, assumes js is part of this, side is which side of this
void RotoPath::reconcileJoint(Joints* js, const int side, const PathV* flexible) {
  const Vec2f& rploc = _bez->getEnd(side);
  Vec2f tot = rploc, oloc;
  int code = 0; // 0 for not different, 1 for different & flexible, 2 for different & !flexible
  int count = 1;
  Joints::iterator c;

  for (c = js->begin(); c!=js->end(); ++c) {
    if (c->_side==0)
      oloc = c->_rp->getEnd(0);
    else
      oloc = c->_rp->getEnd(1);
    
    if (find(flexible->begin(), flexible->end(), c->_rp) == flexible->end()) { // constraining curve
      code=2;
      tot = oloc;
      break;
    }

    if (oloc != rploc)
      code = 1;
    tot += oloc;

    ++count;
  }

  if (code==0) return;

  if (code==1)
    tot /= float(count); // average location
  
  for (c = js->begin(); c!=js->end(); ++c) {
    //if (c->_side==0)
    // c->_rp->setElement(0, tot);
    //else
    //c->_rp->setElement(c->_rp->getNumElements()-1, tot);
    c->_rp->setEnd(tot, c->_side);
    c->_rp->handleNewBezCtrls();
  }
  //*rploc = tot;
  setEnd(tot, side);
  handleNewBezCtrls();
  
}

void RotoPath::reconcileOneJointToMe(int which) {
  if (which==0)
    reconcileJointToMe(&_bJoints, &(getEnd(0)));
  else if (which==1)
    reconcileJointToMe(&_eJoints, &(getEnd(1)));
  else
    assert(0);
}

void RotoPath::reconcileJointToMe(Joints* js, const Vec2f* rploc) {
  assert(rploc);
  for (Joints::iterator c = js->begin(); c!=js->end(); ++c) {
    c->_rp->setEnd(*rploc, c->_side);
    //if (c->_side==0)
    //c->_rp->setElement(0, *rploc);
    //else
    //c->_rp->setElement(c->_rp->getNumElements()-1, *rploc);
    c->_rp->handleNewBezCtrls();
  }
}

void RotoPath::copyBackJoints(const PathV* subtree) {
  assert(_prevFrame);
  copyBackJoint(&_bJoints, subtree, 0);
  copyBackJoint(&_eJoints, subtree, 1);
}

void RotoPath::copyBackJoint(Joints* js, const PathV* subtree, int side) {
  for (Joints::iterator c = js->begin(); c != js->end(); ++c)
    if (c->_rp->prevC()) {

      Joint j1;
      j1._rp = c->_rp->prevC();
      j1._side = c->_side;
      _prevFrame->insertJoint(j1, side);

      if (find(subtree->begin(), subtree->end(), c->_rp) == subtree->end()) {
	Joint j2;
	j2._rp = _prevFrame;
	j2._side = side;
	c->_rp->prevC()->insertJoint(j2, c->_side);
      }
    }
}

void RotoPath::copyForwardJoints(const PathV* subtree) {
  assert(_nextFrame);
  copyForwardJoint(&_bJoints, subtree, 0);
  copyForwardJoint(&_eJoints, subtree, 1);
}

void RotoPath::copyForwardJoint(Joints* js, const PathV* subtree, int side) {
  for (Joints::iterator c = js->begin(); c != js->end(); ++c)
    if (c->_rp->nextC()) {

      Joint j1;
      j1._rp = c->_rp->nextC();
      j1._side = c->_side;
      _nextFrame->insertJoint(j1, side);

      if (find(subtree->begin(), subtree->end(), c->_rp) == subtree->end()) {
	Joint j2;
	j2._rp = _nextFrame;
	j2._side = side;
	c->_rp->nextC()->insertJoint(j2, c->_side);
      }
    }
}


void RotoPath::copySetToVector(vector<RotoPath*>* v, set<RotoPath*>* s) {
  for (set<RotoPath*>::const_iterator c = s->begin(); c!= s->end(); ++c)
    v->push_back(*c);
}


bool RotoPath::okEnds() const {
  int n = getNumElements();
  if (n<2) return true;
  return ((getElement(0) != getElement(1)) &&
	  (getElement(n-2) != getElement(n-1)));
}


void RotoPath::goBezier(const float tolSqr, const int* numSamples) {  
  _talkFit.clear();
  int count = getNumElements(), i;
  
  Point2* points = new Point2[count];
  //printf("here goes\n");
  for (i=0; i<count; i++) {
    points[i].x = getPointerToElement(i)->x();
    points[i].y = getPointerToElement(i)->y();
    assert(finite(points[i].x) && finite(points[i].y));
    //printf("%f %f\n",points[i].x, points[i].y);
  }
  FitCurves::FitCurve(points, count, tolSqr);
  assert(_talkFit.getSize() > 0);
  //printf("size %d\n",_talkFit.getSize());
  
  //_bez = new BezSpline(&_talkFit, 2., numSamples);
  _bez = new BezSpline(&_talkFit);
  _bez->redoSamples(2., numSamples);
  
  int newCount = _bez->getDiscreteCount();
  printf("From %d samples to %d, %d control points\n",count, newCount, _bez->numCtrls());

  fillFromBez();

  delete[] points;
}

float RotoPath::findClosestT(const Vec2f& loc) {
  assert(_bez);
  return _bez->findClosestT(loc);
 }

void RotoPath::splitSegment(float t) {
  assert(_bez);
  _bez->split(t);
  handleNewBezCtrls();
  //assert(!_nextFrame);
  if (_nextFrame) {
    delete _nextCont; _nextCont = NULL; // G!
    
    //  assert(_nextCont);
    //  _nextCont->splitSamples(t);
  }
  if (_prevFrame) {
    assert(_prevCont && _prevFrame->_nextCont);
    _prevFrame->_nextCont->splitSamples(t);
    //_prevCont->splitDomain();//G!
  }

  buildTouched(true);
}

void RotoPath::fillFromBez() {
  int newCount = _bez->getDiscreteCount();
  if (_tangents) delete[] _tangents;
  _tangents = new Vec2f[newCount];

  resetElements();
  addNElements(newCount);
  for (int i=0; i<newCount; i++) {
    setElement(i,_bez->getDiscreteLoc(i));
    //if (i>0)
    //assert(getElement(i) != getElement(i-1)); 
    _tangents[i] = _bez->getDiscreteTangent(i);
    _tangents[i].Normalize();
  }
  _samplingDirty = false;
}

void RotoPath::handleNewBezCtrls() {
  _bez->calcAllEval2();
  _bez->redoSamples(2.,NULL);
  fillFromBez();
}

void RotoPath::takeBez(const BezSpline* b) {
  if (_bez) delete _bez;
  _bez = new BezSpline(b);
  //*_bez = *b;
  assert(!_bez->hasZMap());
  assert(!_bez->hasvarMap());
  fillFromBez();
  assert(!_bez->hasZMap());
  assert(!_bez->hasvarMap());
}



void RotoPath::buildCan(bool first) {
  _bez = new BezSpline();

  /*
  if (first)
    _bez->setNumSegs(1);
  else
    _bez->setNumSegs(2);
  _bez->addCtrl(Vec2f(40,40));
  if (first) {
    _bez->addCtrl(Vec2f(60,70));
    _bez->addCtrl(Vec2f(80,10));
  }
  else {
    _bez->addCtrl(Vec2f(50,25));
    _bez->addCtrl(Vec2f(60,32.5));
    _bez->addCtrl(Vec2f(70,40));
    _bez->addCtrl(Vec2f(80,47.5));
    _bez->addCtrl(Vec2f(90,55));
  }
  _bez->addCtrl(Vec2f(100,40));
  */

  
  _bez->setNumSegs(1);
  _bez->addCtrl(Vec2f(40,40));
  if (first) {
    _bez->addCtrl(Vec2f(60,70));
    _bez->addCtrl(Vec2f(80,10));
  }
  else {
    _bez->addCtrl(Vec2f(60,10));
    _bez->addCtrl(Vec2f(80,70));
  }
  _bez->addCtrl(Vec2f(100,40));
  

  _bez->calcAllEval2();      
  _bez->redoSamples(2., NULL);
  fillFromBez();
}

void RotoPath::setCan(int t) {
  Vec2f b(60,70-t*10);
  Vec2f c(80,10+t*10);
  _bez->setControl(b,1);
  _bez->setControl(c,2);
  _bez->calcAllEval2();      
  _bez->redoSamples(2., NULL);
  fillFromBez();
}


void RotoPath::addCtrl(const Vec2f& v) {
  assert(_bez);
  _bez->addCtrl(v);
}

void RotoPath::startManualCreation() {
  _bez = new BezSpline();
}

void RotoPath::finishManualCreation() {
  //fillFromBez();
  handleNewBezCtrls();
  //printf("finish manual creation");

  buildTouched(true);
}

void RotoPath::buildTouched(const bool val) {
  _touched.clear();
  for (int i=0; i<_bez->numCtrls(); ++i)
    _touched.push_back(val);
}

void RotoPath::setCtrlTouched(const int which, const bool val) {
  assert(which < (int)_touched.size());
  _touched[which] = val;
}

void RotoPath::setTouchedAll(const bool val) {
  assert(_touched.size() > 0);
  for (std::vector<bool>::iterator c = _touched.begin();
       c != _touched.end(); ++c)
    (*c) = val;
}

void RotoPath::estimateTouched() {
  assert(_bez && _bez->numCtrls() > 0);
  if (_fixed) 
    buildTouched(true);
  else
    buildTouched(false);

  if (endFixed(0)) setCtrlTouched(0, true);
  if (endFixed(1)) setCtrlTouched(_bez->numCtrls()-1, true);
  
    for (std::set<int>::const_iterator c = _internalsFixed.begin(); 
	 c != _internalsFixed.end(); ++c)
      setCtrlTouched(*c, true);
}

void RotoPath::startFinishManualCreation() {
  _bez->finishBuilding();
}



void RotoPath::buildBackReconcileJoints(PathV& paths, int bFrame, int aFrame) {
  int i,j;
  PathV currPaths = paths;
  PathV currPrevPaths = paths;
  PathV::iterator pc;
  for (i=bFrame-1; i>aFrame; i--) { // iterate over inbetween frames

    for (pc = currPrevPaths.begin(), j=0; pc != currPrevPaths.end(); ++pc, ++j) {
      currPrevPaths[j] = (*pc)->prevC();
      assert(currPrevPaths[j]);
    }

    for (pc = currPaths.begin(); pc != currPaths.end(); ++pc) {  // iterate over current curves, build & reconcile
      (*pc)->copyBackJoints(&currPaths);
      assert((*pc)->prevC());
      (*pc)->prevC()->reconcileJoints(&currPrevPaths);  // changed 11/07/03, wierd, should have surfaced before
    }
    
    for (pc = currPaths.begin(), j=0; pc != currPaths.end(); ++pc, ++j)   // step current curves back in time
      currPaths[j] = (*pc)->prevC();
    
  }
}

void RotoPath::buildForwardReconcileJoints(PathV& paths, int aFrame, int bFrame) {
  int i,j;
  PathV currPaths = paths;
  PathV currNextPaths = paths;
  PathV::iterator pc;
  for (i=aFrame; i<bFrame; ++i) { // iterate over inbetween frames

    for (pc = currNextPaths.begin(), j=0; pc != currNextPaths.end(); ++pc, ++j) {
      currNextPaths[j] = (*pc)->nextC();
      assert(currNextPaths[j]);
    }

    for (pc = currPaths.begin(); pc != currPaths.end(); ++pc) {  // iterate over current curves, build & reconcile
      (*pc)->copyForwardJoints(&currPaths);
      assert((*pc)->nextC());
      (*pc)->nextC()->reconcileJoints(&currNextPaths);  // changed 11/07/03, wierd, should have surfaced before
    }
    
    for (pc = currPaths.begin(), j=0; pc != currPaths.end(); ++pc, ++j)   // step current curves forwards in time
      currPaths[j] = (*pc)->nextC();
    
  }
}


void RotoPath::translate(const Vec2f& delta) {
  AbstractPath::translate(delta);
  _bez->translate(delta);
}


RotoPathPair RotoPath::split(const int c) {
  assert(_bez);
  RotoPath *a = new RotoPath(), *b = new RotoPath();
  a->startManualCreation();
  b->startManualCreation();
  for (int i=0; i<_bez->numCtrls(); ++i) {
    if (i<c) { // add to a
      a->addCtrl(*(_bez->getCtrl(i)));
    }
    else if (i>c) { // add to b
      b->addCtrl(*(_bez->getCtrl(i)));
    }
    else {  // split point, add to both
      a->addCtrl(*(_bez->getCtrl(i)));
      b->addCtrl(*(_bez->getCtrl(i)));
    }
  }
  
  a->startFinishManualCreation();
  b->startFinishManualCreation();
  a->finishManualCreation();
  b->finishManualCreation();

  a->_bJoints = _bJoints;
  replaceMeInJoints(a, 0);
  b->_eJoints = _eJoints;
  replaceMeInJoints(b, 1);
  
  return RotoPathPair(a,b);
}


bool jointsOk(const Joints& j) {
  if (j.size() >= 20)
    printf("fuck\n");
  
  //assert(j.size() < 20);
  return (j.size() < 20);
}

bool RotoPath::allJointsOk() {
  Joints::const_iterator c, found;
  for (int s=0; s<2; ++s) {
    Joint toFind(this, s);
    const Joints& jnt = joints(s);
    for (c=jnt.begin(); c!=jnt.end(); ++c) {
      const Joints& oj = c->_rp->joints(c->_side);
      found = find(oj.begin(), oj.end(), toFind);
      if (found == oj.end()) {
	assert(0);
	return false;
      }
    }
  }


  return true;
}

void RotoPath::forgetRegion(const RotoRegion* r) {
  if (_regions[0] == r)
    _regions[0] = NULL;
    if (_regions[1] == r)
    _regions[1] = NULL;
}

void RotoPath::calcLerp() {
  assert(_bez && !_bez->empty());
  int numC = _bez->numCtrls();
  for (int i=0; i<numC; ++i) {
   
    if (touched(i)) continue;  // skip this point
    RotoPath *first=this, *last=this;
    int c0=0, c1=0;

    // get first
    while(first->prevC() && first->prevC()->getNumControls() == numC) {
      first = first->prevC();
      ++c0;
      if (first->touched(i)) break;
    } 
    if (first == this) continue;   // skip this point

    // get last
    while(last->nextC() && last->nextC()->getNumControls() == numC) {
      last = last->nextC();
      ++c1;
      if (last->touched(i)) break;
    } 
    if (last == this) continue;   // skip this point

    assert(c0 > 0 && c1 > 0);
    Vec2f firstV = *(first->_bez->getCtrl(i));
    Vec2f lastV  = *(last->_bez->getCtrl(i));
    float l = float(c0) / float(c0+c1);
    assert(l>0 && l<1);
    Vec2f newVal;
    Vec2f_LinInterp(newVal, firstV, lastV, l);
    _bez->setControl(newVal, i);

  } // loop over points

  handleNewBezCtrls();
}


Vec2i RotoPath::getStats() const {
  Vec2i res(0,0);
  // handle 0
  // iterate over joints, if my rotoPath pointer is lowest in memory, use, otherwise skip
  const Joints& jts = joints(0);
  Joints::const_iterator c;
  bool touched = jointTouched(0);
  bool least = true;
  for (c=jts.begin(); c!=jts.end(); ++c) {
    //printf("%x %x %d\n",this, c->_rp, c->_rp < this);
    if (c->_rp < this) {
      least = false; break; 
    }
  }
  if (touched && least)
    res.Inc(1,1);
  else if (!touched && least)
    res.Inc(0,1);
    
  // handle end
  const Joints& jts1 = joints(1);
  touched = jointTouched(1);
  least = true;
  for (c=jts1.begin(); c!=jts1.end(); ++c) {
    if (c->_rp < this) {
      least = false;break;
    }
  }
  if (touched && least)
    res.Inc(1,1);
  else if (!touched && least)
    res.Inc(0,1);

  //  handle middle
  for (int i=1; i<_touched.size()-1; ++i) {
    if (this->touched(i))
      res.Inc(1,1);
    else
      res.Inc(0,1);
  }
		
  return res;
}


bool RotoPath::touched(const int i) const {
  assert(_bez->numCtrls() == (int)_touched.size());
  assert(i>=0 && i<(int)_touched.size());
  if (i==0) return (_touched[i] || jointTouched(0));
  if (i==_touched.size()-1) return (_touched[i] || jointTouched(1));
  return (_touched[i]);
}

bool RotoPath::endTouched(const int side) const {
  if (side==0) return _touched[0];
  else return _touched[_touched.size()-1];
}

bool RotoPath::jointTouched(const int side) const {
  const Joints& jts = joints(side);
  Joints::const_iterator c;
  bool touched = false;
  for (c=jts.begin(); c!=jts.end(); ++c) {
    touched = c->_rp->endTouched(c->_side);
    if (touched) return true;
  }
  return touched;
}

int RotoPath::pickCtrl(const float x, const float y) {
  Vec2f pt(x,y);
  int which;
  double dist  = distanceToCtrls2(pt, &which);
  if (dist < 25)
    return which;
  else
    return -1;
}

////   transformSelected iterators   /////////////

SplineCtrlIterator RotoPath::transformSelectedCtrlIterator() {
  SplineCtrlIterator it;
  it._ctrls = &(_bez->getControls());
  it._selected = &(_transformSelected);
  if (_transformSelected.size() == 0)
    it._i = it._ctrls->size(); // at end
  else {
    it._j = it._selected->begin();
    it._i = 3 * *(it._j);
  }
  return it;
}

SplineCtrlIterator& SplineCtrlIterator::operator++() {
  if ( _i<_ctrls->size() ) {
    ++_i;
    if (_i%3==0 && _i < _ctrls->size()-1) {
      ++_j;
      if (_j==_selected->end())
	_i = _ctrls->size(); // put at end
      else
	_i = 3 * *_j;
    }
  }
  return *this;
}

bool SplineCtrlIterator::end() const {
  if (_i>= _ctrls->size()) return true;
  else return false;
}

const Vec2f& SplineCtrlIterator::ctrl() const {
  assert(_i < _ctrls->size());
  return (*_ctrls)[_i];  
}

Vec2f& SplineCtrlIterator::ctrl() {
  assert(_i < _ctrls->size());
  return (*_ctrls)[_i];
}

SplineSampleIterator RotoPath::transformSelectedSampleIterator() {
  SplineSampleIterator it;
  it._samples = &(_bez->getSamples());
  it._selected = &(_transformSelected);
  if (_transformSelected.size() == 0)
    it._i = it._samples->size(); // at end
  else {
    it._j = it._selected->begin();
    it._i = -1;
    ++it;
  }

  return it;
}

bool SplineSampleIterator::end() const {
  if (_i >= (int)_samples->size()) return true;
  else return false;
}

SplineSampleIterator& SplineSampleIterator::operator++() {
  if (_i >= (int)_samples->size()) return *this;

  bool startedOk = (int((*_samples)[_i]._baset) == *_j);
  ++_i;

  // stepped out of range, increment _j
  if (startedOk  && 
      int((*_samples)[_i]._baset) != *_j)
    ++_j;
  
  // if not in selected, keep going
  while (int((*_samples)[_i]._baset) != *_j &&
	 _i < (int)_samples->size())
    ++_i;  

  return *this;
}
const DSample& SplineSampleIterator::sample() const {
  assert(!end());
  return (*_samples)[_i];
}
DSample& SplineSampleIterator::sample() {
  assert(!end());
  return (*_samples)[_i];
}


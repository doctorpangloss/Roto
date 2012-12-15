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



#include "transformer.h"

Transformer::Transformer() { 

}

// adds segment to selected
void Transformer::addSplineSegment(RotoPath* rp, float t) {
  assert(rp);
  int base = (int) t;
  rp->transformSelect(base);
  _selected.insert(rp);
  //printf("added %d to transformed\n", (int)t);

  base *= 3;
  for (int i=0; i<4; ++i)
    _bbox.includePoint(rp->getCtrlRef(base+i));

  _centroid.Set((_bbox.lower.x() + _bbox.upper.x())/2.,
		(_bbox.lower.y() + _bbox.upper.y())/2.);
}

void Transformer::addSpline(RotoPath* rp) {
  assert(rp);
  _selected.insert(rp);
  
  int ns = rp->getNumSegs(), base;
  for (int s=0; s<ns; ++s) {
    rp->transformSelect(s);
    base = 3*s;
    for (int i=0; i<4; ++i)
      _bbox.includePoint(rp->getCtrlRef(base+i));
  }

  _centroid.Set((_bbox.lower.x() + _bbox.upper.x())/2.,
		(_bbox.lower.y() + _bbox.upper.y())/2.);
}
  
void Transformer::recalcBBox() {
  _bbox.clear();
  for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c) {
    SplineCtrlIterator c2 = (*c)->transformSelectedCtrlIterator();
    for (; !c2.end(); ++c2) {
      //printf("bboc %d\n",c2._i);
      _bbox.includePoint(c2.ctrl());
    }
  }
}

// true if there is something to transform
bool Transformer::transformable() const { 
  return (!_selected.empty()); }

// true if currently transforming something
bool Transformer::transforming() const { 
  return _transforming; }


// used to determine transform mode (rotate, scale), starting config, start transform
void Transformer::startTransform(const Vec2f& loc) {
  _transforming = true;
  if (_centroid.distanceTo2(loc) < 25) {
    _transformPoint = _centroid;
    _mode = M_TRANSLATE;
  }

  else if (_bbox.upper.distanceTo2(loc) < 25) {
    _transformPoint.Set(_bbox.upper.x(), _bbox.upper.y());
    _mode = M_UNIFORMSCALE;
    _last = 1;
  }

  else {
    _transformPoint = loc;
    _mode = M_ROTATE;
    _last = 0;
  }

  //printf("loc: %f %f, lowerRight: %f %f\n",loc.x(), loc.y(), lowerRight.x(), lowerRight.y());
  printf("mode: %d\n",_mode);
}


// during transformation, mouse moved to new loc
void Transformer::handleNewLoc(const Vec2f& loc) {
  if (_mode==M_UNIFORMSCALE) {
    float mag=0;
    //Vec2f lowerRight(_bbox.upper.x(), _bbox.upper.y());
    Vec2f pc(loc, _centroid);
    Vec2f lrc(_transformPoint, _centroid);
    mag = pc.Dot2(lrc) / lrc.Len2();
    mag = MAX(mag,0);
    //printf("magFactor %f\n",mag);

    for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
      (*c)->uniformScaleSelectedAbout(1. + mag-_last, _centroid);
    _last = mag;
    recalcBBox();
  }

  else if (_mode==M_TRANSLATE) {
    Vec2f delta(loc, _transformPoint);
    for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
      (*c)->translateSelected(delta);
    _transformPoint = loc;
    recalcBBox();
  }

  else if (_mode == M_ROTATE) {
    Vec2f a(loc, _centroid);
    Vec2f b(_transformPoint, _centroid);
    float costheta = a.Dot2(b) / (a.Len() * b.Len());
    float sintheta = a.Cross2(b) / (a.Len() * b.Len());
    float theta = atan2(sintheta, costheta);
    for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
      (*c)->rotateSelectedAbout(-(theta-_last), _centroid);
    _last = theta;
    recalcBBox();
  }
  
}

void Transformer::endTransform(const Vec2f& loc) {
  _transforming = false;
  for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
    (*c)->handleNewBezCtrls();

  recalcBBox();
  _centroid.Set((_bbox.lower.x() + _bbox.upper.x())/2.,
		(_bbox.lower.y() + _bbox.upper.y())/2.);

  for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
    (*c)->setTouchedAll(true);
}

void Transformer::clear() { 
  for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
    (*c)->clearTransformSelected();
  _selected.clear(); 
  _transforming = false; 
  _bbox.clear();
}

void Transformer::render() {
  //printf("transformer::render\n");
  for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
    (*c)->renderTransformSelected(); 

  // print bounding box
  if (_bbox.empty()) return;
  glColor3f(0,1,0);
  glBegin(GL_LINE_STRIP);
  glVertex2f(_bbox.lower.x(), _bbox.lower.y());
  glVertex2f(_bbox.upper.x(), _bbox.lower.y());
  glVertex2f(_bbox.upper.x(), _bbox.upper.y());
  glVertex2f(_bbox.lower.x(), _bbox.upper.y());
  glVertex2f(_bbox.lower.x(), _bbox.lower.y());
  glEnd();

  if (!transforming()) {
    glPushAttrib(GL_POINT_BIT);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    glVertex2f(_centroid.x(), _centroid.y());
    glEnd();
    glPopAttrib();
  }
  
}

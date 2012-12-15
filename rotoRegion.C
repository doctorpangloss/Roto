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



#include "rotoRegion.h"

RotoRegion::RotoRegion(const PathV& v) {
  for (PathV::const_iterator c = v.begin(); c!= v.end(); ++c)
    _rotoPaths.push_back(*c);
  assert(_rotoPaths.size() > 1);

  _sides.push_back(RR_RIGHT);  // we start first curve 0->1
  const RotoPath* last = _rotoPaths.front();
  _rotoPaths.front()->setRegion(this, RR_RIGHT);

  int i=1;
  PathV::const_iterator c = _rotoPaths.begin();
  ++c; // move one up
  float d0, d1;
  for (; c!= _rotoPaths.end(); ++c, ++i) {
    const Vec2f endPt   = (*c)->getEnd( getLast(i-1)  );
    d0 = (*c)->_bez->distanceToEnd2(endPt, 0);
    d1 = (*c)->_bez->distanceToEnd2(endPt, 1);
    if (d0 < d1)
      _sides.push_back(RR_RIGHT);
    else
      _sides.push_back(RR_LEFT);
    (*c)->setRegion(this, _sides[i]);

    last = *c;
  }
  assert(_sides.size() == _rotoPaths.size());

  createMask();
}

void RotoRegion::createMask() {
  _mask = new unsigned char[_w*_h];
  _listNum = glGenLists(1);
  assert(_listNum != 0);
  calculateMask();
}

RotoRegion::RotoRegion() {
  _mask = NULL;
  _listNum = 0;
}

RotoRegion::~RotoRegion() {
  PathV::const_iterator c;
  for (c = _rotoPaths.begin(); c!= _rotoPaths.end(); ++c  )
    (*c)->forgetRegion(this);
  if (_mask) delete[] _mask;
}


void RotoRegion::render() const {
  glColor3f(0,0,1);
  glCallList(_listNum);
}

void rrfuck(GLenum errno) {
  const GLubyte *errString = gluErrorString(errno);
  printf("OpenGL error: %s\n",errString);
  assert(0);
}
void rrfuckCombine(GLdouble coords[3], GLdouble* vd[4], GLfloat weight[4], void **outData) {
  GLdouble *vertex = (GLdouble *) malloc(3*sizeof(GLdouble));
  memcpy(vertex, coords, 3*sizeof(GLdouble));
  *outData = vertex;

  //printf("vertex data\n");
  //for (int i=0; i<4; i++)
  //printf("%f %f %f\n",vd[i][0], vd[i][1], vd[i][2]);
  //assert(0);
}

void RotoRegion::calculateMask() {

  int numPoints = 0, d;
  PathV::const_iterator c;
  for (c = _rotoPaths.begin(); c!= _rotoPaths.end(); ++c  ) {
    d = (*c)->_bez->getDiscreteCount();
    assert(d>0);
    numPoints += d;
  }

  GLdouble *v = new GLdouble[numPoints*3];
  Vec2f loc;
  Vec2f endPt;
  int gi, l;
  float d0, d1;
  std::vector<short>::const_iterator c2 = _sides.begin();
  for (l=0, gi=0, c = _rotoPaths.begin(); c!= _rotoPaths.end(); ++c, ++c2, ++l ) {


    
    d0 = (*c)->_bez->distanceToEnd2(endPt, 0);
    d1 = (*c)->_bez->distanceToEnd2(endPt, 1);
    
    if ((l !=0 && d0 < d1) || (l==0 && *c2 == RR_RIGHT)) {
      for (int i=0; i<(*c)->_bez->getDiscreteCount(); ++i, gi+=3) {
	loc = (*c)->_bez->getDiscreteLoc(i);
	v[gi] = loc.x(); v[gi+1] = loc.y(); v[gi+2] = 0;
	printf("Right: %f %f\n",loc.x(), loc.y());
	if ((*c)->_bez->getDiscreteCount()-1 == i)
	  endPt = loc;
      }
    }

    else {
      for (int i=(*c)->_bez->getDiscreteCount()-1; i>=0; --i, gi+=3) {
	loc = (*c)->_bez->getDiscreteLoc(i);
	v[gi] = loc.x(); v[gi+1] = loc.y(); v[gi+2] = 0;
	printf("Left: %f %f\n",loc.x(), loc.y());
	if (i==0)
	  endPt = loc;
      }
    }

  }

  // render tesselated version, make it render into stencil buffer but fail depth test
  // then, read stencil buffer into mask
  assert(glGetError() == GL_NO_ERROR);
  glEnable(GL_STENCIL_TEST);
  glEnable(GL_DEPTH_TEST);
  glClear(GL_STENCIL_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glStencilFunc(GL_ALWAYS, 1, 0x1);
  glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);
  glDepthFunc(GL_NEVER);
  assert(glGetError() == GL_NO_ERROR);

  glNewList(_listNum, GL_COMPILE_AND_EXECUTE);

  GLUtesselator* tobj = gluNewTess();
  gluTessCallback(tobj, GLU_TESS_VERTEX, 
		  (GLvoid (*) ()) &glVertex3dv);
  gluTessCallback(tobj, GLU_TESS_BEGIN, 
		  (GLvoid (*) ()) &glBegin);
  gluTessCallback(tobj, GLU_TESS_END, 
		  (GLvoid (*) ()) &glEnd);
  gluTessCallback(tobj, GLU_TESS_ERROR, 
		  (GLvoid (*) ()) &rrfuck);
  gluTessCallback(tobj, GLU_TESS_COMBINE, 
		  (GLvoid (*) ()) &rrfuckCombine);

  gluTessBeginPolygon(tobj, NULL);
  gluTessBeginContour(tobj);
  
  int i;
  for (i=0; i<numPoints; i++)
    gluTessVertex(tobj, v+i*3, v+i*3);

  gluTessEndContour(tobj);
  gluTessEndPolygon(tobj);
  gluDeleteTess(tobj);

  glEndList();
  assert(glGetError() == GL_NO_ERROR);
  delete[] v;

  glReadPixels(0,0,_w, _h, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, _mask);

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_STENCIL_TEST);

  // flip the stencil upside down
  int j=0, bj=_h-1;
  while (j<bj) {
    for (i=0; i<_w; ++i) 
      std::swap(_mask[j*_w + i], _mask[bj*_w + i]);
    ++j; --bj;
  }

  /*
  QImage im (_w, _h, 32);
  int index=0;
  for (j=0; j<_h; ++j)
    for (i=0; i<_w; ++i, ++index)
      im.setPixel(i, j,qRgb(_mask[index]*255,0,0));
  im.save("out.png", "PNG");
  */
}

int RotoRegion::getLast(const int i) {
  if (_sides[i] == RR_RIGHT)
    return 1;
  else
    return 0;
}

int RotoRegion::getFirst(const int i) {
  if (_sides[i] == RR_RIGHT)
    return 0;
  else
    return 1;
}

bool RotoRegion::withinRegion(const float x, const float y) {
  const int xx = int(x), yy = int(y);
  return (_mask[yy*_w + xx]);
}

void RotoRegion::renderMask(uchar* bits) const {
  int index=0, i, j, i4=0;
  for (j=0; j<_h; ++j)
    for (i=0; i<_w; ++i, ++index, i4+=4) {
      if (_mask[index]) {
	bits[i4]   = 0;
	bits[i4+1] = 0;
	bits[i4+2] = 0;
      }
    }
}

RotoRegion* RotoRegion::tryToCopyForward() {
  RotoRegion* t = new RotoRegion();

  PathV::const_iterator c;
  int i=0;
  for (c = _rotoPaths.begin(); c!= _rotoPaths.end(); ++c, ++i  ) {  
    if ((*c)->nextC()) {
      t->_rotoPaths.push_back((*c)->nextC());
      (*c)->nextC()->setRegion(t, _sides[i]);
    }
    else {
      delete t;
      return NULL;
    }
  }
  
  t->_sides = _sides;
  assert(t->_sides.size() == t->_rotoPaths.size());
  t->createMask();

  return t;
}

unsigned char* RotoRegion::copyMask() const {
  unsigned char* res = new unsigned char[_w*_h];
  memcpy(res, _mask, _w*_h*sizeof(unsigned char));
  return res;
}

void RotoRegion::orMask(unsigned char* omask) const {
  int index=0, i, j;
  for (j=0; j<_h; ++j)
    for (i=0; i<_w; ++i, ++index) {
      omask[index] = ( omask[index] || _mask[index] );
    }
}

int RotoRegion::_w;
int RotoRegion::_h;

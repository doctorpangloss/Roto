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



#include "drawPatch.h"
#include "bBoxf2D.h"
//#include <math.h>

// this better only used for loading
DrawPatch::DrawPatch() {
   _glMatrix = NULL;
   _justDrawn = 0;
   _stroke=NULL; 
   _dlist = 0;
}

DrawPatch::DrawPatch(DrawPatch* other, const double* trans = NULL) {
  _x = other->_x; _y = other->_y;
  _w = other->_w; _h = other->_h;
  _stroke = new Stroke(*(other->_stroke));

  // DEBUG
  printf("DrawPatch.C constructor : made new stroke \n") ; 
  

  _glMatrix = NULL;
  _dlist = 0;

  _nextFrame = _prevFrame = NULL;
  _lastTranslation = other->_lastTranslation;
  _justDrawn = 1; // this is not accurate for copy construction
  if (trans) {

    _trans[0] = (1.+trans[0])*other->_trans[0] + trans[1]*other->_trans[2];
    _trans[1] = (1.+trans[0])*other->_trans[1] + trans[1]*other->_trans[3];
    _trans[2] = trans[2]*other->_trans[0] + (1.+trans[3])*other->_trans[2];
    _trans[3] = trans[2]*other->_trans[1] + (1.+trans[3])*other->_trans[3];
    _trans[4] = other->_trans[4]+trans[4]; 
    _trans[5] = other->_trans[5]+trans[5];
  }
  else {               
    memcpy(_trans, other->_trans, sizeof(double)*6);
  }
}

DrawPatch::DrawPatch(const int x, const int y, const int w, const int h, Vec3f color) :
  _w(w), _h(h), _color(color) {
  _x = float(x);
  _y = float(y);
  _stroke = NULL;
  _glMatrix = NULL;
  //_mask = NULL;
  //_maskGL = NULL;
  _nextFrame = _prevFrame = NULL;
  _justDrawn = 1; 
  _trans[0] = _trans[3] = 1.;
  _trans[1] = _trans[2] = 0.;
  _trans[4] = _trans[5] = 0;
  _dlist = 0;

}

DrawPatch::~DrawPatch() {
  if (_glMatrix) delete[] _glMatrix;
}

void DrawPatch::setTrans(const double* trans) {
  memcpy(_trans, trans, sizeof(double)*6);
  _trans[0] += 1.;
  _trans[3] += 1.;
  buildGLTransMatrix();
}

void DrawPatch::initDummyStroke() {
  _stroke = new Stroke();
  
  /** may have to init texture info here... see drawPath.h ***/

  printf("DrawPatch.C initDummmyStroke : made new stroke \n") ; 

  _stroke->radius() = 3;
  _stroke->addControlPoint(_x,  _y, 3);
  _stroke->addControlPoint(_x+_w,  _y, 3);
  _stroke->addControlPoint(_x+_w,  _y+_h, 3);
  _stroke->addControlPoint(_x,  _y+_h, 3);
   _stroke->addControlPoint(_x,  _y, 3);
}

void DrawPatch::buildGLTransMatrix() const {
  if (!_glMatrix)
    _glMatrix = new double[16];
  memset(_glMatrix, 0, 16*sizeof(double));
  _glMatrix[0] =  _trans[0];
  _glMatrix[1] =  _trans[1];
  _glMatrix[4] =  _trans[2]; 
  _glMatrix[5] =  _trans[3]; 
  _glMatrix[15] = 1;
  _glMatrix[12] = _trans[4];
  _glMatrix[13] = _trans[5];
}

void DrawPatch::render() const {

  if (_dlist == 0) 
    calculateDisplayList();
  glColor4f(_color.r(), _color.g(), _color.b(), _dAlpha);
  glCallList(_dlist);

}

void DrawPatch::calculateDisplayList() const {

  if (_dlist == 0)
    _dlist = glGenLists(1);
  assert(glGetError() == GL_NO_ERROR);
  assert(_dlist != 0);
  
  glNewList(_dlist, GL_COMPILE);

  Vec2f myCenter = center();
  glPushMatrix();

  glTranslatef(myCenter.x(), myCenter.y(),0);

  if (!_glMatrix)
    buildGLTransMatrix();
  glMultMatrixd(_glMatrix);  

  glTranslatef(-myCenter.x(), - myCenter.y(),0);

  _stroke->render();

  /*
  // render crosshairs
  glColor3f(0,0,1);
  glTranslatef(myCenter.x(), myCenter.y(),0);
  glBegin(GL_LINES);
  glVertex3i(-8,0,0);
  glVertex3i(8,0,0);
  glVertex3i(0,-8,0);
  glVertex3i(0,8,0);
  glEnd();
  */

  glPopMatrix();

  glEndList();
  assert(glGetError() == GL_NO_ERROR);  
}

void DrawPatch::renderSelected() const {
  Vec2f myCenter = center();
  glPushMatrix();
  
  glTranslatef(myCenter.x(), myCenter.y(),0);

  if (!_glMatrix)
    buildGLTransMatrix();
  glMultMatrixd(_glMatrix);  
  
  glTranslatef(-myCenter.x(), - myCenter.y(),0);
  

  Vec2f loc(_x,_y);
  glRectf(loc.x()-2, loc.y()-2, 
	  loc.x()+2, loc.y()+2);
  loc += Vec2f(_w,_h);
  glRectf(loc.x()-2, loc.y()-2, 
	  loc.x()+2, loc.y()+2);

  glPopMatrix();
}

void DrawPatch::gutMe(const DrawPatch* other) {
  _x = other->_x; _y = other->_y;
  _w = other->_w; _h = other->_h;
  if (_stroke) delete _stroke;
  _stroke = new Stroke(*(other->_stroke));

  /***** may have to init texture info stuff here - see drawPath.h ***/

  _glMatrix = NULL;
}


void DrawPatch::save(FILE *fp) const {
  fwrite(&_w, sizeof(int), 1, fp); // w
  fwrite(&_h, sizeof(int), 1, fp); // h
  fwrite(&_x, sizeof(int), 1, fp); // x
  fwrite(&_y, sizeof(int), 1, fp); // y
  fwrite(&_color, sizeof(Vec3f), 1, fp); //_color
  fwrite(_trans, sizeof(double), 6, fp); // _trans

  if (_prevFrame)    // _prevFrame
    _prevFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  if (_nextFrame)    // _nextFrame
    _nextFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  assert(_stroke);
  _stroke->save(fp);
}

void DrawPatch::load(FILE *fp) {
  int res, dptr[2];
  res = fread(&_w, sizeof(int), 1, fp); // w 
  assert(res==1);
  res = fread(&_h, sizeof(int), 1, fp); // h
  assert(res==1);
  res = fread(&_x, sizeof(int), 1, fp); // x
  assert(res==1);
  res = fread(&_y, sizeof(int), 1, fp); // y
  assert(res==1);
  res = fread(&_color, sizeof(Vec3f), 1, fp); //_color
  assert(res==1);
  res = fread(_trans, sizeof(double), 6, fp); // _trans
  assert(res==6);

  res = fread(dptr,sizeof(int),2,fp); // _prevFrame pointer
  assert(res == 2);
  if (dptr[0] == -1)
    _prevFrame = NULL;
  else
    _prevFrame = _is->resolveDrawPatchWritePtr(dptr);
  
  res = fread(dptr,sizeof(int),2,fp); // _nextFrame pointer
  assert(res == 2);
  if (dptr[0] == -1)
    _nextFrame = NULL;
  else
    _nextFrame = _is->resolveDrawPatchWritePtr(dptr);

  _stroke = new Stroke(fp);

  printf("DrawPatch.C load : made new stroke \n") ; 

  /*** may have to init texture data here
   *** see drawPath.C ... */


}


void DrawPatch::writePtr(FILE* fp) const {
  fwrite(_ptrDoc,sizeof(int),2,fp);
}

void DrawPatch::writeNullPtr(FILE* fp) {
  int dummy[2] = {-1,-1};
  fwrite(&dummy,sizeof(int),2,fp);
}

int DrawPatch::_globalh;
float DrawPatch::_dAlpha;
ImgSequence* DrawPatch::_is;

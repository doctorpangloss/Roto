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



#ifndef DRAWPATCH_H
#define DRAWPATCH_H

#ifndef uchar
typedef unsigned char uchar;
#endif

#include <stdio.h>
#include <GL/gl.h>
#include <qimage.h>
#include "hertzStrokes/Stroke.h"
#include "imgSequence.h"
#include "jl_vectors.h"

class DrawPatch {
  
 public:
  DrawPatch(); // only for loading from files

  DrawPatch(const int x, const int y, const int w, const int h, Vec3f color);

  DrawPatch(DrawPatch* other, const double* trans = NULL);
  // If trans is given, will set out trans to transformed matrix
  // of other, by premultiplying other's _trans with new trans (after adding identity to trans);
  
  void takeStroke(Stroke* st) { _stroke = st; }

  void setTrans(const double* trans);
  const double* getTrans() const { return _trans; }
  const void getTrans(double* here) const { memcpy(here,_trans,6*sizeof(double)); }

  ~DrawPatch();

  /*
  void prepareCapture() {
    printf("preparing to cpature\n");
    glEnable(GL_DEPTH_TEST);
    
    }*/

  //void capture();

  //void saveImage(char* filename);

  void render() const;

  void renderSelected() const;

  //bool onPatch(const int x, const int y) const;

  float x() const { return _x;}
  float y() const { return _y; }
  int width() const { return _w;}
  int height() const { return _h; }
  Vec2f center() const { return Vec2f(_x+float(_w)/2.f, _y+float(_h)/2.f); }
  //Vec2f glCenter() const { return Vec2f(_glx+float(_glw)/2.f, _gly+float(_glh)/2.f); }

  int justDrawn() const { return _justDrawn; }
  
  void setJustDrawn(int val) { _justDrawn = val; }

  /*
  void translate(const float dx, const float dy, const bool setLast=true) { 
    _x+=dx; _y+=dy;
    _glx += dx; _gly -= dy;
    if (setLast)
      _lastTranslation.Set(dx,dy);
  }
  */
  
  void setLastTranslate(const float dx, const float dy) {
     _lastTranslation.Set(dx,dy);
  }

  void setTransTranslate(const float dx, const float dy) {
    _trans[4] = dx; _trans[5] = dy;
    buildGLTransMatrix();
  }

  void setXY(const float x, const float y) {
    _x = x; _y = y;
  }
  void setCenter(const float x, const float y) {
    _x = x - _w/2.f;
    _y = y - _h/2.f;
  }

  void save(FILE *fp) const;
  void load(FILE *fp);

  Vec2f getLastTranslation() const {
    return _lastTranslation; }

  //void initToBlack();

  void initDummyStroke();

  void gutMe(const DrawPatch* other);

  // pointer documentation
  int* getPtrDoc() { return _ptrDoc;}
  const int* getPtrDoc() const { return _ptrDoc;}
  void setPtrDoc(const int frame, const int num) {
    _ptrDoc[0] = frame; _ptrDoc[1] = num; }
  void writePtr(FILE* fp) const;

  DrawPatch* nextP() { return _nextFrame; }
  void setNextP(DrawPatch* next) {
    _nextFrame = next; }
  DrawPatch* prevP() { return _prevFrame; }
  void setPrevP(DrawPatch* prev) {
    _prevFrame = prev; }

  DrawPatch *prev, *next;

  static int _globalh;
  static float _dAlpha;
  static ImgSequence* _is;

 private:
  void calculateDisplayList() const;

  void transform(Vec2f* out, const Vec2f loc) const {
    out->set_x( _trans[0]*loc.x() + _trans[2]*loc.y() );
    out->set_y( _trans[1]*loc.x() + _trans[3]*loc.y() );
  }
  //void buildGL();
  void buildGLTransMatrix() const;

  static void writeNullPtr(FILE* fp);

  int _w, _h;
  float _x, _y;

  Vec3f _color;
  DrawPatch *_prevFrame, *_nextFrame; 
  double _trans[6]; // affine transformation around center
                    // careful, this is column-major 
  Stroke *_stroke;
  
  mutable double *_glMatrix;

  Vec2f _lastTranslation;
  int _justDrawn;
  // pointer documentation
  int _ptrDoc[2]; // frame number, num in list
  mutable int _dlist;
};

#endif

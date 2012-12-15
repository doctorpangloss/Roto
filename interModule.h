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



#ifndef INTERMODULE_H
#define INTERMODULE_H

#ifdef __APPLE__
using namespace std;
#endif
#include <qwidget.h>
#include <string>
#include <map>
#include <GL/gl.h>
#include <GL/glu.h>
#include "jl_vectors.h"
#include "dynArray.h"

typedef std::map<int,std::string> CmdMap;
typedef std::map<int, bool> EnabledMap;

class InterModule : public QObject {
  //abstract class

Q_OBJECT

 public:

  virtual ~InterModule() {}

  virtual void paintGL() = 0;
  
  virtual void mousePressEvent( QMouseEvent *e ) = 0;
  virtual void mouseReleaseEvent( QMouseEvent *e ) = 0;
  virtual void mouseMoveEvent( QMouseEvent *e ) = 0;
  virtual void keyPressEvent ( QKeyEvent *e ) = 0;

  virtual void frameChange(int i) = 0;

  virtual CmdMap* getPopupFunctions(const Vec2i& loc, DynArray<int,5>* enab) = 0;
  virtual void callPopup(const int which) = 0;

  void mutualInit();
  void toggleShow(int state);

  virtual void reportPressure(const int p) {}

  virtual void runCan() = 0;

  bool getRange(int& a, int& b);
  
  bool getMinimalRenderMode() const {
    return _minimalRenderMode;
  }

  Vec2f unproject(const int x, const int y, const int h) const; // need to reverse y

  void setPropMode(int s) { _propMode = s;}

  void shiftChanged(const bool in) { _shift = in; }
  void controlChanged(const bool in) { _control = in; }

 signals:
  void selectChanged();

 protected:

  bool _showAllCorr;
  bool _minimalRenderMode;
  int _propMode; // 0 for forward, 1 for keyframe
  bool _shift; // true, it's down, false it's up
  bool _control; // true, it's down, false it's up
};



#endif

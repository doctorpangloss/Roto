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



#include "interModule.h"
#include "rangeDialog.h"


void InterModule::mutualInit() {
  _showAllCorr = false;
  _minimalRenderMode = 0;
  _propMode = 1;
  _shift = false;
}

void InterModule::toggleShow(int state) {
  if (state==2)
    _showAllCorr = true;
  else
    _showAllCorr = false;
}


Vec2f InterModule::unproject(const int xx, const int yy, const int h) const {
  float x = float(xx) + .5, y = float(yy) + .5;
  GLdouble proj[16], model[16];
  GLint view[4];
  GLdouble ox, oy, oz;
  Vec2f result;

  glGetDoublev(GL_PROJECTION_MATRIX, proj);
  glGetDoublev(GL_MODELVIEW_MATRIX, model);
  glGetIntegerv(GL_VIEWPORT, view);
  
  int res = gluUnProject(x, h-y, 0, model, proj, view, &ox, &oy, &oz);
  assert(res == GL_TRUE);
  result.Set(ox, oy);
  return result;  
}


bool InterModule::getRange(int& a, int& b) {
  Vec2i range;
  RangeDialog rd(NULL, a,b);
  if (rd.exec() == QDialog::Accepted) {
    range = rd.getRange();
    a = range.x();
    b = range.y();
    return true;
  }
  else
    return false;
}


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



#ifndef BBOXF2D_H
#define BBOXF2D_H

#include <limits.h>
#include <float.h>
#include <stdio.h>
#include "jl_vectors.h"

class Bboxf2D {
public:

  Bboxf2D() {
    lower.Set(FLT_MAX,FLT_MAX);
    upper.Set(FLT_MIN, FLT_MIN);
  }

  Bboxf2D(const Bboxf2D& other) {
    upper = other.upper; lower = other.lower;
  }

  void clear() {
    lower.Set(FLT_MAX,FLT_MAX);
    upper.Set(FLT_MIN, FLT_MIN);
  }

  bool empty() { 
    return (lower.x() > upper.x() || lower.y() > upper.y());
  }

  void setToPoint(const Vec2f& pt) {
    lower = pt; upper = pt;
  }

  void setToUnion(const Bboxf2D& a, Bboxf2D& b) {
    lower.Set(MIN(a.lower.x(),b.lower.x()),
	      MIN(a.lower.y(),b.lower.y()));
    upper.Set(MAX(a.upper.x(),b.upper.x()),
	      MAX(a.upper.y(),b.upper.y()));
  }

  void setToUnion(Bboxf2D& b) {
    lower.Set(MIN(lower.x(),b.lower.x()),
	      MIN(lower.y(),b.lower.y()));
    upper.Set(MAX(upper.x(),b.upper.x()),
	      MAX(upper.y(),b.upper.y()));
  }

  void includePoint(float x, float y) {
    if (x < lower.x()) lower.set_x(x);
    if (y < lower.y()) lower.set_y(y);

    if (x>upper.x()) upper.set_x(x);
    if (y>upper.y()) upper.set_y(y);
  }

  void includePoint(const Vec2f& pt) {
    if (pt.x() < lower.x()) lower.set_x(pt.x());
    if (pt.y() < lower.y()) lower.set_y(pt.y());
    
    if (pt.x()>upper.x()) upper.set_x(pt.x());
    if (pt.y()>upper.y()) upper.set_y(pt.y());
  }

  void print() const {
    printf("from %f,%f to %f,%f\n",lower.x(),lower.y(),upper.x(),upper.y());
  }

  bool intersect(const Bboxf2D& other) const {
    if (upper.x() < other.lower.x() || other.upper.x() < lower.x() ||
	upper.y() < other.lower.y() || other.upper.y() < lower.y())
      return 0;
    return 1;
  }

  float width() const {return upper.x()-lower.x(); }
  float height() const {return upper.y()-lower.y(); }


  
  // ---- data
  Vec2f lower;
  Vec2f upper;

};

#endif

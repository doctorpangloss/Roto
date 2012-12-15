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



#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include "bBoxf2D.h"
#include "rotoPath.h"

enum TransformMode { M_ROTATE, M_UNIFORMSCALE, M_TRANSLATE };

class Transformer {
 public:

  Transformer();

  // adds segment to selected
  void addSplineSegment(RotoPath* rp, float t);

  // adds whole spline to selected
  void addSpline(RotoPath* rp);
  
  // true if there is something to transform
  bool transformable() const;

  // true if currently transforming something
  bool transforming() const;

  // used to determine transform mode (rotate, scale), starting config, start transform
  void startTransform(const Vec2f& loc);

  // during transformation, mouse moved to new loc
  void handleNewLoc(const Vec2f& loc);

  void endTransform(const Vec2f& loc);

  void clear();

  void render();

 private:

  void recalcBBox();

  bool _transforming;
  TransformMode _mode;
  std::set<RotoPath*> _selected;
  Bboxf2D _bbox;
  Vec2f _centroid;
  float _last;
  Vec2f _transformPoint;

};

#endif

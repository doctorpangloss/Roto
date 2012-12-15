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



#ifndef ROTOREGION_H
#define ROTOREGION_H



//enum RotoRegionSide { RR_LEFT, RR_RIGHT };

#define RR_LEFT 0
#define RR_RIGHT 1

#include "rotoPath.h"
class RotoPath;
typedef std::vector<RotoPath*> PathV;


class RotoRegion {

 public:

  RotoRegion(const PathV& v);
  RotoRegion();
  ~RotoRegion();

  int getFirst(const int i);
  int getLast(const int i);

  void createMask();
  void calculateMask();

  unsigned char* copyMask() const ;
  void orMask(unsigned char* omask) const;

  void render() const;

  void renderMask(uchar* bits) const;

  bool withinRegion(const float x, const float y);

  RotoRegion* tryToCopyForward();

  static int _w, _h;
 
  PathV _rotoPaths;
  std::vector<short> _sides;
  unsigned char *_mask;
  int _listNum;
};




#endif

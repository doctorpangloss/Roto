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



#ifndef __POINT_HH__
#define __POINT_HH__

#include <math.h>
#include <assert.h>

struct Point
{
  float x,y;

  Point() { };
  Point(float x,float y) { this->x = x; this->y = y; }
  Point(const Point & v) { x = v.x;y=v.y; }
  Point & operator=(const Point & v) { x = v.x;y=v.y; return *this; }

  Point operator+(const Point & v) const { return Point(x+v.x,y+v.y); }
  Point operator-(const Point & v) const { return Point(x-v.x,y-v.y); }
  Point operator*(float v) const { return Point(x*v,y*v); }
  Point operator/(float v) const { return Point(x/v,y/v); }
  bool operator==(const Point & v) const { return x == v.x && y == v.y; }
  Point operator-() const { return Point(-x,-y); }

  Point & operator+=(const Point & v) { x += v.x;y+=v.y; return *this; }
  float L2sqr() const { return x*x+y*y; }
  float L2() const { return sqrt(x*x+y*y); }
  Point round() const { return Point(floor(x+.5),floor(y+.5)); }
  float dot(const Point & p) const { return x*p.x+y*p.y; }

  void normalize() { float mag=L2(); if (mag != 0) {x/=mag;y/=mag;} }
};


struct PointR : public Point
{
  float r;

  PointR() { };
  PointR(float x,float y,float r) : Point(x,y) { this->r = r; }
  PointR(const PointR & v) : Point(v.x,v.y), r(v.r) {  }
  Point& operator=(const PointR & v) { x = v.x;y=v.y;r=v.r; return *this; }

  PointR operator+(const PointR & v) const { return PointR(x+v.x,y+v.y,r+v.r);}
  PointR operator-(const PointR & v) const { return PointR(x-v.x,y-v.y,r+v.r);}
  PointR operator*(float v) const { return PointR(x*v,y*v,r*v); }
  PointR operator/(float v) const { return PointR(x/v,y/v,r/v); }
  PointR operator-() const { return PointR(-x,-y,-r); }
  bool operator==(const PointR & v) const { return x == v.x && y == v.y &&
                                           r == v.r; }

  PointR & operator+=(const PointR & v) { x += v.x;y+=v.y;r+=v.r;return *this;}
  float L2sqr() const { return x*x+y*y; }
  float L2() const { return sqrt(x*x+y*y); }
};



#endif

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



#ifndef CUBICEVAL_H
#define CUBICEVAL_H

#include <stdio.h>


class CubicEval {

public:
  
  CubicEval(const float a, const float b, const float c, const float d) :
    _a(a), _b(b), _c(c), _d(d) {}

  CubicEval() {}

  void init(const float a, const float b, const float c, const float d) {
    _a = a; _b = b; _c = c; _d = d;
    //printf("%f %f %f %f\n",a,b,c,d);
  }

  void print() const {
    printf("cubic Eval: %f %f %f %f\n",_a,_b,_c,_d);
  }
  
  float f(const float t) const {
    if (t==0) return _d;
    if (t==1) return _a+_b+_c+_d;
    return (_a*t*t*t + _b*t*t + _c*t + _d); }

  float deriv_1(const float t) const {
    if (t==0) return _c;
    return (3.f *_a * t*t + 2.f *_b * t + _c);
  }

  float deriv_2(const float t) const {
    return (6.f*_a*t + 2.f *_b);
  }

  float deriv_3() const {
    return (6.f*_a); }

  float curvatureIntegral(const CubicEval& other) const {
    float ax = 6.*_a;
    float bx = 2.*_b;
    float ay = 6.*other._a;
    float by = 2.*other._b;
    return (ax*ax + ay*ay)/3.f +  ax*bx + ay*by  +  bx*bx + by*by;
  }
  
  int derivRoots(const int start, float* roots) const {
    if (_a==0) return 0; // linear
    float determ = 4.f*_b*_b - 12.f*_a*_c;
    if (determ < 0) return 0;
    if (determ==0) {
      roots[start] = -2.f*_b/(6.f*_a);
      return 1;
    }
    
    roots[start] = (-2.f*_b + sqrt(determ)) / (6.f*_a);
    roots[start+1] = (-2.f*_b - sqrt(determ)) / (6.f*_a);
    return 2;
  }


private:
  float _a,_b,_c,_d;


};


struct CubicEval2 {
  CubicEval _x, _y;

};




#endif

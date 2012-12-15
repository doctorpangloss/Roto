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



#ifndef CORRESPONDER_H
#define CORRESPONDER_H

#include "dynArray.h"
#include "jl_vectors.h"
#include "abstractPath.h"

//define UNUSED  __attribute__((__unused__))

// STUB: normalization!
// SPEED: do every other point


enum DIR { HOR, VER, DIA };

class C_Node {
 public:

  C_Node() { W[0]=W[1]=W[2] = -1;}

  double W[3];
  DIR BT[3];  
  Vec2i _last[3];
};

// virtual class!  Should subclass with cost functions,
// CI,W_s,W_b

class Corresponder {

 public:

  Corresponder(AbstractPath* X, AbstractPath* Y);
  virtual ~Corresponder() { delete[] _graph; }

  void getSolution(int* xints, int* yints) const;
  
  void calculate();

 protected:

  virtual double CI(const int i, const int j){ 
    return 0; //fabs(double(i)/double(_w) - double(j)/double(_h));
  }
    
  
  virtual double W_s(const int i1, const int j1, const int i2, const int j2) const {
    return 0; }
  
  virtual double W_b(const int i1, 
		     const int j1, 
		     const int i2, 
		     const int j2,
		     const int i3, 
		     const int j3) {
    return 0; }
  
  /* W_b MUST have beginning if clause here, or else!
     Not abstracting for efficiency reasons
  double _W_b(const int i1, const int j1, 
	      const int i2, const int j2,
	      const int i3, const int j3) {
    if (i1<0 || j1<0)
      return _myInf;
  }
  */

  AbstractPath *_X, *_Y;
  int _w,_h;
  double _myInf;
  C_Node* _graph;
};


#endif

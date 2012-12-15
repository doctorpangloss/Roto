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



#ifndef GEIGERCORRESPONDER_H
#define GEIGERCORRESPONDER_H

#include "dynArray.h"
#include "jl_vectors.h"
#include "abstractPath.h"
#include "gcnode.h"
#ifdef __APPLE__
#define _X _X
#define _Y _Y
#endif

//typedef DynArray<Vec2f,200> AbstractPath;

class GeigerCorresponder {

 public:

  // X should be driving axis
  GeigerCorresponder(AbstractPath* X, AbstractPath* Y);
  virtual ~GeigerCorresponder() { delete[] _graph; }

  void getSolution(int* xints) const;
  
  double calculate();

 protected:

  virtual double linkCost(const int i, const int j, 
		  const int ip, const int jp) const;

  virtual double internalCost(const int i, const int j) const;
		    
  AbstractPath *_X, *_Y;
  int _w,_h;
  GC_Node* _graph;
  double _myInf;
  int _maxJump;
};



#endif

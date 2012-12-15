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



#ifndef DRAWCORRESPONDER_H
#define DRAWCORRESPONDER_H

#ifdef __APPLE__
using namespace std;
#endif
#include "rotoCurves.h"
#include "gcnode.h"
#include <map>
#ifdef __APPLE__
#define _X _X
#endif

class DrawCorresponder {

 public:

  DrawCorresponder(const AbstractPath* X, RotoCurves* rotos, const int imw, const int imh);
  ~DrawCorresponder();

  void getSolution(int* xints, RotoPath** ptrs) const;
  
  double calculate();
  

 private:

  double internalCost(const int i, const int j) const;
  double linkCost(const int i, const int j, const int ip, const int jp) const;


  const AbstractPath* _X;
  const RotoCurves* _rotos;
  RotoPath** _rptrs;
  int _w, _h, _bw, _bh;
  GC_Node* _graph;
  map<RotoPath*, int>* _intmap;
  mutable double _myInf;
  DynArray<int, 4>* _buckets;
};

#endif

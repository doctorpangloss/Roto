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



#include "drawCorresponder.h"

#define GRAPH(a,b) _graph[(b)*_w+(a)]

#define GC_1 1.
#define GC_2 .2
//#define GC_3 1.

#define _BSIZE_ 12
#define _FBSIZE_ 12.

DrawCorresponder::DrawCorresponder(const AbstractPath* X, RotoCurves* rotos, const int imw, const int imh) :
  _X(X), _rotos(rotos) {
  
  int i, j;
  _h = 0;
  RotoPathList::const_iterator c;
  for (c=rotos->begin(); c!= rotos->end(); ++c) 
    _h += (*c) ->getNumElements();
  
  _rptrs = new RotoPath *[_h];
  _intmap = new map<RotoPath*, int>();
  
  _bw = imw/_BSIZE_ + 1;
  _bh = imh/_BSIZE_ + 1;
  _buckets = new DynArray<int, 4>[_bw*_bh];
  

  j=0;
  Vec2f loc;
  Vec2i buckLoc;
  for (c=rotos->begin(); c!= rotos->end(); ++c) {
    _intmap->insert(make_pair((*c),j));
    for (i=0; i<(*c)->getNumElements(); i++, j++) {
      loc = (*c)->getElement(i);
      buckLoc.Set((int)(loc.x()/_FBSIZE_), (int)(loc.y()/_FBSIZE_));
      buckLoc.Clamp(0,_bw-1,_bh-1);
      _buckets[buckLoc.y()*_bw + buckLoc.x()].add(j);
      _rptrs[j] = (*c);
    }
  }
  
  _w = _X->getNumElements();
  _graph = new GC_Node[_w*_h];
}

DrawCorresponder::~DrawCorresponder() {
  delete[] _rptrs;
  delete[] _graph;
  delete[] _buckets;
  delete _intmap;
}

void DrawCorresponder::getSolution(int* xints, RotoPath** ptrs) const {
  int i, bptr;
  //int w = _w, h = _h;

  double minCost = DBL_MAX;
  int iMinCost = -1;
  for (i=0; i<_h; i++)
    if (GRAPH(_w-1,i)._cost < minCost) {
      minCost = GRAPH(_w-1,i)._cost;
      iMinCost = i;
    }
  
  bptr = iMinCost; i = _w-1;

  while (i>-1) {
    assert(bptr>-1 && bptr<_h);
    xints[i] = bptr - (*_intmap)[_rptrs[bptr]];
    ptrs[i] = _rptrs[bptr];
    assert(xints[i] < ptrs[i]->getNumElements());
    bptr = GRAPH(i,bptr)._bptr;
    i--;
  }
}
  
double DrawCorresponder::calculate() {
  _myInf = DBL_MAX / 10.;  // to avoid overflow  
  int i,j,k, bi,bj,bp,buckindex;
  Vec2f jLoc;
  Vec2i buckLoc;

  double minCost, cost, internCost;
  int iminCost;
  for (i=0; i<_w; i++) 
    for (j=0; j<_h; j++) { 
      
      minCost = _myInf;
      iminCost = -1;
      
      internCost = internalCost(i,j);
      //offset = (_h-1)*_w;
      jLoc = _rptrs[j]->getElement(j - (*_intmap)[_rptrs[j]]);
      buckLoc.Set((int)(jLoc.x()/_FBSIZE_), (int)(jLoc.y()/_FBSIZE_));
      buckLoc.Clamp(0,_bw-1,_bh-1);
      // _buckets[buckLoc.y()*_bw + buckLoc.x()].add(j);

      //for (k = _h-1; k>-1; k--, offset-=_w) { // matching (i,j) and (i-1,k)

      for (bj=MAX(buckLoc.y()-2,0); bj <= MIN(_bh-1, buckLoc.y()+2); bj++)
	for (bi=MAX(buckLoc.x()-2,0); bi <= MIN(_bw-1, buckLoc.x()+2); bi++) {
	  buckindex = bj*_bw+bi;
	  for (bp=0; bp<_buckets[buckindex].getNumElements(); bp++) {
	    k = _buckets[buckindex].getElement(bp);   // matching (i,j) and (i-1,k)
	    //cost = _graph[offset+i-1]._cost + linkCost(i,j,i-1,k) + internCost;  
	    cost = GRAPH(i-1,k)._cost + linkCost(i,j,i-1,k) + internCost;  // SPEED: GRAPH sucks
	    if (cost < minCost) {
	      minCost = cost;
	      iminCost = k;
	    }
	  }
	}
      
      assert(iminCost != -1 && iminCost<_h);
      GRAPH(i,j)._cost = minCost;
      GRAPH(i,j)._bptr = iminCost;
    }

  minCost = DBL_MAX;
  for (j=0; j<_h; j++)
    if (GRAPH(_w-1,j)._cost < minCost)
      minCost = GRAPH(_w-1,j)._cost;

  assert(minCost < _myInf);
  return minCost;
  
}


double DrawCorresponder::linkCost(const int i, const int j, 
				  const int ip, const int jp) const {
  if (ip<0) return 0;
  double cost = 0;
  
  Vec2f D(_rptrs[j]->getElement(j - (*_intmap)[_rptrs[j]]), _X->getElement(i));
  D -= _rptrs[jp]->getElement(jp - (*_intmap)[_rptrs[jp]]);
  D += _X->getElement(ip);
  cost += /*GC_1 **/ D.Len();

  if (_rptrs[j] != _rptrs[jp]) cost += 20.;

  //if (j==jp)
  //cost += 10;

  //cost += GC_2 * (j-jp)*(j-jp);

  return cost;
}

double DrawCorresponder::internalCost(const int i, const int j) const {
  double cost; //,shit;
  //shit = _X->tangent(i).Dot2(_Y->tangent(j));
  //shit = MIN(1,MAX(-1,shit)); // since vecs are floats...
  //cost = GC_3 * acos(shit);
  
  cost = GC_2 * _X->getElement(i).distanceTo(_rptrs[j]->getElement(j - (*_intmap)[_rptrs[j]]));
  
  return cost;
}

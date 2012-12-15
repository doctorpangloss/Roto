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



#include <float.h>
#include <errno.h>
#include "geigerCorresponder.h"

#define GRAPH(a,b) _graph[(b)*_w+(a)]

#define GC_1 1.
#define GC_2 1.
#define GC_3 1.

GeigerCorresponder::GeigerCorresponder(AbstractPath* X, AbstractPath* Y) : _X(X), _Y(Y) {
  assert(X && Y);
  _w = _X->getNumElements();
  _h = _Y->getNumElements();
  _graph = new GC_Node[_w*_h];
  _maxJump = (int) ceil(3*double(_h)/double(_w));
  printf("max jump %d\n", _maxJump);
}

void GeigerCorresponder::getSolution(int* xints) const {
  int i = _w-1, bptr = _h-1;
  //int w = _w, h = _h;
  
  while (i>-1) {
    assert(bptr>-1);
    xints[i] = bptr;
    bptr = GRAPH(i,bptr)._bptr;
    i--;
  }
  printf("\n");

}

// SPEED: if we stay only using internalCost, we can speed this up:
// compile for each column a column that stores minimum i-1 cost for this
// j (and that index), use this.  n^3 -> n^2 algorithm.
double GeigerCorresponder::calculate() {
  _myInf = DBL_MAX / 10.;  // to avoid overflow  

  int i,j,k;
  _graph[0]._cost = 0;
  for (j=1; j<_h; j++)
    GRAPH(0,j)._cost = _myInf;

  double minCost, cost, internCost;
  int iminCost;
  for (i=1; i<_w; i++) 
    for (j=0; j<_h; j++) { 

      double dia = fabs(double(i)/double(_w-1) - double(j)/double(_h-1));
      if (dia > .3) {
	GRAPH(i,j)._cost = _myInf;
	continue;
      }

      minCost = DBL_MAX;
      iminCost = -1;

      internCost = internalCost(i,j);
      for (k = j; k>-1; k--) { // matching (i,j) and (i-1,k)
	cost = GRAPH(i-1,k)._cost + linkCost(i,j,i-1,k) + internCost;
	if (cost < minCost) {
	  minCost = cost;
	  iminCost = k;
	}	
	if (j-k > _maxJump)
	  break;
      }	

      assert(iminCost != -1);
      GRAPH(i,j)._cost = minCost;
      GRAPH(i,j)._bptr = iminCost;
  }
  assert(minCost < _myInf);
  return GRAPH(_w-1, _h-1)._cost;
}


double GeigerCorresponder::linkCost(const int i, const int j, 
		const int ip, const int jp) const {
  double cost = 0;

  Vec2f D(_Y->getElement(j), _X->getElement(i));
  D -= _Y->getElement(jp);
  D += _X->getElement(ip);
  cost += GC_1 * D.Len();
  if (j==jp)
    cost += 10;

  //cost += GC_2 * (j-jp)*(j-jp);

  return cost;
}

double GeigerCorresponder::internalCost(const int i, const int j) const {
  double cost; //,shit;
  //shit = _X->tangent(i).Dot2(_Y->tangent(j));
  //shit = MIN(1,MAX(-1,shit)); // since vecs are floats...
  //cost = GC_3 * acos(shit);

  cost = _X->getElement(i).distanceTo(_Y->getElement(j));

  return cost;
}

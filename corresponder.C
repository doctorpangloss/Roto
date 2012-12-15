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
#include <stdio.h>
#include "corresponder.h"

#define WVER(a,b) _graph[(b)*_w+(a)].W[VER]
#define WHOR(a,b) _graph[(b)*_w+(a)].W[HOR]
#define WDIA(a,b) _graph[(b)*_w+(a)].W[DIA]
#define GRAPH(a,b) _graph[(b)*_w+(a)]
#define MIN3(a,b,c) (MIN(MIN((a),(b)),(c)))
#define Kb .5

// taken from Sederberg & Greenwood

Corresponder::Corresponder(AbstractPath* X, AbstractPath* Y) : _X(X), _Y(Y) {
  _w = X->getNumElements(); _h = Y->getNumElements();
  _myInf = DBL_MAX / 10.;  // to avoid overflow
  _graph = new C_Node[_w*_h];
}

void Corresponder::calculate() {

  int i,j;
  WVER(0,0) = 0;
  WHOR(0,0) = 0;
  WDIA(0,0) = 0;
  double ci, tmph, tmpv, tmpd;

  for (i=0; i<_w; i++) {
    for (j=0; j<_h; j++) {
      if (i==0 && j==0) continue;
      
      /*
      double dia = fabs(double(i)/double(_w) - double(j)/double(_h));
      if (dia > .3) {
	WHOR(i,j) = WVER(i,j) = WDIA(i,j) = _myInf;
	continue;
      }
      */

      ci = CI(i,j);

      // horizontal
      if (i!=0) {
	tmph = WHOR(i-1,j) + Kb*W_b(i-2,j,i-1,j,i,j);
	tmpd = WDIA(i-1,j) + Kb*W_b(i-2,j-1,i-1,j,i,j);
	if (tmph < tmpd) {
	  GRAPH(i,j)._last[HOR].Set(i-1, GRAPH(i-1,j)._last[HOR].y());
	  WHOR(i,j) = ci + W_s(i-1,j,i,j) + tmph;
	  GRAPH(i-1,j).BT[HOR] = HOR;
	}
	else {
	  GRAPH(i,j)._last[HOR].Set(i-1, GRAPH(i-1,j)._last[DIA].y());
	  WHOR(i,j) = ci + W_s(i-1,j,i,j) + tmpd;
	  GRAPH(i-1,j).BT[HOR] = DIA;
	}
      }
      else 
	WHOR(0,j) = _myInf;
      assert(WHOR(i,j) >= 0);

      // vertical
      if (j!=0) {
	tmpv = WVER(i,j-1) + Kb*W_b(i,j-2,i,j-1,i,j);
	tmpd = WDIA(i,j-1) + Kb*W_b(i-1,j-2,i,j-1,i,j);
	if (tmpv < tmpd) {
	  GRAPH(i,j)._last[VER].Set(GRAPH(i,j-1)._last[VER].x(), j-1);
	  WVER(i,j) = ci + W_s(i,j-1,i,j) + tmpv;
	  GRAPH(i,j-1).BT[VER] = VER;
	}
	else {
	  GRAPH(i,j)._last[VER].Set(GRAPH(i,j-1)._last[DIA].x(), j-1);
	  WVER(i,j) = ci + W_s(i,j-1,i,j) + tmpd;
	  GRAPH(i,j-1).BT[VER] = DIA;
	}
      }
      else 
	WVER(i,0) = _myInf;
      assert(WVER(i,j) >= 0);

      // diagonal
      if (i!=0 && j!=0) {
	tmph = WHOR(i-1,j-1) + Kb*W_b(i-2,j-1,i-1,j-1,i,j);
	tmpv = WVER(i-1,j-1) + Kb*W_b(i-1,j-2,i-1,j-1,i,j);
	tmpd = WDIA(i-1,j-1) + Kb*W_b(i-2,j-2,i-1,j-1,i,j);
	double tmp = MIN3(tmph,tmpv,tmpd);
	GRAPH(i,j)._last[DIA].Set(i-1,j-1);
	if (tmp == tmph) {
	  WDIA(i,j) = ci + W_s(i-1,j-1,i,j) + tmph;
	  GRAPH(i-1,j-1).BT[DIA] = HOR;
	}
	else if (tmp == tmpv) {
	  WDIA(i,j) = ci + W_s(i-1,j-1,i,j) + tmpv;
	  GRAPH(i-1,j-1).BT[DIA] = VER;
	}
	else {
	  assert(tmp==tmpd);
	  WDIA(i,j) = ci + W_s(i-1,j-1,i,j) + tmpd;
	  GRAPH(i-1,j-1).BT[DIA] = DIA;
	}
      }
      else
	WDIA(i,j) = _myInf;
      assert(WDIA(i,j) >= 0);

    }
  }  
}

void Corresponder::getSolution(int* xints, int* yints) const {
  // now backtrack!
  int i = _w-1, j = _h-1;
  DIR currDir;
  double best = MIN3(WHOR(i,j),WVER(i,j),WDIA(i,j));
  printf("min work %f\n",best);
  if (WHOR(i,j) == best)
    currDir = HOR;
  else if (WVER(i,j) == best)
    currDir = VER;
  else
    currDir = DIA;

  while (i>=0 && j>=0) {
    assert(i>-1 && j>-1);
    if (xints) xints[i] = j;
    if (yints) yints[j] = i;
    
    switch (currDir) {
    case HOR:
      i--; break;
    case VER:
      j--; break;
    case DIA:
      i--; j--; break;
    default:
      assert(0);
    }  

    currDir = GRAPH(i,j).BT[currDir];    
  } // end while

  
  
  /* // DEBUG
  int k;
  if (xints)
    for (k=0; k<_w; k++)
      assert(xints[k] >=0 && xints[k] < _h);
  if (yints)
    for (k=0; k<_h; k++)
      assert(yints[k] >=0 && yints[k] < _w);
  */
  /*
  if (_w ==100) {
    int k,x=0,y=0;
    double shit = 0;
    for (k=0; k<99; k++) {
      shit += W_s(x,y,x+1,y+1);
      x++; y++;
      shit += W_s(x,y,x,y+1);
      y++;
    }
    //assert(x==100 & y==199);
    printf("yippee %f\n",shit);
  }
  */
}




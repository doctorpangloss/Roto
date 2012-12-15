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



#include "rotoCorresponder.h"

double RotoCorresponder::linkCost(const int i, const int j, 
				  const int ip, const int jp) const {
  double cost = 0;
  
  Vec2f D(_Y->getElement(j), _X->getElement(i));
  D -= _Y->getElement(jp);
  D += _X->getElement(ip);
  cost += 1. * D.Len2();
  if (j==jp)
    cost += 3.;
  
  cost += 1. * (j-jp)*(j-jp);
  
  return cost;
}

double RotoCorresponder::internalCost(const int i, const int j) const {
  //printf("right on\n");
  double cost,shit;
  shit = _X->tangent(i).Dot2(_Y->tangent(j));
  shit = MIN(1,MAX(-1,shit)); // since vecs are floats...
  cost = 1. * acos(shit);
  return cost;
}

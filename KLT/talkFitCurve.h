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



#ifndef TALK_H
#define TALK_H

#include <assert.h>
#include <math.h>
#include "llist.h"
#include "jl_vectors.h"


#define T_STEPS 20


class Bez {
public:
  Bez() {lastFrozen=0;}  
  Vec2f _data[4];
  Bez* next;
  Bez* prev;
  bool lastFrozen;
};

class TalkFitCurve : public LList<Bez> {

public:

  TalkFitCurve() : LList<Bez>() {}

  void addBez(float x1,float y1,
	      float x2,float y2,
	      float x3,float y3,
	      float x4,float y4) {

    
    Bez* newB = new Bez();
    assert(finite(x1) && finite(y1) && finite(x2) && finite(y2) &&
	   finite(x3) && finite(y3) && finite(x4) && finite(y4));
    newB->_data[0].Set(x1,y1);
    newB->_data[1].Set(x2,y2);
    newB->_data[2].Set(x3,y3);
    newB->_data[3].Set(x4,y4);
    AddToTail(newB);
  }

};

#endif

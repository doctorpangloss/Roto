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



#include <math.h>
#include "corresponder.h"


class WireCorresponder : public Corresponder {
  
 public:

  WireCorresponder(AbstractPath* X, AbstractPath* Y) : 
    Corresponder(X,Y) {}
  ~WireCorresponder() {}

 protected:

  double W_s(const int i0, const int j0, const int i1, const int j1) const {
    double L0,L1,res;
    L0 = _X->getPointerToElement(i0)->distanceTo(*(_X->getPointerToElement(i1)));
    L1 = _Y->getPointerToElement(j0)->distanceTo(*(_Y->getPointerToElement(j1)));
    res =  2.*(L1-L0)*(L1-L0) / (L0 + L1);
    assert (res >=0 && res < 100000);
    return res;
  }
  
  
  double W_b(const int i1, const int j1, 
	     const int i2, const int j2,
	     const int i3, const int j3) {
    if (i1<0 || j1<0) {
      if (i1==-1 && j1==-1)
	return 0;
      return _myInf;
    }
    Vec2f P1,P2;
    double a1,a2,da;

    // i
    if (i1==i2 && i2==i3)
      a1 = M_PI;
   
    else if (i2==i3) {
      Vec2f_Sub(P1,*(_X->getPointerToElement(i1)),
		*(_X->getPointerToElement(i2)));
    
      if (i3<_w-1) {
	Vec2f_Sub(P2,*(_X->getPointerToElement(i2)),
		  *(_X->getPointerToElement(i3+1)));
	a1 = atan2(P1.Cross2(P2),P1.Dot2(P2)) + M_PI;
	a1 = M_PI_2 + a1/2;
      }
      else 
	a1 = M_PI;
    }
    
    else if (i1==i2) {
      Vec2f_Sub(P2,*(_X->getPointerToElement(i2)),
		*(_X->getPointerToElement(i3)));
      if (i1>0) {
	Vec2f_Sub(P1,*(_X->getPointerToElement(i1-1)),
		  *(_X->getPointerToElement(i2)));
	a1 = atan2(P1.Cross2(P2),P1.Dot2(P2)) + M_PI;
	a1 = M_PI_2 + a1/2;
      }
      else
	a1 = M_PI;
    }
    
    else {
      Vec2f_Sub(P1,*(_X->getPointerToElement(i1)),
		*(_X->getPointerToElement(i2)));
      Vec2f_Sub(P2,*(_X->getPointerToElement(i2)),
		*(_X->getPointerToElement(i3)));  
      a1 = atan2(P1.Cross2(P2),P1.Dot2(P2)) + M_PI;
    }


    // j
    if (j1==j2 && j2==j3)
      a2 = M_PI;
    
    else if (j2==j3) {
      Vec2f_Sub(P1,*(_Y->getPointerToElement(j1)),
		*(_Y->getPointerToElement(j2)));
      
      if (j3<_h-1) {
	Vec2f_Sub(P2,*(_Y->getPointerToElement(j2)),
		  *(_Y->getPointerToElement(j3+1)));
	a2 = atan2(P1.Cross2(P2),P1.Dot2(P2)) + M_PI;
	a2 = M_PI_2 + a2/2;
      }
      else 
	a2 = M_PI;
    }
    
    else if (j1==j2) {
      Vec2f_Sub(P2,*(_Y->getPointerToElement(j2)),
		*(_Y->getPointerToElement(j3)));
      if (j1>0) {
	Vec2f_Sub(P1,*(_Y->getPointerToElement(j1-1)),
		  *(_Y->getPointerToElement(j2)));
	a2 = atan2(P1.Cross2(P2),P1.Dot2(P2)) + M_PI;
	a2 = M_PI_2 + a2/2;
      }
      else
	a2 = M_PI;
    }
    
    else {
      Vec2f_Sub(P1,*(_Y->getPointerToElement(j1)),
		*(_Y->getPointerToElement(j2)));
      Vec2f_Sub(P2,*(_Y->getPointerToElement(j2)),
		*(_Y->getPointerToElement(j3)));  
      a2 = atan2(P1.Cross2(P2),P1.Dot2(P2)) + M_PI;
    }
    
    /*
    Vec2f_Sub(P1,*(_Y->getPointerToElement(j1)),
	      *(_Y->getPointerToElement(j2)));
    Vec2f_Sub(P2,*(_Y->getPointerToElement(j2)),
	      *(_Y->getPointerToElement(j3)));
    a2 = atan2(P1.Cross2(P2),P1.Dot2(P2)) + M_PI;
    */

    da = fabs(a2-a1);
    da = MIN(2*M_PI-da,da);
    assert(da>=0 && da < 100000);

    return da;
  }
  

};

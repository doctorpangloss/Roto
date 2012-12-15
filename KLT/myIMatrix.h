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



#ifndef MYIMATRIX_H
#define MYIMATRIX_H

#include <qimage.h>
#include "linearSolver.h"

class MyIMatrix : public implicitMatrix {

 public:

  MyIMatrix(const double *Z, const int numFrames, const int l2) : _Z(Z), _numFrames(numFrames), _l2(l2) {
    _d=_l2*_numFrames;
  }


  int numVar() const { return _numFrames*_l2; }
  
  virtual void matVecMult(const double x[], double r[]) const{
    //printf("\n");
    int i,offset=0, col, t, tl2, n,low,high;

    for (t=0,i=0,tl2=0; t<_numFrames; t++, tl2+=_l2)
      for (n=0; n<_l2; n++, i++,offset+=_d) {
	r[i]=0;
	
	// main diagonal
	low=MAX(tl2,i-5); high=MIN(tl2+_l2,i+6);
	for (col=low; col<high; col++) {
	  r[i] += _Z[offset+col]*x[col];
	}

	// left strip
	if (tl2>0) {
	  low=MAX(tl2-_l2,i-_l2-5); high=MIN(tl2,i-_l2+6);
	  for (col=low; col<high; col++) {
	    r[i] += _Z[offset+col]*x[col];
	  }
	}

	// right strip
	low=MAX(tl2+_l2, i+_l2-5); high=MIN(tl2+_l2+_l2, i+_l2+6);
	if (low < _d) {
	  for (col=low; col<high; col++) {
	    r[i] += _Z[offset+col]*x[col];
	    //printf("%7.2f ",_Z[offset+col]);
	  }
	}
	//printf("\n");
      } // end n
    
    /*    
    int i,j, row=0, col;
    for (i=0; i<_d; i++, row+=_d) {
      r[i]=0;

      // main diagonal
      for (j=-5; j<=5; j++) {
	col = i + j;
	if (col>-1 && col <_d) {
	  r[i] += _Z[row+col]*x[col];
	}
      }

      // left strip
      if (i-_l2 >= 0) {
	for (j=-5; j<=5; j++) {
	  col = i-_l2 + j;
	  if (col>-1) {
	    r[i] += _Z[row+col]*x[col];
	  }
	}
      }

      // right strip
      if (i+_l2 < _d) {
	for (j=-5; j<=5; j++) {
	  col = i+_l2 + j;
	  if (col <_d) {
	    r[i] += _Z[row+col]*x[col];
	  }
	}
      }
	
    }
    */

    //for (t=0; t<_numFrames*_l2; t++)
    //printf("%7.2f ",r[t]);
    //printf("\n\n");
  }

 private:
  
  const double *_Z;
  int _d;
  const int _numFrames,_l2;
};


#endif

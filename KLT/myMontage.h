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



#ifndef MYMONTAGE_H
#define MYMONTAGE_H

#include <qimage.h>

class MyMontage {
  
 public:

  MyMontage(QImage *ims, int num, int mag, int numCols) {
    int tw = ims[0].width(), th = ims[0].height();
    int numRows = num/numCols+1, xoffset,yoffset;
    _mont = QImage(numCols*tw*mag, numRows*th*mag, QImage::Format_ARGB32);
    int whichImage=0;
    for (int j=0; j<numRows; j++) {
      yoffset = j*th*mag;
      for (int i=0; i<numCols; i++) {
	xoffset = i*tw*mag;
	
	int magw=ims[whichImage].width()*mag, magh = ims[whichImage].height()*mag;
	for (int imj = 0; imj<magh; imj++)
	  for (int imi = 0; imi<magw; imi++) {
	    assert(imi/mag < ims[whichImage].width() &&
		   imj/mag < ims[whichImage].height());
	    _mont.setPixel(xoffset+imi, yoffset+imj,
			   ims[whichImage].pixel(imi/mag,imj/mag));
	  }

	whichImage++;
	if (whichImage >= num) return;
      }
    }

  }

  void save(char* filename) {
    _mont.save(filename,"PNG");
  }

 private:
  QImage _mont;

};

#endif

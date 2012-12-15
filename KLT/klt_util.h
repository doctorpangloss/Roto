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



#ifndef _KLT_UTIL_H_
#define _KLT_UTIL_H_

class QImage;
#include <qimage.h>
#include "kernels.h"

class  KLT_FloatImage {
 public:

  KLT_FloatImage(int w, int h);
  KLT_FloatImage();
  KLT_FloatImage(const QImage im); // currently assumes this is a color image to convert to grayscale

  // this constructor takes color image derivatives, combines into one
  KLT_FloatImage(int w, int h, 
		 const KLT_FloatImage* drx, const KLT_FloatImage* dry,
		 const KLT_FloatImage* dgx, const KLT_FloatImage* dgy,
		 const KLT_FloatImage* dbx, const KLT_FloatImage* dby);

  ~KLT_FloatImage() { delete[] data; }  

  KLT_FloatImage* getSmoothed(const Kernels* kern);
  
  void setSize(int w, int h); 
  
  void takeRed(const QImage im);
  void takeBlue(const QImage im);
  void takeGreen(const QImage im);
  
  void computeGradients(Kernels* kern, KLT_FloatImage* gradx, KLT_FloatImage* grady) const;

  float interpolate(const float x, const float y, int* status) const;

  float& g(unsigned int i) { return data[i]; }
  float v(unsigned int i) const { return data[i]; }

  int ncols;
  int nrows;
  float *data;

 private: 
  void convolveImageVert(const ConvolutionKernel* kernel,
			 KLT_FloatImage* imgout) const;
  
  void convolveImageHoriz(const ConvolutionKernel* kernel,
			  KLT_FloatImage* imgout) const;
  
  void convolveSeparate(const ConvolutionKernel* horiz_kernel,
			const ConvolutionKernel* vert_kernel,
			KLT_FloatImage* imgout) const;
};



#endif



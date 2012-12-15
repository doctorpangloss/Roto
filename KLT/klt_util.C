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



/*********************************************************************
 * klt_util.c
 *********************************************************************/

/* Standard includes */
#include <assert.h>
#include <stdlib.h>  /* malloc() */

/* Our includes */
#include "base.h"
#include "error.h"
//#include "pnmio.h"
//#include "klt.h"
#include "klt_util.h"


/*
float _KLTComputeSmoothSigma(  // add to tracking context class
  KLT_TrackingContext tc)
{
  return (tc->smooth_sigma_fact * max(tc->window_width, tc->window_height));
}
*/



KLT_FloatImage::KLT_FloatImage(int w, int h) : ncols(w), nrows(h) {
  data = new float[ncols*nrows];
}

KLT_FloatImage::KLT_FloatImage() {
  data = NULL; ncols = -1; nrows = -1;
}

void KLT_FloatImage::setSize(int w, int h) {
  assert(data==NULL);
  ncols = w; nrows = h;
  data = new float[ncols*nrows];
}


KLT_FloatImage::KLT_FloatImage(const QImage im) {
  ncols = im.width(); nrows = im.height();
  data = new float[ncols*nrows];
  for (int j=0; j<nrows; j++)
    for (int i=0; i<ncols; i++) {
      QRgb color = im.pixel(i,j);
      data[j*ncols+i] = (.3*float(qRed(color)) + 
			 .59*float(qGreen(color))
			 +.11*float(qBlue(color))); // varies up to 255      
    }
}
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define _ME_ 70.0f 
#define fsqrt(a) (float) sqrt(a)
KLT_FloatImage::KLT_FloatImage(int w, int h, 
	       const KLT_FloatImage* drx, const KLT_FloatImage* dry,
	       const KLT_FloatImage* dgx, const KLT_FloatImage* dgy,
	       const KLT_FloatImage* dbx, const KLT_FloatImage* dby) {

  int index=0;
  float min = 100, tmp;
  ncols = w; nrows = h;
  data = new float[ncols*nrows];
  for (int j=0; j<nrows; ++j)
    for (int i=0; i<ncols; ++i, ++index) {
      data[index] = _ME_ - MAX(-_ME_, 
			       tmp = (  .3f* sqrt( drx->v(index)*drx->v(index) + dry ->v(index)*dry->v(index) ) + 
					.59f*sqrt( dgx->v(index)*dgx->v(index) + dgy ->v(index)*dgy->v(index) ) + 
					.11f*sqrt( dbx->v(index)*dbx->v(index) + dby ->v(index)*dby->v(index) )     )   );
      if (data[index] < min)
	min = data[index];
    }
  printf("min edge is %f\n",min);
}


// SPEED: could do these all in one loop?
void KLT_FloatImage::takeRed(const QImage im) {
  assert(im.width() == ncols && im.height() == nrows);
  for (int j=0; j<nrows; j++)
    for (int i=0; i<ncols; i++) {
      QRgb color = im.pixel(i,j);
      data[j*ncols+i] = float(qRed(color)); // varies up to 255      
    }
}

void KLT_FloatImage::takeBlue(const QImage im) {
  assert(im.width() == ncols && im.height() == nrows);
  for (int j=0; j<nrows; j++)
    for (int i=0; i<ncols; i++) {
      QRgb color = im.pixel(i,j);
      data[j*ncols+i] = float(qBlue(color)); // varies up to 255      
    }
}
void KLT_FloatImage::takeGreen(const QImage im) {
  assert(im.width() == ncols && im.height() == nrows);
  for (int j=0; j<nrows; j++)
    for (int i=0; i<ncols; i++) {
      QRgb color = im.pixel(i,j);
      data[j*ncols+i] = float(qGreen(color)); // varies up to 255      
    }
}

/*
 

void _KLTPrintSubFloatImage(
  _KLT_FloatImage floatimg,
  int x0, int y0,
  int width, int height)
{
  int ncols = floatimg->ncols;
  int offset;
  int i, j;

  assert(x0 >= 0);
  assert(y0 >= 0);
  assert(x0 + width <= ncols);
  assert(y0 + height <= floatimg->nrows);

  fprintf(stderr, "\n");
  for (j = 0 ; j < height ; j++)  {
    for (i = 0 ; i < width ; i++)  {
      offset = (j+y0)*ncols + (i+x0);
      fprintf(stderr, "%6.2f ", *(floatimg->data + offset));
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}
*/


/*
 

void _KLTWriteFloatImageToPGM(
  _KLT_FloatImage img,
  char *filename)
{
  int npixs = img->ncols * img->nrows;
  float mmax = -999999.9, mmin = 999999.9;
  float fact;
  float *ptr;
  uchar *byteimg, *ptrout;
  int i;

  // Calculate minimum and maximum values of float image 
  ptr = img->data;
  for (i = 0 ; i < npixs ; i++)  {
    mmax = max(mmax, *ptr);
    mmin = min(mmin, *ptr);
    ptr++;
  }
	
  // Allocate memory to hold converted image 
  byteimg = (uchar *) malloc(npixs * sizeof(uchar));

  // Convert image from float to uchar 
  fact = 255.0 / (mmax-mmin);
  ptr = img->data;
  ptrout = byteimg;
  for (i = 0 ; i < npixs ; i++)  {
    *ptrout++ = (uchar) ((*ptr++ - mmin) * fact);
  }

  // Write uchar image to PGM 
  pgmWriteFile(filename, byteimg, img->ncols, img->nrows);

  // Free memory 
  free(byteimg);
}

*/


KLT_FloatImage* KLT_FloatImage::getSmoothed(const Kernels* kern) {
   
  KLT_FloatImage *output = new KLT_FloatImage(ncols, nrows);
  
  
  /* Compute kernel, if necessary; gauss_deriv is not used */
  //if (fabs(sigma - sigma_last) > 0.05)
  //_computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
  
  convolveSeparate(kern->gauss(), kern->gauss(), output);
  return output;
}


void KLT_FloatImage::convolveSeparate(const ConvolutionKernel* horiz_kernel,
				      const ConvolutionKernel* vert_kernel,
				      KLT_FloatImage* imgout) const {
  assert(imgout->ncols == ncols && imgout->nrows == nrows);

  // Create temporary image 
  KLT_FloatImage* tmpimg = new KLT_FloatImage(ncols, nrows);
  
  // Do convolution 
  convolveImageHoriz(horiz_kernel, tmpimg);

  tmpimg->convolveImageVert(vert_kernel, imgout);

  // Free memory 
  delete tmpimg;
}



void KLT_FloatImage::convolveImageHoriz(const ConvolutionKernel* kernel,
					KLT_FloatImage* imgout) const {

  float *ptrrow = data;           /* Points to row's first pixel */
  float *ptrout = imgout->data, /* Points to next output pixel */
    *ppp;
  float sum;
  int radius = kernel->width / 2;
  int i, j, k;

  /* Kernel width must be odd */
  assert(kernel->width % 2 == 1);

  /* Must read from and write to different images */
  assert(this != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= ncols);
  assert(imgout->nrows >= nrows);

  /* For each row, do ... */
  for (j = 0 ; j < nrows ; j++)  {

    /* Zero leftmost columns */
    for (i = 0 ; i < radius ; i++)
      *ptrout++ = 0.0;

    /* Convolve middle columns with kernel */
    for ( ; i < ncols - radius ; i++)  {
      ppp = ptrrow + i - radius;
      sum = 0.0;
      for (k = kernel->width-1 ; k >= 0 ; k--)
        sum += *ppp++ * kernel->data[k];
      *ptrout++ = sum;
    }

    /* Zero rightmost columns */
    for ( ; i < ncols ; i++)
      *ptrout++ = 0.0;

    ptrrow += ncols;
  }
}




void KLT_FloatImage::convolveImageVert(const ConvolutionKernel* kernel,
				       KLT_FloatImage* imgout) const {
  float *ptrcol = data;          /* Points to row's first pixel */
  float *ptrout = imgout->data,  /* Points to next output pixel */
    *ppp;
  float sum;
  int radius = kernel->width / 2;
  int i, j, k;

  /* Kernel width must be odd */
  assert(kernel->width % 2 == 1);

  /* Must read from and write to different images */
  assert(this != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= ncols);
  assert(imgout->nrows >= nrows);

  /* For each column, do ... */
  for (i = 0 ; i < ncols ; i++)  {

    /* Zero topmost rows */
    for (j = 0 ; j < radius ; j++)  {
      *ptrout = 0.0;
      ptrout += ncols;
    }

    /* Convolve middle rows with kernel */
    for ( ; j < nrows - radius ; j++)  {
      ppp = ptrcol + ncols * (j - radius);
      sum = 0.0;
      for (k = kernel->width-1 ; k >= 0 ; k--)  {
        sum += *ppp * kernel->data[k];
        ppp += ncols;
      }
      *ptrout = sum;
      ptrout += ncols;
    }

    /* Zero bottommost rows */
    for ( ; j < nrows ; j++)  {
      *ptrout = 0.0;
      ptrout += ncols;
    }

    ptrcol++;
    ptrout -= nrows * ncols - 1;
  }
}


void KLT_FloatImage::computeGradients(Kernels* kern, KLT_FloatImage* gradx, 
				      KLT_FloatImage* grady) const {
  assert(gradx->ncols >= ncols);
  assert(gradx->nrows >= nrows);
  assert(grady->ncols >= ncols);
  assert(grady->nrows >= nrows);

  convolveSeparate(kern->gaussDeriv(), kern->gauss(), gradx);
  convolveSeparate(kern->gauss(), kern->gaussDeriv(), grady);
}

/*********************************************************************
 * _interpolate
 * 
 * Given a point (x,y) in an image, computes the bilinear interpolated 
 * gray-level value of the point in the image.  
 */

float KLT_FloatImage::interpolate(const float x, const float y, int* status) const {
  int xt = (int) x;  /* coordinates of top-left corner */
  int yt = (int) y;
  float ax = x - xt;
  float ay = y - yt;
  float *ptr = data + (ncols*yt) + xt;
  *status=1;
  
  if (xt<1 || yt<1 || xt>=ncols-4 || yt>=nrows-4) {  // This is set for handling gradients, may break older code 
      *status = 0;
      //fprintf(stderr, "(xt,yt)=(%d,%d)  imgsize=(%d,%d)\n"
      //      "(x,y)=(%f,%f)  (ax,ay)=(%f,%f)\n",
      //      xt, yt, ncols, nrows, x, y, ax, ay);
      //fflush(stderr); 
      return 0;
  }  

  float d1,d2=0;
  if (ay!=0) {
    d1 = *ptr + ay * (*(ptr+ncols) - *ptr);
    if (ax!=0) d2 = *(ptr+1) + ay * (*(ptr+ncols+1) - *(ptr+1));
  }
  else {
    d1 = *ptr;
    if (ax!=0) d2 = *(ptr+1);
  }
  
  if (ax!=0) {
    return (d1 + ax*(d2-d1));
  }
  else {
    return d1;
  }


  /*
  assert (xt >= 0 && yt >= 0 && xt <= ncols - 2 && yt <= nrows - 2);

  return ( (1-ax) * (1-ay) * *ptr +
           ax   * (1-ay) * *(ptr+1) +
           (1-ax) *   ay   * *(ptr+(ncols)) +
           ax   *   ay   * *(ptr+(ncols)+1) );
  */
}


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
 * pyramid.c
 *
 *********************************************************************/

/* Standard includes */
#include <assert.h>
#include <stdlib.h>		/* malloc() ? */
#include <string.h>		/* memset() ? */
#include <math.h>		/* */
#include <stdio.h>
#include <limits.h>

/* Our includes */
#include "base.h"
#include "error.h"
#include "pyramid.h"

#define max(a,b)	((a) > (b) ? (a) : (b))
#define min(a,b)	((a) < (b) ? (a) : (b))
/*
KLT_Pyramid::KLT_Pyramid() {
  img = NULL; ncols = NULL; nrows = NULL;
  _inited = false;
}

void KLT_Pyramid::setParams(int basecols, int baserows, int nsubsampling, int nnlevel) {
  _basencols = basecols;
  _basenrows = baserows;
  subsampling = nsubsampling;
  nLevels = nnlevel;
}
*/


KLT_Pyramid::KLT_Pyramid(int basencols,int basenrows,int nsubsampling,int nnlevels, bool initMe) :
  subsampling(nsubsampling), nLevels(nnlevels), _basencols(basencols), _basenrows(basenrows) {

  if (subsampling != 2 && subsampling != 4 && 
      subsampling != 8 && subsampling != 16 && subsampling != 32)
    KLTError("(_KLTCreatePyramid)  Pyramid's subsampling must "
             "be either 2, 4, 8, 16, or 32");

  if (initMe)
    setupStructure(basencols, basenrows);
  else {
    img = NULL; ncols = NULL; nrows = NULL;
    _inited = false;
  }
}

void KLT_Pyramid::initMe() {
  setupStructure(_basencols, _basenrows);
}

void KLT_Pyramid::setupStructure(int basencols, int basenrows) {
  img = new KLT_FloatImage[nLevels];
  ncols = new int[nLevels];
  nrows = new int[nLevels];

  for (int i = 0 ; i < nLevels ; i++)  {
    img[i].setSize(basencols, basenrows);
    ncols[i] = basencols; nrows[i] = basenrows;
    basencols /= subsampling; basenrows /= subsampling;
  }
  
  _inited = true;
}

void KLT_Pyramid::dumpPyramid() {
  delete[] img; img = NULL;
  delete[] ncols; ncols = NULL;
  delete[] nrows; nrows = NULL;
  _inited = false;
}
/*********************************************************************
 *
 */

KLT_Pyramid::~KLT_Pyramid() {
  if (_inited) {
    delete[] img;
    delete[] ncols;
    delete[] nrows;
  }
}


/*********************************************************************
 *
 */

void KLT_Pyramid::computePyramid(KLT_FloatImage* imgIn, float sigma_fact) {
  assert(_inited);
  KLT_FloatImage *currimg;
  int inncols = imgIn->ncols, innrows = imgIn->nrows;
  int subhalf = subsampling / 2;
  float sigma = subsampling * sigma_fact;  /* empirically determined */
  Kernels kern(sigma);
  int oldncols;
  int i, x, y;
	

  assert(ncols[0] == inncols);
  assert(nrows[0] == innrows);

  /* Copy original image to level 0 of pyramid */
  memcpy(img[0].data, imgIn->data, inncols*innrows*sizeof(float));

  currimg = imgIn;
  
  for (i = 1 ; i < nLevels ; i++)  {
    KLT_FloatImage* tmpimg = currimg->getSmoothed(&kern);

    /* Subsample */
    oldncols = inncols;
    inncols /= subsampling;  innrows /= subsampling;
    for (y = 0 ; y < innrows ; ++y)
      for (x = 0 ; x < inncols ; ++x)
        img[i].data[y*inncols+x] = 
          tmpimg->data[(subsampling*y+subhalf)*oldncols +
                      (subsampling*x+subhalf)];

    /* Reassign current image */
    currimg = img+i;
    delete tmpimg;
  }
}


/*********************************************************************
 *
 */

void KLT_Pyramid::write(FILE* fp) {
  int i;//, npixs, k;
  //int basecols = ncols[0];
  //int baserows = nrows[0];
  //float *ptr;
  //float mmax, mmin, fact;
  //ushort *byteimg, *ptrout;
  
  assert(fp);

  fwrite(ncols,sizeof(int),1,fp);
  fwrite(nrows,sizeof(int),1,fp);
  fwrite(&(subsampling),sizeof(int),1,fp);
  fwrite(&(nLevels),sizeof(int),1,fp);

  for (i = 0 ; i < nLevels ; i++)
    fwrite(img[i].data, sizeof(float), ncols[i]*nrows[i], fp); 

  /*
  for (i = 0 ; i < nLevels ; i++) {
    //printf("%d %d\n",baserows,basecols);

      npixs = basecols * baserows;
      mmax = -999999.9; mmin = 999999.9;

      // Calculate minimum and maximum values of float image 
      ptr = img[i].data;
      for (k = 0 ; k < npixs ; k++)  {
	  mmax = max(mmax, *ptr);
	  mmin = min(mmin, *ptr);
	  ptr++;
      }

      // write max,min data so values can be reconstructed 
      assert(mmax != mmin);
      //printf("max, min: %f %f\n",mmin,mmax);
      fwrite(&(mmin),sizeof(float),1,fp);
      fwrite(&(mmax),sizeof(float),1,fp);

      // Allocate memory to hold converted image 
      byteimg = (ushort *) malloc(npixs * sizeof(ushort));
      // Convert image from float to ushort
      fact = ((float)USHRT_MAX) / (mmax-mmin);
      ptr = img[i].data;
      ptrout = byteimg;
      for (k = 0 ; k < npixs ; k++)  {
	  *ptrout++ = (ushort) ((*ptr++ - mmin) * fact);
      }	 
 
      fwrite(byteimg, sizeof(ushort),baserows*basecols, fp);
      // Free memory 
      free(byteimg);
      baserows /= subsampling; basecols /= subsampling;
  }  
*/
}

/*********************************************************************
 *
 */

KLT_Pyramid::KLT_Pyramid(FILE* fp) {
  int basecols,baserows, i;//, k, tmp, npixs;
  //ushort *byteimg, *ptr;
  //float fact, *ptrout, mmin, mmax;
  
  fread(&basecols,sizeof(int),1,fp);
  fread(&baserows,sizeof(int),1,fp);
  fread(&subsampling,sizeof(int),1,fp);
  fread(&nLevels,sizeof(int),1,fp);
  
  setupStructure(basecols, baserows);
  
  for (i = 0 ; i < nLevels; i++) 
    fread(img[i].data, sizeof(float), ncols[i]*nrows[i], fp);

  /*
  for (i = 0 ; i < nLevels; i++) {
    fread(&(mmin),sizeof(float),1,fp);
    fread(&(mmax),sizeof(float),1,fp);
    npixs = basecols * baserows;
    // Allocate memory to hold inputted image 
    byteimg = (ushort *) malloc(npixs * sizeof(ushort));
    tmp = fread(byteimg,sizeof(ushort),npixs,fp);    
    if (tmp != baserows*basecols) 
      KLTError("Pyramid file read error\n");

    fact = ((float)USHRT_MAX) / ((float)(mmax-mmin));
    ptr = byteimg;
    ptrout = img[i].data;
    for (k = 0 ; k < npixs ; k++)  {
      *ptrout++ = ((float)(*ptr++) + mmin*fact) / fact;
    }
    
    // Free memory 
    free(byteimg);
    baserows /= subsampling; basecols /= subsampling;
  }
  */
}


KLT_ColorPyramid::KLT_ColorPyramid(int basecols, int baserows, int subsampling, int nlevel) {
  _r = new KLT_Pyramid(basecols, baserows, subsampling, nlevel);
  _g = new KLT_Pyramid(basecols, baserows, subsampling, nlevel);
  _b = new KLT_Pyramid(basecols, baserows, subsampling, nlevel);
}

KLT_ColorPyramid::~KLT_ColorPyramid() {
  delete _r; delete _g; delete _b;
}

void KLT_ColorPyramid::smoothAndComputePyramid(const QImage im, const Kernels* kern, float sigma_fact) {
  KLT_FloatImage floatImg(im.width(), im.height()), *smoothImg;
  
  floatImg.takeRed(im);
  smoothImg = floatImg.getSmoothed(kern);
  _r->computePyramid(smoothImg, sigma_fact);
  delete smoothImg;

  floatImg.takeGreen(im);
  smoothImg = floatImg.getSmoothed(kern);
  _g->computePyramid(smoothImg, sigma_fact);
  delete smoothImg;

  floatImg.takeBlue(im);
  smoothImg = floatImg.getSmoothed(kern);
  _b->computePyramid(smoothImg, sigma_fact);
  delete smoothImg;
}

/*
// SPEED: put all into one procedure, flow is same for r,g,b
int KLT_ColorPyramid::color(const float x, const float y, const int level, Vec3f& putHere) const {

  int status;
  putHere.set_r(_r->getFImage(level)->interpolate(x,y,&status));
  if (!status) return -1;
  putHere.set_g(_g->getFImage(level)->interpolate(x,y,&status));
  putHere.set_b(_b->getFImage(level)->interpolate(x,y,&status));

  return 0;
} 
*/

// SPEED: roll in gradient accesses, to
int KLT_ColorPyramid::color(const float x, const float y, const int level, Vec3f& putHere) const {
  if (x<0 || y<0) return -1;
  int xt = (int) x;  /* coordinates of top-left corner */
  int yt = (int) y;
  float ax = x - xt;
  float ay = y - yt;
  int ncols = _r->getNCols(level), nrows = _r->getNRows(level);
  int offset = (ncols*yt) + xt;
  float *ptr_r = _r->getFImage(level)->data + offset,
    *ptr_g = _g->getFImage(level)->data + offset,
    *ptr_b = _b->getFImage(level)->data + offset;
  
  
  if (xt<0 || yt<0 || xt>(ncols-2) || yt>(nrows-2)) {

    //fprintf(stderr, "(xt,yt)=(%d,%d)  imgsize=(%d,%d)\n"
    //    "(x,y)=(%f,%f)  (ax,ay)=(%f,%f)\n",
    //    xt, yt, ncols, nrows, x, y, ax, ay);
    //fflush(stderr); 
    return -1;
  }  

  Vec3f d1,d2;
  if (ay!=0) {
    d1.Set(*ptr_r + ay * (*(ptr_r+ncols) - *ptr_r),
	   *ptr_g + ay * (*(ptr_g+ncols) - *ptr_g),
	   *ptr_b + ay * (*(ptr_b+ncols) - *ptr_b));
    if (ax!=0) 
      d2.Set(*(ptr_r+1) + ay * (*(ptr_r+ncols+1) - *(ptr_r+1)),
	     *(ptr_g+1) + ay * (*(ptr_g+ncols+1) - *(ptr_g+1)),
	     *(ptr_b+1) + ay * (*(ptr_b+ncols+1) - *(ptr_b+1)));
  }
  else {
    d1.Set(*ptr_r, *ptr_g, *ptr_b);
    if (ax!=0)
      d2.Set(*(ptr_r+1),*(ptr_g+1),*(ptr_b+1));
  }
  if (ax!=0)
    Vec3f_Lerp(putHere,d1,d2,ax);
  else
    putHere = d1;

  return 0;
} 


void KLT_ColorPyramid::write(FILE* fp) const {
  assert(_r && _g && _b);
  _r->write(fp);
  _g->write(fp);
  _b->write(fp);
}

KLT_ColorPyramid::KLT_ColorPyramid(FILE* fp) {
  _r = new KLT_Pyramid(fp);
  _g = new KLT_Pyramid(fp);
  _b = new KLT_Pyramid(fp);
}

#define BRACK(a) min(max(((int)(a)),0), 255)
void KLT_ColorPyramid::writeImages(char* name) const {
  int nlevels = _r->getNLevels();
  for (int t=0; t<nlevels; t++) {
    int w = _r->getNCols(t), h = _r->getNRows(t);
    QImage im(w,h,32);
    KLT_FloatImage *rf = _r->getFImage(t), *gf = _g->getFImage(t),  *bf = _b->getFImage(t);
    for (int j=0; j<h; j++)
      for (int i=0; i<w; i++) {
	int offset = j*w+i;
	im.setPixel(i,j,qRgb(BRACK(rf->g(offset)),
			     BRACK(gf->g(offset)),
			     BRACK(bf->g(offset))));
      }
    
    char myName[250];
    sprintf(myName,"%s_p%d.png",name,t);
    im.save(myName, "PNG");
  } 
}

void KLT_ColorPyramid::writeDerivImages(char* name) const {
  int nlevels = _r->getNLevels();
  for (int t=0; t<nlevels; t++) {
    int w = _r->getNCols(t), h = _r->getNRows(t);
    QImage im(w,h,32);
    KLT_FloatImage *rf = _r->getFImage(t), *gf = _g->getFImage(t),  *bf = _b->getFImage(t);
    for (int j=0; j<h; j++)
      for (int i=0; i<w; i++) {
	int offset = j*w+i;
	im.setPixel(i,j,qRgb(BRACK(rf->g(offset)/2. + 125.),
			     BRACK(gf->g(offset)/2. + 125.),
			     BRACK(bf->g(offset)/2. + 125.)));
      }
    
    char myName[250];
    sprintf(myName,"%s_p%d.png",name,t);
    im.save(myName, "PNG");
  } 
}

//---------

void KLT_Pyramid::writeImages(char* name) const {
  int nlevels = getNLevels();
  for (int t=0; t<nlevels; t++) {
    int w = getNCols(t), h = getNRows(t);
    QImage im(w,h,32);
    const KLT_FloatImage *f = getFImage(t);
    for (int j=0; j<h; j++)
      for (int i=0; i<w; i++) {
	int offset = j*w+i;
	int val = BRACK(3. * f->v(offset));
	im.setPixel(i,j,qRgb(val, val, val));
      }
    
    char myName[250];
    sprintf(myName,"%s_edges_p%d.png",name,t);
    im.save(myName, "PNG");
  } 
}

void KLT_Pyramid::writeDerivImages(char* name) const {
  int nlevels = getNLevels();
  for (int t=0; t<nlevels; t++) {
    int w =  getNCols(t), h = getNRows(t);
    QImage im(w,h,32);
    const KLT_FloatImage *f = getFImage(t);
    for (int j=0; j<h; j++)
      for (int i=0; i<w; i++) {
	int offset = j*w+i;
	int val = BRACK(3.*f->v(offset)/2. + 125.);
	im.setPixel(i,j,qRgb(val, val, val));
      }
    
    char myName[250];
    sprintf(myName,"%s_edges_p%d.png",name,t);
    im.save(myName, "PNG");
  } 
}


int KLT_ColorPyramid::getNLevels() const {
  assert((_r->getNLevels() == _g->getNLevels())  && 
	 (_g->getNLevels()== _b->getNLevels()));
  return _r->getNLevels();
}









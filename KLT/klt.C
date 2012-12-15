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




#include <qdatetime.h>
#include <math.h>    
#include <stdlib.h>  
//#include "ccmath.h"
#include "klt.h"
#include "error.h"
#include "myMontage.h"
#include "myIMatrix.h"
//#include "myulb/setulb.h"

//#include "kltSpline.C"
#include "myassert.h"

KLT_FullPyramid::KLT_FullPyramid() {
  img = gradx = grady = NULL;
  nPyramidLevels = -1;
}

KLT_FullPyramid::KLT_FullPyramid(const QImage im, const KLT_TrackingContext* tc) {
  nPyramidLevels = -1;
  initMe(im, tc);
}

void KLT_FullPyramid::initMe(const QImage im, const KLT_TrackingContext* tc) {
  assert(!img);
  nPyramidLevels = tc->nPyramidLevels;
  img = new KLT_Pyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  Kernels kernSmooth(tc->smooth_sigma_fact);
  KLT_FloatImage floatImg(im);
  KLT_FloatImage* smoothImg = floatImg.getSmoothed(&kernSmooth);
  img->computePyramid(smoothImg, tc->pyramid_sigma_fact);
  delete smoothImg;
  gradx = new KLT_Pyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  grady = new KLT_Pyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  Kernels kern(tc->grad_sigma);
  for (int i = 0 ; i < tc->nPyramidLevels ; i++)
    img->getFImage(i)->computeGradients(&kern,gradx->getFImage(i),grady->getFImage(i));

  assert(img && gradx && grady);
}

void KLT_FullPyramid::initMeFromEdges(const QImage im, const KLT_TrackingContext* tc) {
  assert(!img);
  nPyramidLevels = tc->nPyramidLevels;
  img = new KLT_Pyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  
  KLT_FloatImage imgr(im.width(), im.height()),
    imgg(im.width(), im.height()), imgb(im.width(), im.height());
  imgr.takeRed(im); imgg.takeGreen(im); imgb.takeBlue(im);
  
  Kernels kern(tc->grad_sigma);
  KLT_FloatImage imgdrx(im.width(), im.height()),
    imgdgx(im.width(), im.height()), imgdbx(im.width(), im.height());
  KLT_FloatImage imgdry(im.width(), im.height()),
    imgdgy(im.width(), im.height()), imgdby(im.width(), im.height());
  imgr.computeGradients(&kern, &imgdrx, &imgdry);
  imgg.computeGradients(&kern, &imgdgx, &imgdgy);
  imgb.computeGradients(&kern, &imgdbx, &imgdby);
  
  KLT_FloatImage combineE(im.width(), im.height(), 
			  &imgdrx, &imgdry,
			  &imgdgx, &imgdgy,  // ONE, dby
			  &imgdbx, &imgdby); // dgx

  img->computePyramid(&combineE, tc->pyramid_sigma_fact);
  gradx = new KLT_Pyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  grady = new KLT_Pyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  for (int i = 0 ; i < tc->nPyramidLevels ; i++)
    img->getFImage(i)->computeGradients(&kern,gradx->getFImage(i),grady->getFImage(i));
  
  assert(img && gradx && grady);
}

KLT_FullPyramid::~KLT_FullPyramid() {
  if (img) {
    delete img; delete gradx; delete grady;
  }
}

//-----------------------------------------------------------------


KLT_FullCPyramid::KLT_FullCPyramid() {
  img = gradx = grady = NULL;
  nPyramidLevels = -1;
}

KLT_FullCPyramid::KLT_FullCPyramid(const QImage im, const KLT_TrackingContext* tc) {
  img = gradx = grady = NULL;
  nPyramidLevels = -1;
  initMe(im, tc);
}

void KLT_FullCPyramid::initMe(const QImage im, const KLT_TrackingContext* tc) {
  assert(!img);
  img = new KLT_ColorPyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  Kernels kernSmooth(tc->smooth_sigma_fact);
  //KLT_FloatImage floatImg(im);
  //KLT_FloatImage* smoothImg = floatImg.getSmoothed(&kernSmooth);
  img->smoothAndComputePyramid(im, &kernSmooth, tc->pyramid_sigma_fact);
  //img->computePyramid(smoothImg, tc->pyramid_sigma_fact);
  //delete smoothImg;
  gradx = new KLT_ColorPyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  grady = new KLT_ColorPyramid(im.width(), im.height(), tc->subsampling, tc->nPyramidLevels);
  Kernels kern(tc->grad_sigma);
  for (int i = 0 ; i < tc->nPyramidLevels ; i++) {
    img->r()->getFImage(i)->computeGradients(&kern,gradx->r()->getFImage(i),grady->r()->getFImage(i));
    img->g()->getFImage(i)->computeGradients(&kern,gradx->g()->getFImage(i),grady->g()->getFImage(i));
    img->b()->getFImage(i)->computeGradients(&kern,gradx->b()->getFImage(i),grady->b()->getFImage(i));
  }
  nPyramidLevels = tc->nPyramidLevels;
  assert(img && gradx && grady);
}

KLT_FullCPyramid::~KLT_FullCPyramid() {
  if (img) {
    delete img; delete gradx; delete grady;
  }
}

void KLT_FullCPyramid::write(FILE* fp) const {
  assert(img && gradx && grady && fp);
  assert(nPyramidLevels>0);
  assert(gradx->getNLevels() >0 && grady->getNLevels()>0 && img->getNLevels()>0);
  fwrite(&nPyramidLevels, sizeof(int), 1, fp);
  img->write(fp);
  gradx->write(fp);
  grady->write(fp);
}

bool KLT_FullCPyramid::load(FILE* fp, const KLT_TrackingContext* tc) {
  assert(!img && !gradx && !grady && fp);
  fread(&nPyramidLevels, sizeof(int), 1, fp);
  if (nPyramidLevels != tc->nPyramidLevels) 
    return false;
  img = new KLT_ColorPyramid(fp);
  gradx = new KLT_ColorPyramid(fp);
  grady = new KLT_ColorPyramid(fp);
  return true;
}

void KLT_FullCPyramid::writeImages(char* imgname, char* gxname, char* gyname) {
  img->writeImages(imgname);
  gradx->writeDerivImages(gxname);
  grady->writeDerivImages(gyname);
}

void KLT_FullPyramid::writeImages(char* imgname, char* gxname, char* gyname) {
  img->writeImages(imgname);
  gradx->writeDerivImages(gxname);
  grady->writeDerivImages(gyname);
}

void KLT_FullPyramid::write(FILE* fp) const {
  assert(img && gradx && grady && fp);
  fwrite(&nPyramidLevels, sizeof(int), 1, fp);
  img->write(fp);
  gradx->write(fp);
  grady->write(fp);
}

bool KLT_FullPyramid::load(FILE* fp, const KLT_TrackingContext* tc) {
  assert(!img && !gradx && !grady && fp);
  fread(&nPyramidLevels, sizeof(int), 1, fp);
  if (nPyramidLevels != tc->nPyramidLevels) 
    return false;
  img = new KLT_Pyramid(fp);
  gradx = new KLT_Pyramid(fp);
  grady = new KLT_Pyramid(fp);
  return true;
}


//-----------------------------------------------------------------


KLT_TrackingContext::KLT_TrackingContext() {
  mutualInit();
}



KLT_TrackingContext::~KLT_TrackingContext() {
  //if (pyramid_last)
  //delete pyramid_last;
}

void KLT_TrackingContext::mutualInit() {


  /* Set values to default values */
  max_iterations = 15;
  min_displacement = 0.1;
  max_residue = 6.0;
  grad_sigma = 1.0;
  smooth_sigma_fact = 1.3; // usually 1.3
  pyramid_sigma_fact = 0.9;
  //pyramid_last = NULL;
  constantWindow = true;
  pinLast = true;
  useNormalEq = true;

  nPyramidLevels = 4;
  subsampling = 2;

  smoothAlpha = .01; // usually equals 0.01 for patches, should differentiate
  smooth2Deriv =40000.; //1  //40000.;
  smooth1Deriv = 500.; ///100.; //500. for earlier term
  smooth0Deriv = 0.1; //10.;
  shape2Deriv = 5000;
  edgeWeight = 2000.;
  useImage = true;
  useEdges = true;
  useDiffScale = false;
  usePseudo = false;
  printSingulars = false;
  dumpWindows = false;
  D2mode = 0;  //  0 for original, 1 for new
  // checkWindow(); // not necessary while window is 13
  _stateOk=true;
  //_A = NULL; // DEBUG
  steinum = 0; //G!
}

void KLT_TrackingContext::copySettings(const KLT_TrackingContext* o) {


  /* Set values to default values */
  max_iterations = o->max_iterations;
  min_displacement = o->min_displacement;
  max_residue = o->max_residue;
  grad_sigma = o->grad_sigma;
  smooth_sigma_fact = o->smooth_sigma_fact; 
  pyramid_sigma_fact = o->pyramid_sigma_fact;
  //pyramid_last = NULL;
  constantWindow = o->constantWindow;
  pinLast = o->pinLast;
  useNormalEq = o->useNormalEq;

  nPyramidLevels = o->nPyramidLevels;
  subsampling = o->subsampling;

  smoothAlpha = o->smoothAlpha; 
  smooth2Deriv = o->smooth2Deriv;
  smooth1Deriv = o->smooth1Deriv; 
  smooth0Deriv = o->smooth0Deriv;
  shape2Deriv = o->shape2Deriv;
  edgeWeight = o->edgeWeight;
  useImage = o->useImage;
  useEdges = o->useEdges;
  useDiffScale = o->useDiffScale;
  usePseudo = o->usePseudo;
  printSingulars = o->printSingulars;
  dumpWindows = o->dumpWindows;
  D2mode = o->D2mode;  //  0 for original, 1 for new
  // checkWindow(); // not necessary while window is 13
  _stateOk=o->_stateOk;
  //_A = NULL; // DEBUG
}


void printDoubleArray(FILE* fp, const double* a, const int nrows, const int ncols) {
  for (int j=0; j<nrows; j++) {
    for (int i=0; i<ncols; i++) {
      fprintf(fp,"%.2f ", a[j*ncols+i]);
    }
    fprintf(fp,"\n");
  }
}

void printDoubleArray(FILE* fp, const double* a, const int nrows, const int ncols, const int offset, const int span) {
  for (int j=0; j<nrows; j++) {
    for (int i=0; i<ncols; i++) {
      fprintf(fp,"%.2f ", a[offset + j*span+i]);
    }
    fprintf(fp,"\n");
  }
}


void KLT_TrackingContext::checkWindow(int& window_width, int& window_height) {
  // Check window size (and correct if necessary) 
  if (window_width % 2 != 1) {
    window_width = window_width+1;
    //KLTWarning("(KLTChangeTCPyramid) Window width must be odd.  "
    //         "Changing to %d.\n", window_width);
  }
  if (window_height % 2 != 1) {
    window_height = window_height+1;
    //KLTWarning("(KLTChangeTCPyramid) Window height must be odd.  "
    //"Changing to %d.\n", window_height);
  }
  if (window_width < 9) {
    window_width = 9;
    //KLTWarning("(KLTChangeTCPyramid) Window width must be at least three.  \n"
    //            "Changing to %d.\n", window_width);
  }
  if (window_height < 9) {
    window_height = 9;
    //KLTWarning("(KLTChangeTCPyramid) Window height must be at least three.  \n"
    //         "Changing to %d.\n", window_height);
  }

}


/* 
float KLT_TrackingContext::pyramidSigma() const {
  return (pyramid_sigma_fact * subsampling);
}
*/




/*
void KLTCreateFile(KLT_PixelType *img, int ncols, int nrows, KLT_TrackingContext tc,
		   FILE* fp1, FILE* fp2, FILE* fp3) {
      _KLT_FloatImage tmpimg, floatimg;
      _KLT_Pyramid pyramid, pyramid_gradx, pyramid_grady;
      int i;
      
      floatimg = _KLTCreateFloatImage(ncols, nrows);
      tmpimg = _KLTCreateFloatImage(ncols, nrows);
      _KLTToFloatImage(img, ncols, nrows, tmpimg);
      _KLTComputeSmoothedImage(tmpimg, _KLTComputeSmoothSigma(tc), floatimg);
      _KLTWriteFloatImageToPGM(floatimg,"test.pgm"); 
      pyramid = _KLTCreatePyramid(ncols, nrows, tc->subsampling, tc->nPyramidLevels);
      _KLTComputePyramid(floatimg, pyramid, tc->pyramid_sigma_fact);
      pyramid_gradx = _KLTCreatePyramid(ncols, nrows, tc->subsampling, tc->nPyramidLevels);
      pyramid_grady = _KLTCreatePyramid(ncols, nrows, tc->subsampling, tc->nPyramidLevels);
      for (i = 0 ; i < tc->nPyramidLevels ; i++)
	  _KLTComputeGradients(pyramid->img[i], tc->grad_sigma, 
			       pyramid_gradx->img[i],
			       pyramid_grady->img[i]);

      _KLTWritePyramid(fp1,pyramid);
      _KLTWritePyramid(fp2,pyramid_gradx);
      _KLTWritePyramid(fp3,pyramid_grady);
}


void KLT_StorePyramids(KLT_TrackingContext tc, FILE* fp1, FILE* fp2, FILE* fp3) {
    tc->pyramid_last = _KLTReadPyramid(fp1);
    tc->pyramid_last_gradx = _KLTReadPyramid(fp2);
    tc->pyramid_last_grady = _KLTReadPyramid(fp3);
}
*/


//-----------------------------------------------------------------------------------------------



// make sure P is appropriately time shifted
int KLT_TrackingContext::calculateMtNorm(const int n, const int i, const int j, const int l, const Vec2f* P,
		 double *G, double *KS, double *MS) {
  if (n==0) {
    calculateKStripNorm0(KS,i,j,P);
    rmmult(MS, G, KS, 3, 2, 6);
    return 0;
  }
  else if (n==l-1) {
    calculateKStripNorml(KS,n,i,j,P);
    rmmult(MS, G, KS, 3, 2, 6);
    return (2*(n-2));
  }
  else {
    calculateKStripNorm(KS,n,i,j,P);
    rmmult(MS, G, KS, 3, 2, 6);
    return (2*(n-1));
  }
}


void KLT_TrackingContext::calculateKStripNorm(double* K, const int n, const double i, const double j, const Vec2f* P) const {
  memset(K,0, 12*sizeof(double));
  double fx,fy,g;

  // calculate fx,fy
  double dx,dy,root;
  dx = P[n+1].x() - P[n-1].x(); dy = P[n+1].y() - P[n-1].y();
  root = sqrt(dx*dx+dy*dy);
  if (root < .001) { 
    printf("fuck\n");
  return; }
  root = 1./root;
  fx = fy = root;
  root = root*root*root;
  fx -= dx*dx*root;
  fy -= dy*dy*root;
  
  g = -dx*dy*root; 
  
  K[2] = K[9] = 1.; 
  K[0] = -( K[4]  = i*fx - j*g );
  K[1] = -( K[5]  = i*g - j*fy );
  K[6] = -( K[10] = i*g + j*fx );
  K[7] = -( K[11] = i*fy + j*g );
}
void KLT_TrackingContext::calculateKStripNorm0(double* K, const double i, const double j, const Vec2f* P) const {
  memset(K,0, 12*sizeof(double));
  double fx,fy,g;
  
  // calculate fx,fy
  double dx,dy,root;
  dx = P[1].x() - P[0].x(); dy = P[1].y() - P[0].y();
  root = sqrt(dx*dx+dy*dy);
  if (root < .001)  { 
    printf("fuck\n");
  return; }
  root = 1./root;
  fx = fy = root;
  root = root*root*root;
  fx -= dx*dx*root;
  fy -= dy*dy*root;
  
  g = -dx*dy*root; 

  K[0] = 1. - ( K[2] = i*fx - j*g );
  K[1] =    - ( K[3] = i*g - j*fy );
  K[6] =    - ( K[8] = i*g + j*fx );
  K[7] = 1. - ( K[9] = i*fy + j*g );
}
void KLT_TrackingContext::calculateKStripNorml(double* K, const int n, const double i, const double j, const Vec2f* P) const {
  memset(K,0, 12*sizeof(double));
  double fx,fy,g;
 
  // calculate fx,fy
  double dx,dy,root;
  dx = P[n].x() - P[n-1].x(); dy = P[n].y() - P[n-1].y();
  root = sqrt(dx*dx+dy*dy);
  if (root < .001)  { 
    printf("fuck\n");
  return; }
  root = 1./root;
  fx = fy = root;
  root = root*root*root;
  fx -= dx*dx*root;
  fy -= dy*dy*root;
  
  g = -dx*dy*root;  
  
  K[4]  = 1. - ( K[2] = -(i*fx - j*g) );
  K[5]  =    - ( K[3] = -(i*g - j*fy) );
  K[10] =    - ( K[8] = -(i*g + j*fx) );
  K[11] = 1. - ( K[9] = -(i*fy + j*g) );
}

/*
void KLT_TrackingContext::calculateKStrip(double* K, const double i, const double j, const double invc) const {
  double ic = i*invc*.5, 
    jc = j*invc*.5;
   memset(K,0, 12*sizeof(double));
   K[0]=-ic; K[1]=jc; K[2]=1; K[4]=ic; K[5]=-jc;
   K[6]=-jc; K[7]=-ic; K[9]=1; K[10]=jc; K[11]=ic;
}
void KLT_TrackingContext::calculateKStrip0(double* K, const double i, const double j, const double invc) const {
  double ic = i*invc, // SPEED!  pre-divide by 2, multiply by 2 at 0 and l
    jc = j*invc;
   memset(K,0, 12*sizeof(double));
   K[0]=1.-ic; K[1]=jc; K[2]=ic; K[3]=-jc;
   K[6]=-jc; K[7]=1.-ic; K[8]=jc; K[9]=ic;
}
void KLT_TrackingContext::calculateKStripl(double* K, const double i, const double j, const double invc) const {
  double ic = i*invc, // SPEED!  Pass it in
    jc = j*invc;
   memset(K,0, 12*sizeof(double));
   K[2]=-ic; K[3]=jc; K[4]=1+ic; K[5]=-jc;
   K[8]=-jc; K[9]=-ic; K[10]=jc; K[11]=1+ic;
}
*/


void KLT_TrackingContext::calculateKStrip(double* K, const double i, const double j, const double invsubs) const {
  double ic = i*.25, jc = j*.25;
   memset(K,0, 12*sizeof(double));
   K[0]=-ic; K[1]=jc; K[2]=1; K[4]=ic; K[5]=-jc;
   K[6]=-jc; K[7]=-ic; K[9]=1; K[10]=jc; K[11]=ic;
   for (int i=0; i<12; i++)
     K[i] *= invsubs;
}
void KLT_TrackingContext::calculateKStrip0(double* K, const double i, const double j, const double invsubs) const {
  double ic = i*.5,  jc = j*.5;
   memset(K,0, 12*sizeof(double));
   K[0]=1.-ic; K[1]=jc; K[2]=ic; K[3]=-jc;
   K[6]=-jc; K[7]=1.-ic; K[8]=jc; K[9]=ic;
   for (int i=0; i<10; i++)
     K[i] *= invsubs;
}
void KLT_TrackingContext::calculateKStripl(double* K, const double i, const double j, const double invsubs) const {
  double ic = i*.5, jc = j*.5;
   memset(K,0, 12*sizeof(double));
   K[2]=-ic; K[3]=jc; K[4]=1+ic; K[5]=-jc;
   K[8]=-jc; K[9]=-ic; K[10]=jc; K[11]=1+ic;
   for (int i=0; i<12; i++)
     K[i] *= invsubs;
}

#define _K(n,a,b) K[(b)*(2*l) + (n)*2 + (a)] // n, (0,1) for (x,y) of n, (0,1) for (x,y) row
void KLT_TrackingContext::calculateK(double* K, const int n, const double i, 
				     const double j, const double invc, const int l) const {
  double ic, jc;
  memset(K,0, 4*l*sizeof(double));
  if (n==0) {
    ic = i*invc;
    jc = j*invc;
    
    _K(0,0,0) = _K(0,1,1)  = 1.-ic;
    _K(0,1,0) = jc;
    _K(0,0,1) = -jc;

    _K(1,0,0) = _K(1,1,1)  = ic;
    _K(1,1,0) = -jc;
    _K(1,0,1) = jc;
  }
  else if (n==l-1) {
    ic = i*invc;
    jc = j*invc;

    _K(l-2,0,0) = _K(l-2,1,1) = -ic;
    _K(l-2,1,0) = jc;
    _K(l-2,0,1) = -jc;

    _K(l-1,0,0) = _K(l-1,1,1) = 1+ic;
    _K(l-1,1,0) = -jc;
    _K(l-1,0,1) = jc;

  }
  else {
    ic = i*invc*.5;
    jc = j*invc*.5;

    _K(n,0,0) = _K(n,1,1) = 1;
    
    _K(n-1,0,0) = _K(n-1,1,1) = -ic;
    _K(n-1,1,0) = jc;
    _K(n-1,0,1) = -jc;
    
    _K(n+1,0,0) = _K(n+1,1,1) = ic;
    _K(n+1,0,1) = jc;
    _K(n+1,1,0) = -jc;
  }
}


/*
void KLT_TrackingContext::h(const int n, const float i, const float j, const int l, const float invc,
       const Vec2f* P, Vec2f& result) const { 
  Vec2f Tn, Nn;
  if (n==l-1) {
    Vec2f_Sub(Tn, P[l-1],P[l-2]);
    Tn *= invc;
  }
  else if (n==0) {
    Vec2f_Sub(Tn, P[1],P[0]);
    Tn *= invc;
  }
  else {
    Vec2f_Sub(Tn, P[n+1],P[n-1]);
    Tn *= invc;
    Tn *= 0.5f;
  }
  Nn.Set(-Tn.y(), Tn.x());

  result = P[n];
  result.Inc(i*Tn.x(), i*Tn.y());
  result.Inc(j*Nn.x(), j*Nn.y());
}
*/

void KLT_TrackingContext::h(const int n, const float i, const float j, const int l, const float invsubs,
       const Vec2f* P, Vec2f& result) const { 
  Vec2f Tn, Nn;
  if (n==l-1) {
    Vec2f_Sub(Tn, P[l-1],P[l-2]);
    Tn *= .5f;   // .5 for spacing, 
  }
  else if (n==0) {
    Vec2f_Sub(Tn, P[1],P[0]);
    Tn *= .5f;
  }
  else {
    Vec2f_Sub(Tn, P[n+1],P[n-1]);
    Tn *= .25f;
  }
  Nn.Set(-Tn.y(), Tn.x());

  result = P[n];
  result.Inc(i*Tn.x(), i*Tn.y());
  result.Inc(j*Nn.x(), j*Nn.y());
  result *= invsubs; // for subsampling
}

// consider only normalizing Nn
void KLT_TrackingContext::hnorm(const int n, const float i, const float j, const int l,
       const Vec2f* P, Vec2f& result) const { 
  Vec2f Tn, Nn;
  if (n==l-1) {
    Vec2f_Sub(Tn, P[l-1],P[l-2]);
    Tn.Normalize();
  }
  else if (n==0) {
    Vec2f_Sub(Tn, P[1],P[0]);
    Tn.Normalize();
  }
  else {
    Vec2f_Sub(Tn, P[n+1],P[n-1]);
    Tn.Normalize();
  }
  Nn.Set(-Tn.y(), Tn.x());

  result = P[n];
  result.Inc(i*Tn.x(), i*Tn.y());
  result.Inc(j*Nn.x(), j*Nn.y());
}

/*void KLT_TrackingContext::rmmultplus(double *rm,double *a,double *b,int n,int m,int l) const { 
  double z,*q0,*p,*q; int i,j,k;
  q0=(double *)calloc(m,sizeof(double));
  for(i=0; i<l ;++i,++rm){
    for(k=0,p=b+i; k<m ;p+=l) q0[k++]= *p;
    for(j=0,p=a,q=rm; j<n ;++j,q+=l){
      for(k=0,z=0.; k<m ;) z+= *p++ * q0[k++];
      *q+=z;
    }
  }
  free(q0);
  }*/

void KLT_TrackingContext::visMatrix(const double *A, const int rows, const int cols) const {
  int F = 1, val;
  double e = .00001;
  QImage im(F*cols, F*rows, QImage::Format_Indexed8);
  im.setNumColors(2);
  im.setColor(0,qRgb(255,255,255));
  im.setColor(1,qRgb(0,0,0));
  for (int r=0; r<rows; r++)
    for (int c=0; c<cols; c++) {
      if (fabs(A[r*cols + c]) > e) 
	val = 1;
      else
	val = 0;
      for (int i=0; i<F; i++)
	for (int j=0; j<F; j++)
	  im.setPixel(c*F + i, r*F + j, val);
    }
  im.save("vis.png","PNG");
}

void KLT_TrackingContext::assertSym(const double *A, const int cols) const {
  double e = .00001;
  for (int r=0; r<cols; r++)
    for (int c=0; c<cols; c++) 
      assert(fabs(A[r*cols + c] - A[c*cols+r]) < e);
}

//-----------------------------------------------------------------------------------------------


void KLT_TrackingContext::calculateNewK(double* K, const double i, const double j, const double invsubs) const {
  double ic = i*.5, jc = j*.5;
  // this is for non-expanding normal window
  /*
   K[0]=1.-ic; K[1]=jc; K[2]=ic; K[3]=-jc;
   K[4]=-jc; K[5]=1.-ic; K[6]=jc; K[7]=ic;
   for (int i=0; i<8; i++)   // SPEED: unroll this, two times as many multiplies as necessary
     K[i] *= invsubs;
  */

  //K[0]=invsubs*(1.-ic); K[1]=jc;              K[2]=invsubs*ic; K[3]=-jc;
  //K[4]=-jc;             K[5]=invsubs*(1.-ic); K[6]=jc;         K[7]=invsubs*ic;

  K[0]=invsubs*(1.-ic); K[1]=jc;   K[2]=invsubs*ic; K[3]=-jc;
  K[4]=K[3];            K[5]=K[0]; K[6]=jc;         K[7]=K[2];
}


void KLT_TrackingContext::hnew(const int n, const float i, const float j, const float invsubs, 
       const Vec2f* P, Vec2f& result) const { 
  Vec2f Tn(P[n+1],P[n]);
  Tn *= .5f; //invc;
  result = P[n];
  result.Inc(i * Tn.x(),  i * Tn.y());
  result *= invsubs;
  result.Inc(j * -Tn.y(), j * Tn.x());
  //result *= invsubs;  this is for non-expanding normal window, comment previous instance
}

void KLT_TrackingContext::hnew2(const int cp, const int c, const int t, const float i, const float j, const float invsubs, 
       const Vec2f* Z, Vec2f& result) const { 

  int index = _mt->packed_index(c, cp, t);
  Vec2f Tn;

  if (t==0 || t==_mt->_crs[c]._nPoints-2) {     // don't know that two points adjacent in memory
    int index2 = _mt->packed_index(c, cp+1, t);
    Vec2f_Sub(Tn,Z[index2], Z[index]);
    result = Z[index];
  }
  else {                                       // we know they're adjacent
    index = _mt->packed_index(c, cp, t);
    Vec2f_Sub(Tn,Z[index+1], Z[index]);
    result = Z[index];
  }
  

  Tn *= .5f; //invc;
  if (i!=0)
    result.Inc(i * Tn.x(),  i * Tn.y());
  result *= invsubs;
  if (j!=0)
    result.Inc(j * -Tn.y(), j * Tn.x());
  //result *= invsubs;  this is for non-expanding normal window, comment previous instance
}


double KLT_TrackingContext::createTestSolution(const Vec2f* Z, const double* sol, const int numVec, const int l, Vec2f* Z2) {
  memcpy(Z2, Z, numVec*l*sizeof(Vec2f));
  int offset=0, varI=0, t, n;
  double maxStep = 0, tmp;
  for (t=0; t<numVec; t++)
    for (n=0; n<l; n++, offset++, varI+=2) {
      Z2[offset].Inc(sol[varI], sol[varI+1]);
      if ((tmp=fabs(sol[varI])) > maxStep)
	maxStep = tmp;
      if ((tmp=fabs(sol[varI+1])) > maxStep)
	maxStep = tmp;
    }
  return maxStep;
}




// d is search direction, p is current location
void KLT_TrackingContext::projectToTR(double *x, const double* d, const double* p, double trustRadius, int n) {

  double a,b,c, det, r1, tau;
  a = vecSqrLen(n,d);
  b = 2. * vecDot(n,d,p);
  c = vecSqrLen(n,p) - trustRadius*trustRadius;

  det = b*b-4.*a*c;
  assert(det >=0);
  det = sqrt(det);

  r1 = (-b + det) / (2.*a);
  if (r1>0) 
    tau = r1;
  else {
    tau = (-b - det) / (2.*a);
    assert(tau > 0);
  }

  assert(tau >0 && tau <1);
  vecAssign(n,x,d);
  vecTimesScalar(n,x,tau);
  vecAddEqual(n,x,p);
}

#define ABSMIN(a,b) ((fabs(a))<(fabs(b)) ? (a) : (b))
void KLT_TrackingContext::projectToTR2(double *x, const double* d, const double* p, double trustRadius, int n) {
  int i;
  double min = DBL_MAX, k1, k2;
  
  for (i=0; i<n; i++) {

    k1 = (trustRadius-p[i]) / d[i];
    k2 = (-trustRadius-p[i]) / d[i];
    if (k1>0 && k1<min)
	min = k1;
    if (k2>0 && k2<min)
      	min = k2;
    min = MIN(trustRadius, min);
    /*
    k1 = (trustRadius-p[i]) / d[i];
    k2 = (-trustRadius-p[i]) / d[i];
    k1 = ABSMIN(k1,k2);
    k2 = fabs(k1);
    if (k2 < absmin) {
      min_i = i;
      absmin = k2; min = k1;
    }
    */
  }

  
  vecAssign(n,x,d);
  vecTimesScalar(n,x,min);
  vecAddEqual(n,x,p);
}


// PASSED in error is pointer.  NULL for nonlinear.  Otherwise, the actual error
#define DELETE_SS delete[] p; delete[] r; delete[] d; delete[] Bd; delete[] pn; 
int KLT_TrackingContext::steihaugSolver(GenKeeper* B, double* x, const double trustRadius, bool* boundaryHit, double* error2ptr) {

  QTime time;
  time.start();
  
  int maxIter = 5000;
  int n = B->numVar(), j=0;
  double *p = new double[n], *r = new double[n], *d = new double[n], *Bd = new double[n], *pn = new double[n];
  const double* g = B->g();
  double dBd, rr, alpha, beta, rrn, tmp;
  //double error2; 
  double gradnorm2 = vecSqrLen(n,g);
  assert(!isnan(gradnorm2));
  double error2;
  if (error2ptr == NULL) {
    error2 = MIN(.1, pow(gradnorm2,.25)); 
    error2 *= error2*sqrt(gradnorm2);
    //error2 = MIN(1., error2);  
    //error2 = MAX(error2,.000001);
    error2ptr = &error2;  // to NOP later assignment of this
  }
  else
    error2 = *error2ptr;  //ACCURACY just set error very small
  //error2 = .00000001;  //G!
  printf("%d variables, error threshold = %.7f, grad norm : %.5f\n",n,error2, sqrt(gradnorm2));

  //char name[200]; 
  //sprintf(name,"case%.3d.txt",steinum++);
  //FILE* fpe = fopen(name, "w");

  *boundaryHit = false;
  memset(p,0,n*sizeof(double));
  vecAssign(n, r, g);
  vecTimesScalar(n,r,-1.); // CAREFUL: for createMatrices[2,3] ONLY since J'r NOT negated there, it is in createMatrices
  B->handleConstraints(r);
  vecAssign(n, d, r);
  if (B->precond())
    B->precondVec(d);

  //vecAssign(n, r, g);
  //B->matVecMult(x, Bd);
  //vecDiffEqual(n, r, Bd);
  //vecAssign(n, d, r);

  tmp = vecSqrLen(n, r);
  if (tmp <= error2) {
    vecAssign(n,x,p);
    printf("error2 low %f, inner converge\n",tmp);
    DELETE_SS;
    *error2ptr = tmp;
    return 0;
  }
  //printf("gradient length: %f\n",tmp);

  rr = vecDot(n,r,d);  
  while (j<maxIter && _stateOk) { // max iterations of CG 
    j++;
    B->matVecMult(d,Bd);
    dBd = vecDot(n,d,Bd);
    //printf("rr %f dBd %f\n",rr,dBd);
    if (dBd < 0) {
      printf("negative curvature, projecting, inner converge\n");
      assert(0);
      *boundaryHit = true;
      projectToTR2(x,d,p,trustRadius, n);
      DELETE_SS;
      printf("solved in %d  iterations \n",j);
      *error2ptr = tmp;
      return j;
    }

    alpha = rr / dBd;
    vecAssign(n,pn,d);
    vecTimesScalar(n,pn,alpha);
    vecAddEqual(n,pn,p);
    
    //if (vecSqrLen(n,pn) >= trustRadius*trustRadius) {
    
    if (fabs(vecAbsMax(n,pn)) >= trustRadius) {
      printf("Hit trust radius boundary %f, inner converge\n", trustRadius);
      *boundaryHit = true;
      projectToTR2(x,d,p,trustRadius, n);
      //sanityCheck(n,x); 
      DELETE_SS;
      printf("solved in %d  iterations \n",j);
      *error2ptr = tmp;
      return j;
    }
    
    
    vecTimesScalar(n, Bd, alpha); // Bd now alpha*Bd
    vecDiffEqual(n,r,Bd);
    B->handleConstraints(r);
    //printf("iteration %d, error %f\n",j,vecSqrLen(n,r)); 

    //fprintf(fpe,"%.10e %.10e\n",vecSqrLen(n,r), vecSqrLen(n,pn)); 
    tmp=vecSqrLen(n,r);
    if (tmp < error2) { 
      printf("error  near nil (%f), inner converge\n", sqrt(tmp));
      vecAssign(n,x,pn);
      printf("final step, error2 %.4e, 2length %f\n",tmp, vecSqrLen(n,pn));
      //vecPrint(n,x);
      printf("solved in %d  iterations \n",j);
      printf("%d time taken to solve\n",time.elapsed());
      DELETE_SS; 
      *error2ptr = tmp;
      return j;
    }

    
    if (!B->precond()) {
      rrn = vecDot(n,r,r);
      beta = rrn/rr;
      rr = rrn;
      vecTimesScalar(n,d,beta);
      vecAddEqual(n, d, r);
    }
    else {
      vecAssign(n,Bd,r); // we use Bd to mean s in Shewchuk
      B->precondVec(Bd);
      rrn = vecDot(n,r,Bd);
      beta = rrn/rr;
      rr = rrn;
      vecTimesScalar(n,d,beta);
      vecAddEqual(n, d, Bd);      
    }
    
    vecAssign(n,p,pn);
  }
  

  vecAssign(n,x,pn);

  //B->outputMat("B.dat");   
  //B->outputg("g.txt");  
  //std::exit(0);

  printf("%d iterations hit, leaving\n", maxIter);
  printf("error2 %f\n",tmp); 
  printf("%d time taken to solve\n",time.elapsed());
  //std::exit(0);
  *error2ptr = vecSqrLen(n, x);
  return maxIter;
    //printf("solved in %d  iterations \n",j);
    //if (j > 100) exit(0);
  
  
}



int KLT_TrackingContext::longCurveTrack(const KLT_FullCPyramid** pyrms, const int numFrames, const int numP,
					const int wwidth, const int wlh, const int whh, Vec2f* Plocs, Vec2f* Z) {
  //Z is numFrames*numP(from 1 to numFrames) long, Plocs is numP long


  if ((pinLast && numFrames <2) || (!pinLast && numFrames <1)) {
    printf("internal too short span %d %d\n", pinLast, numFrames);
    return KLT_TRACKED;
  }
  //int n,t,offset;
  int r,val=KLT_TRACKED;
  //float fsubs = float(subsampling);
  int curveSampling = 5; 

  for (r = nPyramidLevels - 1 ; r >= 0 ; r--)  {
    curveSampling *= 2;
  }

  for (r = nPyramidLevels - 1 ; r >= 0 ; r--)  {
    curveSampling /= 2;

    printf("level %d\n",r);
    double invsubs = 1. / pow((double)subsampling,r);

    val = internalLongCurveTrack(pyrms,numFrames,numP,wwidth,wlh, whh, 
    			 Plocs,Z, r, invsubs, curveSampling); // assuming original spacing is 2

   

    if (val==KLT_NOT_FOUND)
      return val;
  }
  
  return val;
}


void KLT_TrackingContext::setupLongCurveTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, MultiTrackData* mt, bool redo) {
  _ttask = LCT2;
  __pyrms = pyrms; __pyrmsE = pyrmsE;
  _mt = mt; __redo = redo;
}

void KLT_TrackingContext::run() {
  if (_ttask==LCT2) {
    longCurveTrack2(__pyrms,__pyrmsE,__redo);
  }
  else if (_ttask==SPLINE) {
    splineTrack(__pyrms, __pyrmsE, __redo);
  }
}

void KLT_TrackingContext::runNoThread() {
  if (_ttask==LCT2) {
    longCurveTrack2(__pyrms,__pyrmsE,__redo);
  }
  else if (_ttask==SPLINE) {
    splineTrack(__pyrms, __pyrmsE, __redo);
  }
}


void KLT_TrackingContext::safeTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE,
		const int level, const double invsubs, const int curveSampling, const int hold, 
		const int maxIterations) {
  bool res=true;
  while (1) {
    res = internalLongCurveTrack2(pyrms,pyrmsE,level,invsubs,curveSampling,hold,maxIterations);

    if (res)
      break;

    else {  // state not ok, probably nudge

      while (1) {
	_mt->_z_mutex->lock();
	assert(_stateOk==false);  // only can make it true here
	_stateOk = true;
	_mt->_z_mutex->unlock();
    _mt->_z_wait->wait(_mt->_z_mutex);
	//_mt->_z_wait->wait();  // make sure dragging is done
	printf("state invalidated, starting smooth run\n");
	bool temp1 = useImage, temp2 = useEdges;
	useImage = false;  
	useEdges = false;
	res = internalLongCurveTrack2(pyrms,pyrmsE,level,invsubs,1,1,3);  // smooth!
	useImage = temp1; useEdges = temp2;
	_mt->_z_mutex->lock();
	if (_stateOk) {
	  _mt->_z_mutex->unlock();	
	  break;
	}
	_mt->_z_mutex->unlock();
      } 

      // if this is a smoothing run anyways, don't need to do again
      _mt->_z_mutex->lock();
      if (curveSampling==1 && _stateOk) break;
      _mt->_z_mutex->unlock();	

    }
  }
}

#define DELETE_LCT2 delete[] pyrms; delete pyrmsE;
void KLT_TrackingContext::longCurveTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, 
					 bool redo) {

  if ((pinLast && _mt->_numFrames <2) || (!pinLast && _mt->_numFrames <1)) {
    printf("internal too short span %d %d\n", pinLast, _mt->_numFrames);
    DELETE_LCT2;
    return;
  }
  QTime time;
  time.start();

  int mylevels = nPyramidLevels;

  //int n,t,offset;
  int r;
  //float fsubs = float(subsampling);

  if (redo) {
    bool temp1 = useImage, temp2 = useEdges;
    useImage = false; 
    useEdges = false;
    safeTrack2(pyrms, pyrmsE, 0, 1., 10,5,10);  // smooth special
    useImage = temp1; useEdges = temp2;
  }
  
  else {  
    bool temp1 = useImage, temp2 = useEdges;
    useImage = false;  
    useEdges = false;
    if (useEdges && useImage) initEdgeMins(pyrmsE,0,1); 
    safeTrack2(pyrms, pyrmsE, 0, 1., 1, 1,3);  // smooth 1,1,3
    useImage = temp1; useEdges = temp2;
    mylevels = 2;
  }

  int curveSampling = 4;  
  curveSampling <<= mylevels;
  for (r = mylevels - 1 ; r >= 0 ; r--)  {
    curveSampling >>= 1;
    
    printf("level %d\n",r);
    double invsubs = 1. / pow((double)subsampling,r);

    if (useEdges && useImage) initEdgeMins(pyrmsE,r,invsubs); 
    safeTrack2(pyrms, pyrmsE, r, invsubs, curveSampling, 1,1000); // track   //  was recently 5, made slower
  }
  
  
  //useImage = false;
  //val = internalLongCurveTrack2(pyrms, mt, 0, 1., 3); 
  printf("%f time taken to track\n",float(time.elapsed())/1000.f);
  DELETE_LCT2;

}

void KLT_TrackingContext::initEdgeMins(const KLT_FullPyramid** pyrmsE, 
				       const int level, const double invsubs) {

  int c, l, hw, n, ut = _mt->_numFrames-1, i, validPoint, di, ict;
  double eval, eval2;
  Vec2f loc;
  float finvsubs = (float) invsubs;

  for (c=0; c < _mt->_nCurves; c++) { // iterate over curves
    CurveRecord* cr = _mt->_crs + c;
    if (!cr->_useEdges) continue;
    l = cr->_nPoints;
    
    if (cr->_edgeMin) 
      delete[] cr->_edgeMin;
    cr->_edgeMin = new double[2*l-1];
    ict = 0;
    for (n=0; n<l-1; n++) { // iterate over curve points

      if (n==l-2)
	hw=2;
      else
	hw=1;

      for (i = 0 ; i <= hw ; ++i, ++ict)  {  

	validPoint = 1;
	assert(cr->_edgeMin);

	hnew2(n,c,0,i,0,finvsubs, _mt->_Z, loc);
	eval = pyrmsE[0]->img->getFImage(level)->interpolate(loc.x(), loc.y(), &di);
	validPoint = MIN( di, validPoint);
	hnew2(n,c,ut,i,0,finvsubs,_mt->_Z, loc);
	eval2 = pyrmsE[ut]->img->getFImage(level)->interpolate(loc.x(), loc.y(), &di);
	validPoint = MIN( di, validPoint);
	if (validPoint==1)
	  cr->_edgeMin[ict] = MAX(eval, eval2);
	else
	  cr->_edgeMin[ict] = -1;
	//assert(cr->_edgeMin[ict] ==-1);
	
      } // iterate over pixels
    } // end iterate over curve points
  } // end iterate over curves
  
}

double KLT_TrackingContext::robust(double r) const {
  return (r / (1. + r));
}


#define DELETE_CM3 delete[] theta; delete[] thetas; delete[] thetas1;  delete[] thetas0; delete[] thetaE;
double KLT_TrackingContext::createMatrices3(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE,
					    Vec2f* Z, const int level, const double invsubs, MultiKeeper* keep,
					    MultiKeeper *imgKeep, MultiKeeper *edgeKeep) {

  int numFrames = _mt->_numFrames;
  double *theta = new double[numFrames], *thetaE = new double[numFrames], *thetas = new double[numFrames], 
    *thetas1 = new double[numFrames], *thetas0 = new double[numFrames]; 
  int  validPoint, di;
  int n, t, i, j, c, l, ict;
  int ul = numFrames-1;
  int hw;
  double G_t[6],  G_t1[6]; // 3x2 Jacobian matrices
  Vec3f col1, col, dummy; 
  Vec2f loc;
  double K[8], grad[12], grad1[12], eval;
  const double D0grad[8] = { 1.,0, -1.,0,
			     0,1., 0,-1. };
  const double D1grad[16] = { 1.,0,-1.,0,  -1.,0,1.,0,  
			      0,1.,0,-1.,  0,-1.,0,1 };
  const double D2grad[24] = { 1.,0,-2.,0,1.,0,   -1.,0,2.,0,-1.,0, 
			      0,1.,0,-2.,0,1.,   0,-1.,0,2.,0,-1. };
  double sjac1[12], sjac2[12];

  double invcw=0, invew=0, invnw2=0, invnw1=0, invnw0=0;
  for (c=0; c < _mt->_nCurves; c++) {
    l = _mt->_crs[c]._nPoints;
    invcw += double((l*2-1)*( _mt->_trackWidths[c].y() - _mt->_trackWidths[c].x() + 1));
    invnw2 += double(l-2); invnw1 += double(l-1); invnw0 += double(l);
    invew += double(l*2-1);
  }
  //printf("original # samples %f\n",invcw);
  invcw  = 1./invcw; invew = 1./invew;  // color, edges
  invnw2 = 1./invnw2; invnw1 = 1./invnw1; invnw0 = 1./invnw0; // smoothness

  double smoothCoeff2 = smooth2Deriv * invnw2;
  double smoothCoeff1 = smooth1Deriv * invnw1;
  double smoothCoeff0 = smooth0Deriv * invnw0;
  double edgeCoeff = edgeWeight; //  * invew;

  float finvsubs = (float) invsubs;
  //FILE* fp;
  bool oldD = false;


  memset(theta,0,numFrames*sizeof(double));  // image term
  memset(thetas,0,numFrames*sizeof(double)); //2nd deriv term
  memset(thetas1,0,numFrames*sizeof(double));
  memset(thetas0,0,numFrames*sizeof(double));
  memset(thetaE,0,numFrames*sizeof(double));
  //  memset(samplesLost,0,numFrames*sizeof(int));  

  /* fp = fopen("/homes/gws/aseem/A.txt", "w");
  keep->printSubs(fp);
  fclose(fp);
  exit(0); */

  int imgCompCount=0, edgeCompCount=0;
  if (imgKeep)
    imgKeep->refresh(); 
  if (edgeKeep)
    edgeKeep->refresh();
  
  for (c=0; c < _mt->_nCurves && _stateOk; c++) { // iterate over curves
    CurveRecord* cr = _mt->_crs + c;
    l = cr->_nPoints;

    for (t=0; t<_mt->_numFrames && _stateOk; ++t) {


      for (n=0, ict=0; n<l; n++) { // iterate over curve points
      
	if (n < l-1) { // image term  // second term for FD, not for CD
	  if (n==l-2)
	    hw=2;
	  else
	    hw=1;
	  for (i = 0 ; i <= hw ; ++i,++ict)  {  // 0 for FD, -hw for CD 
	    if (useImage) {
	      for (j = _mt->_trackWidths[c].x() ; j <= _mt->_trackWidths[c].y() ; ++j) { 
		
		validPoint=0;
		hnew2(n,c,t+1,i,j,finvsubs, Z, loc); 
		di = pyrms[t+1]->img->color(loc.x(),loc.y(),level, col1);
		validPoint = MIN( di, validPoint); 
		pyrms[t+1]->gradx->color(loc.x(),loc.y(), level, dummy); 
		G_t1[0]=dummy.r(); G_t1[2]=dummy.g(); G_t1[4]=dummy.b();
		pyrms[t+1]->grady->color(loc.x(),loc.y(), level, dummy);
		G_t1[1]=dummy.r(); G_t1[3]=dummy.g(); G_t1[5]=dummy.b();
		
		hnew2(n,c,t,i,j,finvsubs, Z, loc);   
		di = pyrms[t]->img->color(loc.x(),loc.y(),level, col);
		validPoint = MIN( di, validPoint); 
		pyrms[t]->gradx->color(loc.x(),loc.y(), level, dummy); 
		G_t[0]=dummy.r(); G_t[2]=dummy.g(); G_t[4]=dummy.b();
		pyrms[t]->grady->color(loc.x(),loc.y(), level, dummy);
		G_t[1]=dummy.r(); G_t[3]=dummy.g(); G_t[5]=dummy.b();
		
		if (validPoint==0) {
		  imgCompCount++;
		  col -= col1;
		  theta[t] += col.Len2();
		  
		  
		  // FD	    
		  calculateNewK(K,i,j,invsubs); 
		  if (t>0)
		    rmmult(grad,G_t,K,3,2,4);      // G_t K_t
		  if (t<ul) {
		    for (di=0; di<8; di++)
		      K[di] *=-1.;
		    rmmult(grad1, G_t1,K,3,2,4);   // - G_(t+1) K_(t+1)
		  }
		  
		  if (t==0) {
		    imgKeep->take4Key0(c, n,grad1,  col.r(),1.); // SPEED!  Make no constant version
		    imgKeep->take4Key0(c, n,grad1+4,col.g(),1.);
		    imgKeep->take4Key0(c, n,grad1+8,col.b(),1.);
		  }
		  else if (t==ul) {
		    imgKeep->take4Key1(c, n,grad,  col.r(),1.);
		    imgKeep->take4Key1(c, n,grad+4,col.g(),1.);
		    imgKeep->take4Key1(c, n,grad+8,col.b(),1.);
		  }
		  else {
		    imgKeep->take4(t-1, c, n, grad,   grad1,   col.r(),1.);
		    imgKeep->take4(t-1, c, n, grad+4, grad1+4, col.g(),1.);
		    imgKeep->take4(t-1, c, n, grad+8, grad1+8, col.b(),1.);
		  }
		  

		}
		
	      } // j
	    } // end if useImage
	    
	    // EDGE TERM
	    CurveRecord* cr = _mt->_crs + c;
	    //if (t>0 && useEdges)
	    //printf("shit\n");
	    if (useEdges && useImage && t>0 && cr->_useEdges &&  cr->_edgeMin[ict] > 0) {
	      assert(cr->_edgeMin);

	      // SPEED: could incorporate into image term to avoid recalculating loc
	      hnew2(n,c,t,i,0,finvsubs,Z,loc);
	      eval = pyrmsE[t]->img->getFImage(level)->interpolate(loc.x(), loc.y(), &validPoint); 
	      G_t[0] = pyrmsE[t]->gradx->getFImage(level)->interpolate(loc.x(), loc.y(), &di);
	      G_t[1] = pyrmsE[t]->grady->getFImage(level)->interpolate(loc.x(), loc.y(), &di);

	      if (validPoint==1) {
		edgeCompCount++;
		double tmp = eval / cr->_edgeMin[ict];
		thetaE[t-1] += tmp*tmp;
		//assert(tmp*tmp < 10000); 

		calculateNewK(K,i,0,invsubs);
		rmmult(grad,G_t, K, 1,2,4);
		edgeKeep->take4Special(t-1,c,n,grad,eval, 
				   edgeCoeff / (cr->_edgeMin[ict]*cr->_edgeMin[ict])  );
	      }
	      
	      /*
	      hnew2(n,c,t,i,2,finvsubs,_mt, Z, loc);  // +2j
	      di = pyrms[t]->img->color(loc.x(),loc.y(),level, col1);
	      validPoint = MIN( di, validPoint); 
	      pyrms[t]->gradx->color(loc.x(),loc.y(), level, dummy); 
	      G_t[0]=dummy.r(); G_t[2]=dummy.g(); G_t[4]=dummy.b();
	      pyrms[t]->grady->color(loc.x(),loc.y(), level, dummy);
	      G_t[1]=dummy.r(); G_t[3]=dummy.g(); G_t[5]=dummy.b();
	      calculateNewK(K,i,2,invsubs); 
	      rmmult(grad1,G_t,K,3,2,4);      // G_1 K_1

	      hnew2(n,c,t,i,-2,finvsubs,_mt, Z, loc);  // -2j
	      di = pyrms[t]->img->color(loc.x(),loc.y(),level, col);
	      validPoint = MIN( di, validPoint); 
	      pyrms[t]->gradx->color(loc.x(),loc.y(), level, dummy); 
	      G_t[0]=dummy.r(); G_t[2]=dummy.g(); G_t[4]=dummy.b();
	      pyrms[t]->grady->color(loc.x(),loc.y(), level, dummy);
	      G_t[1]=dummy.r(); G_t[3]=dummy.g(); G_t[5]=dummy.b();
	      calculateNewK(K,i,-2,invsubs); 
	      rmmult(grad,G_t,K,3,2,4);      // - G_2 K_2

	      if (validPoint==0) {
		col1 -= col;
		thetaE[t-1] -= robust(col1.r()*col1.r() / cr->_edgeMin[ict] + 
				      col1.g()*col1.g() / cr->_edgeMin[ict+1] +
				      col1.b()*col1.b() / cr->_edgeMin[ict+2]);
		
		for (di=0; di<12; ++di)
		  grad1[di] -=grad[di];
		
		keep->take4Special(t-1, c, n, cr->_edgeMin[ict],   grad, col1.r(), edgeCoeff); 
		keep->take4Special(t-1, c, n, cr->_edgeMin[ict+1], grad+4, col1.g(), edgeCoeff);
		keep->take4Special(t-1, c, n, cr->_edgeMin[ict+2], grad+8, col1.b(), edgeCoeff);
	      }
	      */

	    } // end edge term
	    
	  } // i

	} // if imageTerm

	const Vec2f* tptr[3]; 
	const Vec2f* t1ptr[3];
	int index, lnvar;
	if (n==0) {
	  tptr[0] = t1ptr[0] = NULL;
	  index = _mt->packed_index_n(c, n, lnvar);
	  tptr[1] = Z + index + t*lnvar;  t1ptr[1] = tptr[1] + lnvar;
	  index = _mt->packed_index(c, n+1, t);
	  tptr[2] = Z + index;  t1ptr[2] = tptr[2] + cr->_nVars;
	  
	}
	else if (n==1) {
	  index = _mt->packed_index_n(c, n-1, lnvar);
	  tptr[0] = Z + index + t*lnvar;  t1ptr[0] = tptr[0] + lnvar;
	  index = _mt->packed_index(c, n, t);
	  tptr[1] = Z + index;  t1ptr[1] = tptr[1] + cr->_nVars;
	  tptr[2] = tptr[1]+1;  t1ptr[2] = tptr[2] + cr->_nVars;
	}
	else if (n==l-2) {
	  index = _mt->packed_index(c, n-1, t);
	  tptr[0] = Z + index;  t1ptr[0] = tptr[0] + cr->_nVars;
	  tptr[1] = tptr[0]+1;  t1ptr[1] = tptr[1] + cr->_nVars;
	  index = _mt->packed_index_n(c, n+1, lnvar);
	  tptr[2] = Z + index + t*lnvar;  t1ptr[2] = tptr[2] + lnvar;
	}
	else if (n==l-1) {
	  tptr[0] = t1ptr[0] = NULL; // this is never used, though we could calculate it
	  index = _mt->packed_index_n(c, n, lnvar);
	  tptr[1] = Z + index + t*lnvar;  t1ptr[1] = tptr[1] + lnvar;
	  tptr[2] = t1ptr[2] = NULL;
	}
	else {
	  index = _mt->packed_index(c, n-1, t);
	  tptr[0] = Z + index;  t1ptr[0] = tptr[0] + cr->_nVars;
	  tptr[1] = tptr[0]+1;  t1ptr[1] = tptr[1] + cr->_nVars;
	  tptr[2] = tptr[1]+1;  t1ptr[2] = tptr[2] + cr->_nVars;
	}

	
	Vec2f res;

	if (n < l-1) { // D1 smoothness term 
	
	  if (oldD) { // old D1
	    res = *(tptr[1]); res -= *(tptr[2]);
	    res -= *(t1ptr[1]); res += *(t1ptr[2]);
	    thetas1[t] += res.Len2();
	
	
	    if (t==0) {
	      keep->take4Key0(c, n, D1grad+4, res.x(), smoothCoeff1);
	      keep->take4Key0(c, n, D1grad+12, res.y(), smoothCoeff1);
	    }
	    else if (t==ul) {
	      keep->take4Key1(c, n, D1grad, res.x(), smoothCoeff1);
	      keep->take4Key1(c, n, D1grad+8, res.y(), smoothCoeff1);
	    }
	    else {
	      keep->take4(t-1, c, n, D1grad, D1grad+4, res.x(), smoothCoeff1);
	      keep->take4(t-1, c, n, D1grad+8, D1grad+12, res.y(), smoothCoeff1);
	    }	  	  
	  }
	
	  else  { // new D1 
	    Vec2f delta1,delta2;
	    Vec2f_Sub(delta1,*(tptr[1]), *(tptr[2]));
	    Vec2f_Sub(delta2, *(t1ptr[1]), *(t1ptr[2]));
	    double resid = delta1.Len2() - delta2.Len2();
	    thetas1[t] += resid*resid;

	    if (t>0) {
	      sjac1[0] = 2. * delta1.x(); sjac1[1] = 2. * delta1.y();
	      sjac1[2] = -sjac1[0]; sjac1[3] = -sjac1[1]; 
	    }
	    if (t<numFrames) {
	      sjac2[0] = -2. * delta2.x(); sjac2[1] = -2. * delta2.y();
	      sjac2[2] = -sjac2[0]; sjac2[3] = -sjac2[1]; 
	    }

	    if (t==0)
	      keep->take4Key0(c, n, sjac2, resid, smoothCoeff1);
	    else if (t==ul) 
	      keep->take4Key1(c, n, sjac1, resid, smoothCoeff1);
	    else
	      keep->take4(t-1, c, n, sjac1, sjac2, resid, smoothCoeff1);
	  } // end length-preserving term
	
	} // end D1 smoothness term

	if (n>0 && n<l-1) { // D2 smoothness term   

	  if (D2mode==0) {  // original D2 term
	    Vec2f res1, res2;
	    res1 = *(tptr[1]); res1 *=-2.f;
	    res1 += *(tptr[0]); res1 += *(tptr[2]);
	    res2 = *(t1ptr[1]); res2 *=-2.f;
	    res2 += *(t1ptr[0]); res2 += *(t1ptr[2]);
	  
	    res1 -= res2;
	    thetas[t] += res1.Len2();
	  
	    if (t==0) {
	      keep->take6Key0(c, (n-1), D2grad+6,  res1.x(), smoothCoeff2);
	      keep->take6Key0(c, (n-1), D2grad+18, res1.y(), smoothCoeff2);
	    }
	    else if (t==ul) { 
	      keep->take6Key1(c, (n-1), D2grad,    res1.x(), smoothCoeff2);
	      keep->take6Key1(c, (n-1), D2grad+12, res1.y(), smoothCoeff2);
	    
	    }
	    else {
	      keep->take6(t-1, c, (n-1), D2grad, D2grad+6,     res1.x(), smoothCoeff2);
	      keep->take6(t-1, c, (n-1), D2grad+12, D2grad+18, res1.y(), smoothCoeff2);
	    }
	  }
	
	  else if (D2mode==1) {  // new 4th order one
	    Vec2f D1 = *(tptr[1]), D2 = *(t1ptr[1]);
	    D1 *= -2.f; D2 *= -2.f;
	    D1 += *(tptr[2]); D1 += *(tptr[0]);
	    D2 += *(t1ptr[2]); D2 += *(t1ptr[0]);
	  
	    double resid = D1.Len2() - D2.Len2();
	    thetas[t] += resid*resid;

	    if (t>0) {
	      sjac1[0] =  2. * D1.x();  sjac1[1] =  2. * D1.y();
	      sjac1[2] = -4. * D1.x();  sjac1[3] = -4. * D1.y();
	      sjac1[4] =  2. * D1.x();  sjac1[5] =  2. * D1.y();
	    }
	    if (t<numFrames) {
	      sjac2[0] = -2. * D2.x();  sjac2[1] = -2. * D2.y();
	      sjac2[2] =  4. * D2.x();  sjac2[3] =  4. * D2.y();
	      sjac2[4] = -2. * D2.x();  sjac2[5] = -2. * D2.y();
	    }

	    if (t==0)
	      keep->take6Key0(c, (n-1), sjac2, resid, smoothCoeff2);
	    else if (t==ul) 	
      keep->take6Key1(c, (n-1), sjac1, resid, smoothCoeff2);
	    else
	      keep->take6(t-1, c, (n-1), sjac1, sjac2, resid, smoothCoeff2);
	  }
	
	} // end D2 smoothness term
	
	// D0 smoothness term
	res = *(tptr[1]); res -= *(t1ptr[1]);
	thetas0[t] += res.Len2();


	if (t==0) {
	  keep->take2Key0(c, n, D0grad+2, res.x(), smoothCoeff0);
	  keep->take2Key0(c, n, D0grad+6, res.y(), smoothCoeff0);
	}
	else if (t==ul) {
	  keep->take2Key1(c, n, D0grad, res.x(), smoothCoeff0);
	  keep->take2Key1(c, n, D0grad+4, res.y(), smoothCoeff0);
	}
	else {
	  keep->take2(t-1, c, n, D0grad, D0grad+2, res.x(), smoothCoeff0);
	  keep->take2(t-1, c, n, D0grad+4, D0grad+6, res.y(), smoothCoeff0);
	}	 
	
      } // n
    } // t


  } // c

  // here is where you would divide imgKeep, edgeKeep appropriately
  if (imgCompCount>0) {
    invcw = double(_mt->_numFrames) / double(imgCompCount);
    keep->add(imgKeep, invcw);

  }
  if (edgeCompCount>0) {
    invew = double(_mt->_numFrames) / double(edgeCompCount);
    keep->add(edgeKeep, invew);
    edgeCoeff *= invew; // for later theta computation
  }

  
  float thetaSum = 0;
  for (i=0; i<numFrames; ++i) 
    thetaSum += theta[i]*invcw + thetas1[i]*smoothCoeff1 + 
      thetas[i]*smoothCoeff2 + thetas0[i]*smoothCoeff0 + thetaE[i]*edgeCoeff;

  if (_stateOk) {
    fprintf(stdout,"%9f = ",thetaSum);
    
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",theta[i]*invcw);
    fprintf(stdout,"\n            ");
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",thetas[i]*smoothCoeff2);
    fprintf(stdout,"\n            ");
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",thetas1[i]*smoothCoeff1);
    fprintf(stdout,"\n            ");
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",thetas0[i]*smoothCoeff0);
    fprintf(stdout,"\n            ");
    for (i=0; i<numFrames-1; ++i) 
      fprintf(stdout,"%.3f + ",thetaE[i]*edgeCoeff);
    printf("\n");
  }
  


  //keep->fill_e(e);
  //keep->print();
  
  /*
  fp = fopen("/homes/gws/aseem/Znew.txt","w"); 
  keep->printFull(fp);
  fclose(fp);
  fp = fopen("/homes/gws/aseem/enew.txt","w");
  keep->printE(fp);
  fclose(fp);
  //exit(0);
  */

  DELETE_CM3;
  return thetaSum; 


}


#define DELETE_CM2 delete[] theta; delete[] thetas; delete[] thetas1;  delete[] thetas0; 
double KLT_TrackingContext::createMatrices2(const KLT_FullCPyramid** pyrms, const int numFrames, const int l,
				       const int, const int wlh, const int whh, const Vec2f* Plocs, Vec2f* Z,const int level, 
				       const double invsubs, Keeper* keep) {
  assert(0);
  double *theta = new double[numFrames], *thetas = new double[numFrames], *thetas1 = new double[numFrames], *thetas0 = new double[numFrames]; 
  int  validPoint, di;
  int n, t, i, j;
  int ul = numFrames-1;
  //double* ZArray = new double[numFrames*2*l];
  int hw;// = wwidth/2;
  double   invcw = 1./double((l*2-1)*(whh-wlh+1));  // (l-1)*(hw+1) for FD, l*wwidth for CD
  double invnw2 = 1./double(l-2), invnw1 = 1./double(l-1), invnw0=1./double(l);
  double G_t[6],  G_t1[6]; // 3x2 Jacobian matrices
  Vec3f col1, col, dummy; 
  Vec2f loc;
  double K[12], grad[18], grad1[18];  // 8,12,12 for FD
  const double D0grad[8] = { 1.,0, -1.,0,
			     0,1., 0,-1. };
  const double D1grad[16] = { 1.,0,-1.,0,  -1.,0,1.,0,  
			      0,1.,0,-1.,  0,-1.,0,1 };
  const double D2grad[24] = { 1.,0,-2.,0,1.,0,   -1.,0,2.,0,-1.,0, 
			      0,1.,0,-2.,0,1.,   0,-1.,0,2.,0,-1. };
  double sjac1[12], sjac2[12];
  //double D[36], KS[12], M[18], Mt[18], Mm1[18], S[20], vp[2], C[12], U[4] = {1.,0.,0.,1.};
  double smoothCoeff2 = smooth2Deriv * invnw2;
  double smoothCoeff1 = smooth1Deriv * invnw1;
  double smoothCoeff0 = smooth0Deriv * invnw0;
  //float finvc = 1./spacing;  //pow(2,-l+2);
  //double dinvc = 1./spacing;
  float finvsubs = (float) invsubs;
  //FILE* fp;
  bool oldD = false;
  //int D2mode = 0;

  memset(theta,0,numFrames*sizeof(double));  // image term
  memset(thetas,0,numFrames*sizeof(double)); //2nd deriv term
  memset(thetas1,0,numFrames*sizeof(double));
  memset(thetas0,0,numFrames*sizeof(double));
  //  memset(samplesLost,0,numFrames*sizeof(int));  

  /* fp = fopen("/homes/gws/aseem/A.txt", "w");
  keep->printSubs(fp);
  fclose(fp);
  exit(0); */
  

  // CORRECTNESS: visualize this stuff, are double-counting some samples, not others, re-consider wwidth/hw
  // if we later want to divide by number of valid points, need to down separate nt loop for image and smoothness
  for (t=0; t<numFrames; ++t) {

    for (n=0; n<l; n++) {
      
      if (useImage && n < l-1) { // image term  // second term for FD, not for CD
	if (n==l-2)
	  hw=2;
	else
	  hw=1;
	for (j = wlh ; j <= whh ; ++j) { 
	  for (i = 0 ; i <= hw ; ++i)  {  // 0 for FD, -hw for CD 
	    
	    validPoint=0;
	    hnew(n,i,j,finvsubs,Z+t*l, loc); // hnew and no l for FD, h for CD 
	    di = pyrms[t+1]->img->color(loc.x(),loc.y(),level, col1);
	    validPoint = MIN( di, validPoint); 
	    pyrms[t+1]->gradx->color(loc.x(),loc.y(), level, dummy); 
	    G_t1[0]=dummy.r(); G_t1[2]=dummy.g(); G_t1[4]=dummy.b();
	    pyrms[t+1]->grady->color(loc.x(),loc.y(), level, dummy);
	    G_t1[1]=dummy.r(); G_t1[3]=dummy.g(); G_t1[5]=dummy.b();

	    if (t>0)
	      hnew(n,i,j,finvsubs,Z+(t-1)*l, loc);   // hnew and no l for FD, h for CD 
	    else
	      hnew(n,i,j,finvsubs,Plocs, loc);       // hnew and no l for FD, h for CD 
	    di = pyrms[t]->img->color(loc.x(),loc.y(),level, col);
	    validPoint = MIN( di, validPoint); 
	    pyrms[t]->gradx->color(loc.x(),loc.y(), level, dummy); 
	    G_t[0]=dummy.r(); G_t[2]=dummy.g(); G_t[4]=dummy.b();
	    pyrms[t]->grady->color(loc.x(),loc.y(), level, dummy);
	    G_t[1]=dummy.r(); G_t[3]=dummy.g(); G_t[5]=dummy.b();

	    if (validPoint==0) {
	      col -= col1;
	      theta[t] += col.Len2();
	      
	      
	      // FD	    
	      calculateNewK(K,i,j,invsubs); 
	      //for (di=0; di<8; di++)         // scale jacobian
	      //K[di] *=invcw;
	      if (t>0)
		rmmult(grad,G_t,K,3,2,4);      // G_t K_t
	      if (t<ul) {
		for (di=0; di<8; di++)
		  K[di] *=-1.;
		rmmult(grad1, G_t1,K,3,2,4);   // - G_(t+1) K_(t+1)
	      }
	      //col *= float(invcw);           // scale residual
	      if (t==0) {
		keep->take4Key0(2*n,grad1,  col.r(),invcw);
		keep->take4Key0(2*n,grad1+4,col.g(),invcw);
		keep->take4Key0(2*n,grad1+8,col.b(),invcw);
	      }
	      else if (t==ul) {
		keep->take4Key1(2*n,grad,  col.r(),invcw);
		keep->take4Key1(2*n,grad+4,col.g(),invcw);
		keep->take4Key1(2*n,grad+8,col.b(),invcw);
	      }
	      else {
		keep->take4(t-1, 2*n, grad,   grad1,   col.r(),invcw);
		keep->take4(t-1, 2*n, grad+4, grad1+4, col.g(),invcw);
		keep->take4(t-1, 2*n, grad+8, grad1+8, col.b(),invcw);
	      }
	      

	      // CD
	      /*
	      if (n==0) {
		calculateKStrip0(K,i,j,invsubs);
		offset=0;

	      }
	      else if (n==l-1) {
		 calculateKStripl(K,i,j,invsubs);
		 offset = 2*(n-2);
	      }
	      else {
		calculateKStrip(K,i,j,invsubs);
		offset = 2*(n-1);
	      }

	      //for (di=0; di<12; di++)         // scale jacobian
	      //K[di] *=invcw;
	      if (t>0)
		rmmult(grad,G_t,K,3,2,6);      // G_t K_t
 
	      if (t<ul) {
		for (di=0; di<12; di++)
		  K[di] *=-1.;
		rmmult(grad1, G_t1,K,3,2,6);   // - G_(t+1) K_(t+1)
	      }
	
	      if (t==0) {
		keep->take6Key0(offset,grad1,  col.r(),invcw);
		keep->take6Key0(offset,grad1+6,col.g(),invcw);
		keep->take6Key0(offset,grad1+12,col.b(),invcw);
	      }
	      else if (t==ul) {
		keep->take6Key1(offset,grad,  col.r(),invcw);
		keep->take6Key1(offset,grad+6,col.g(),invcw);
		keep->take6Key1(offset,grad+12,col.b(),invcw);
	      }
	      else {
		keep->take6(t-1, offset, grad,   grad1,   col.r(),invcw);
		keep->take6(t-1, offset, grad+6, grad1+6, col.g(),invcw);
		keep->take6(t-1, offset, grad+12, grad1+12, col.b(),invcw);
	      }
	      */

	    }
	    
	  } // i
	} // j
      } // if imageTerm

      const Vec2f *tptr = ((t==0) ? Plocs+n : Z + (t-1)*l + n), *t1ptr = Z + t*l + n;
      Vec2f res;

      if (n < l-1) { // D1 smoothness term
	
	if (oldD) { // old D1
	  res = *tptr; res -= *(tptr+1);
	  res -= *t1ptr; res += *(t1ptr+1);
	  thetas1[t] += res.Len2();
	
	
	  if (t==0) {
	    keep->take4Key0(2*n, D1grad+4, res.x(), smoothCoeff1);
	    keep->take4Key0(2*n, D1grad+12, res.y(), smoothCoeff1);
	  }
	  else if (t==ul) {
	    keep->take4Key1(2*n, D1grad, res.x(), smoothCoeff1);
	    keep->take4Key1(2*n, D1grad+8, res.y(), smoothCoeff1);
	  }
	  else {
	    keep->take4(t-1, 2*n, D1grad, D1grad+4, res.x(), smoothCoeff1);
	    keep->take4(t-1, 2*n, D1grad+8, D1grad+12, res.y(), smoothCoeff1);
	  }	  	  
	}
	
	else { // new D1
	  Vec2f delta1,delta2;
	  Vec2f_Sub(delta1,*tptr, *(tptr+1));
	  Vec2f_Sub(delta2, *t1ptr, *(t1ptr+1));
	  double resid = delta1.Len2() - delta2.Len2();
	  thetas1[t] += resid*resid;

	  if (t>0) {
	    sjac1[0] = 2. * delta1.x(); sjac1[1] = 2. * delta1.y();
	    sjac1[2] = -sjac1[0]; sjac1[3] = -sjac1[1]; 
	  }
	  if (t<numFrames) {
	    sjac2[0] = -2. * delta2.x(); sjac2[1] = -2. * delta2.y();
	    sjac2[2] = -sjac2[0]; sjac2[3] = -sjac2[1]; 
	  }

	  if (t==0)
	    keep->take4Key0(2*n, sjac2, resid, smoothCoeff1);
	  else if (t==ul) 
	    keep->take4Key1(2*n, sjac1, resid, smoothCoeff1);
	  else
	    keep->take4(t-1, 2*n, sjac1, sjac2, resid, smoothCoeff1);
	} // end length-preserving term
	
      } // end D1 smoothness term

      if (n>0 && n<l-1) { // D2 smoothness term 

	if (D2mode==0) {  // original D2 term
	  Vec2f res1, res2;
	  res1 = *tptr; res1 *=-2.f;
	  res1 += *(tptr-1); res1 += *(tptr+1);
	  res2 = *t1ptr; res2 *=-2.f;
	  res2 += *(t1ptr-1); res2 += *(t1ptr+1);
	  
	  res1 -= res2;
	  thetas[t] += res1.Len2();
	  
	  if (t==0) {
	    keep->take6Key0(2*(n-1), D2grad+6,  res1.x(), smoothCoeff2);
	    keep->take6Key0(2*(n-1), D2grad+18, res1.y(), smoothCoeff2);
	  }
	  else if (t==ul) { 
	    keep->take6Key1(2*(n-1), D2grad,    res1.x(), smoothCoeff2);
	    keep->take6Key1(2*(n-1), D2grad+12, res1.y(), smoothCoeff2);
	    
	  }
	  else {
	    keep->take6(t-1, 2*(n-1), D2grad, D2grad+6,     res1.x(), smoothCoeff2);
	    keep->take6(t-1, 2*(n-1), D2grad+12, D2grad+18, res1.y(), smoothCoeff2);
	  }
	}
	
	else if (D2mode==1) {  // new 4th order one
	  Vec2f D1 = *tptr,D2 = *t1ptr;
	  D1 *= -2.f; D2 *= -2.f;
	  D1 += *(tptr+1); D1 += *(tptr-1);
	  D2 += *(t1ptr+1); D2 += *(t1ptr-1);
	  
	  double resid = D1.Len2() - D2.Len2();
	  thetas[t] += resid*resid;

	  if (t>0) {
	    sjac1[0] =  2. * D1.x();  sjac1[1] =  2. * D1.y();
	    sjac1[2] = -4. * D1.x();  sjac1[3] = -4. * D1.y();
	    sjac1[4] =  2. * D1.x();  sjac1[5] =  2. * D1.y();
	  }
	  if (t<numFrames) {
	    sjac2[0] = -2. * D2.x();  sjac2[1] = -2. * D2.y();
	    sjac2[2] =  4. * D2.x();  sjac2[3] =  4. * D2.y();
	    sjac2[4] = -2. * D2.x();  sjac2[5] = -2. * D2.y();
	  }

	  if (t==0)
	    keep->take6Key0(2*(n-1), sjac2, resid, smoothCoeff2);
	  else if (t==ul) 
	    keep->take6Key1(2*(n-1), sjac1, resid, smoothCoeff2);
	  else
	    keep->take6(t-1, 2*(n-1), sjac1, sjac2, resid, smoothCoeff2);
	}
	
      } // end D2 smoothness term
      
      // D0 smoothness term
      res = *tptr; res -= *t1ptr;
      thetas0[t] += res.Len2();
      
      if (t==0) {
	keep->take2Key0(2*n, D0grad+2, res.x(), smoothCoeff0);
	keep->take2Key0(2*n, D0grad+6, res.y(), smoothCoeff0);
      }
      else if (t==ul) {
	keep->take2Key1(2*n, D0grad, res.x(), smoothCoeff0);
	keep->take2Key1(2*n, D0grad+4, res.y(), smoothCoeff0);
      }
      else {
	keep->take2(t-1, 2*n, D0grad, D0grad+2, res.x(), smoothCoeff0);
	keep->take2(t-1, 2*n, D0grad+4, D0grad+6, res.y(), smoothCoeff0);
      }	 
      
    } // n
  } // t


  
  float thetaSum = 0;
  for (i=0; i<numFrames; ++i) 
    thetaSum += theta[i]*invcw + thetas1[i]*smoothCoeff1 + 
      thetas[i]*smoothCoeff2 + thetas0[i]*smoothCoeff0;
  fprintf(stdout,"%9f = ",thetaSum);
  
  for (i=0; i<numFrames; ++i) 
    fprintf(stdout,"%.3f + ",theta[i]*invcw);
  fprintf(stdout,"\n            ");
  for (i=0; i<numFrames; ++i) 
    fprintf(stdout,"%.3f + ",thetas[i]*smoothCoeff2);
  fprintf(stdout,"\n            ");
  for (i=0; i<numFrames; ++i) 
    fprintf(stdout,"%.3f + ",thetas1[i]*smoothCoeff1);
  fprintf(stdout,"\n            ");
  for (i=0; i<numFrames; ++i) 
    fprintf(stdout,"%.3f + ",thetas0[i]*smoothCoeff0);
  
  printf("\n");


  //keep->fill_e(e);
  //keep->print();
  
  /*
  fp = fopen("/homes/gws/aseem/Znew.txt","w"); 
  keep->printFull(fp);
  fclose(fp);
  fp = fopen("/homes/gws/aseem/enew.txt","w");
  keep->printE(fp);
  fclose(fp);
  //exit(0);
  */

  DELETE_CM2;
  return thetaSum;
}


#define ILCT2_DELETE   delete[] x; delete[] Z2; delete keep1; delete keep2; if (imgKeep) delete imgKeep; if (edgeKeep) delete edgeKeep;
bool KLT_TrackingContext::internalLongCurveTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, 
						 const int level, const double invsubs, const int curveSampling, const int hold, 
						 const int maxIterations) {
  
  int iteration=0;
  
  bool toContinue = true;
  //int numVec = _mt->_numFrames-1;
  int sizeZ = _mt->size_Z();

  Vec2f* Z2 = new Vec2f[sizeZ];
  memcpy(Z2, _mt->_Z, sizeZ*sizeof(Vec2f));
  double trustRadius=10.;  // initial trust radius.    
  MultiKeeper *keep1 = new MultiKeeper(_mt, curveSampling, hold), *keep2 = new MultiKeeper(_mt, curveSampling, hold);  
  MultiKeeper *imgKeep = NULL, *edgeKeep = NULL;
  if (useImage)
    imgKeep = new MultiKeeper(keep1);
  if (useEdges)
    edgeKeep = new MultiKeeper(keep1);
  double *x = new double[keep1->numVar()];

  double prec;
  if (level==0)
    prec = .25;
  else
    prec = 1.0;
  double currTheta, newTheta;
  //keep1->printCurves(mt->_Z);
  currTheta = createMatrices3(pyrms, pyrmsE, _mt->_Z, level, invsubs, keep1, imgKeep, edgeKeep);

  
  //FILE* fp = fopen("/homes/gws/aseem/A.dat","w"); 
  //FILE* fp2 = fopen("/homes/gws/aseem/big.dat","w"); 
  //keep1->printSparse(fp);
  //keep1->fprint(fp2);
  //keep1->outputg();
  //fp = fopen("/homes/gws/aseem/A.dat","w");
  //keep1->printSubs(fp);
  //fclose(fp);
  //fclose(fp2);
  //keep1->outputg();
  //keep1->djacs.output();
  //std::exit(0);

  do {
    fprintf(stdout,"starting iteration %d\n", iteration);
    
    double ro, maxStep = 10.; // 10 is just to force into loop initially
    
    bool boundaryHit;
    while (maxStep > .1  && trustRadius >= prec && iteration<=maxIterations) {
      
      int numIter = steihaugSolver(keep1, x, trustRadius, &boundaryHit); 

      /*
      if (!res) {
	FILE* fp = fopen("/homes/gws/aseem/A.dat","w"); 
	keep1->printSparse(fp);
	keep1->outputg();
	fclose(fp);
	std::exit(0);
      }
      */

      /*fp = fopen("/homes/gws/aseem/x.txt","w");
      for (int rt=0; rt<keep1->numVar(); rt++) 
	fprintf(fp,"%.5f\n",x[rt]);
      fclose(fp);
      exit(0);*/
      //exit(0);

      // eval solution
      maxStep = keep1->createTestSol(_mt->_Z, x, Z2); // write to Z2 // change this to a function of _mt?
      printf("max step %f\n",maxStep);
      assert(maxStep < trustRadius + .00001);
      
      if (_stateOk)
	newTheta = createMatrices3(pyrms, pyrmsE, Z2, level, invsubs, keep2, imgKeep, edgeKeep);
      else
	newTheta = DBL_MAX;
      
      ro = keep1->calculateRo(x, newTheta, currTheta);
      printf("ro: %f\n",ro);
      
      printf("old %f, new %f\n", currTheta, newTheta); 
      if (ro < 0) {// .25 otherwise , or 0 
	trustRadius *= .25;
	printf("reducing trustRadius to %f\n",trustRadius);
      }
      else if (ro > .75 &&  boundaryHit) // add ro > .75 &&
	trustRadius = MIN(2.*trustRadius, 10.);

      _mt->_z_mutex->lock(); 
      if (_stateOk) {      
	if (ro <= 0) { // don't take step
	  printf("Not taking step\n");
	  keep2->refresh();
	}
	else { // step is fine
	  iteration++;
	  printf("Taking step, iteration %d\n\n", iteration);
	  MultiKeeper* kswap = keep1;
	  keep1 = keep2; 
	  keep2 = kswap;
	  keep2->refresh();
	  currTheta = newTheta;
	  	  
	  memcpy(_mt->_Z,Z2,sizeZ*sizeof(Vec2f)); // only need to copy data that can change (not doing this currently)	
	}
	_mt->_z_mutex->unlock();
      }
      else {  // don't apply update, get outta here
	_mt->_z_mutex->unlock();
	ILCT2_DELETE
	  return false;
      } // end else, stateOk block

    } // while loop
    
    if (trustRadius < prec || maxStep < .1  || iteration>maxIterations)
      toContinue=false;
    else
      toContinue=true;
    printf("\n\n");
    
    
  } while (toContinue);
  
  ILCT2_DELETE
  return true;
}



int KLT_TrackingContext::internalLongCurveTrack(const KLT_FullCPyramid** pyrms,const int numFrames, const int l,
						const int wwidth, const int wlh, const int whh, const Vec2f* Plocs, Vec2f* Z, 
						const int level, const double invsubs, const int curveSampling) {
  
  int  ok = KLT_TRACKED;
  int iteration=0, numVar;

  bool toContinue = true;
  int numVec = (pinLast) ? numFrames-1 : numFrames; 

  numVar = numVec*2*l;
  double *x = new double[numVar]; // e:1 x 2lt, 
  Vec2f *Z2 = new Vec2f[numFrames*l];
  memcpy(Z2, Z, numFrames*l*sizeof(Vec2f));
  double trustRadius=10.;  // initial trust radius.   
  assert(l>2);

  Keeper *keep1 = new Keeper(numVec, 2*l, curveSampling), *keep2 = new Keeper(numVec, 2*l, curveSampling); 
  double prec = 0.1; 
  double currTheta, newTheta;
  currTheta = createMatrices2(pyrms, numFrames, l, wwidth, wlh, whh, Plocs, Z, 
			      level, invsubs, keep1);
  do {
    fprintf(stdout,"starting iteration %d, %d active frames, %d variables\n",iteration, numVec, numVar);

    double ro, maxStep = 10.; // 10 is just to force into loop initially
    
    bool boundaryHit;
    while (maxStep >= prec && trustRadius >= prec) {

      steihaugSolver(keep1, x, trustRadius, &boundaryHit); // stei, remember to clear x
      
      // eval solution
      maxStep = keep1->createTestSol(Z,x,Z2); // write to Z2
      printf("max step %f\n",maxStep);

      newTheta = createMatrices2(pyrms, numFrames, l, wwidth, wlh, whh, Plocs, Z2, 
				 level, invsubs, keep2); 


      ro = keep1->calculateRo(x, newTheta, currTheta);
      printf("ro: %f\n",ro);

      // calculate ro
      //if (newTheta < currTheta) ro = 1.; else ro = -1.; // bad approx
      
      if (_stateOk) printf("old %f, new %f\n", currTheta, newTheta);
      if (ro < .25) {// .25 otherwise , or 0
	trustRadius *= .25;
	printf("reducing trustRadius to %f\n",trustRadius);
      }
      else if (ro > .75 &&  boundaryHit) // add ro > .75 &&
	trustRadius = MIN(2.*trustRadius, 10.);
      
      if (ro <= 0) { // don't take step
	keep2->refresh();
      }
      else { // step is fine
	iteration++;
	printf("Taking step, iteration %d\n\n", iteration);
	Keeper* kswap = keep1;
	keep1 = keep2; 
	keep2 = kswap;
	keep2->refresh();
	currTheta = newTheta;
	memcpy(Z,Z2,numVec*l*sizeof(Vec2f)); // only need to copy data that can change
	// CORRECTNESS: really unfortunate Z is float
      }
    }
    
    if (trustRadius < prec || maxStep < prec)
      toContinue=false;
    else
      toContinue=true;
    printf("\n\n");

    

    //fp = fopen("/homes/gws/aseem/Z.txt","w");
    //im.printExpanded(fp);
    //fclose(fp);
    
       
    
  } while (toContinue);
  
  
  delete[] x;
  delete keep1; delete keep2; 
  delete[] Z2;
  return ok;
}



//-----------------------------------------------------------------------------------------------

#define INTERNAL_RESHAPE_DELETEALL delete[] e; delete[] A_l; delete[] A_; delete[] x; delete[] thetas; delete [] theta;
int KLT_TrackingContext::longCurveReshape(int numFrames, int numP, Vec2f* Plocs, 
					  Vec2f* Z, Vec2f* Sb, bool pL) {
  if ((pL && numFrames <2) || (!pL && numFrames <1)) {
    printf("internal reshape too short span %d %d\n", pL, numFrames);
    return KLT_TRACKED;
  }
  int numVec = (pL) ? numFrames-1 : numFrames, l = numP;
  int numVar = numVec*2*l, l2 = 2*l;
  double *e = new double[numVar];
  bool ok = KLT_TRACKED;
  MySparseMat *A_l = new MySparseMat[numVec-1], *A_ = new MySparseMat[numVec];
  double *x = new double[numFrames*2*l];   // array version of Z
  double *theta = new double[numFrames], *thetas = new double[numFrames];
  double S[20], vp[20];
  double invn = 1./double(numP);
  int i,t,n, offset, s_offset;

  double *xplocs = new double[2*l];
  for (n=0,offset=0; n<l; n++, offset+=2) {
    xplocs[offset]=Plocs[n].x(); xplocs[offset+1]=Plocs[n].y();
  }

  double smoothCoeff2 = shape2Deriv * invn;

  memset(e,0,numVar*sizeof(double));

  printf("%d variables\n", numVar);
  for (int step=0; step < 2; step++) {         // two iterations prints final results
    memset(thetas,0,numFrames*sizeof(double));
    memset(theta,0,numFrames*sizeof(double));
    offset=0; s_offset=0;
    for (n=0; n<numFrames; ++n)
      for (t=0; t<l; ++t,++offset) {
	x[s_offset++] = Z[offset].x(); x[s_offset++] = Z[offset].y();
      }

    for (t=0; t<numVec; ++t) {  // time, starts at frame after initial frozen
      int lt2 = 2*l*t;
      A_[t].init(l2,false);
      if (t>0)
	A_l[t-1].init(l2,true);
      
      A_[t].addDiagonal(invn); // distance too; simple, huh?
      
      for (n=0; n<l; ++n) {  // space
	
	//---------------- Distance
	
	Vec2f diff(Sb[t*l + n], Z[t*l + n]); 
	//Vec2f diff(Plocs[n], Z[t*l + n]);
	theta[t] += diff.Len2();
	e[lt2 + 2*n]     += invn*diff.x();
	e[lt2 + 2*n + 1] += invn*diff.y(); 
	
	
	//---------------- Smoothness
	// THis section just computes values for thetas
	if (t==0) {
	  if (n<l-1 && n>0) {
	    Vec2f res = Plocs[n + 1];
	    Vec2f res1 = res;
	    Vec2f other = Plocs[n]; 
	    res1 -= other;
	    other *=2.f;
	    res -= other;
	    res += Plocs[n - 1];
	    res -= Z[n + 1];
	    res1 -= Z[n + 1];
	    other = Z[n]; 
	    res1 += other;
	    other *=2.f;
	    res += other;
	    res -= Z[n - 1];
	    thetas[0] += res.Len2(); 
	    //thetas1[0] += res1.Len2();
	  }
	}
	else if (t<numFrames) {
	  if (n<l-1 && n>0) {
	    Vec2f res = Z[(t-1)*l + n + 1];
	    Vec2f res1 = res;
	    Vec2f other = Z[(t-1)*l + n]; 
	    res1 -= other;
	    other *=2.f;
	    res -= other;
	    res += Z[(t-1)*l + n - 1];
	    res -= Z[t*l + n + 1];
	    res1 -= Z[t*l + n + 1];
	    other = Z[t*l+n]; 
	    res1 += other;
	    other *=2.f;
	    res += other;
	    res -= Z[t*l + n - 1];
	    thetas[t] += res.Len2(); 
	    //thetas1[t] += res1.Len2();
	  }
	}
	if (pL && t==numFrames-2) {
	  if (n<l-1 && n>0) {
	    Vec2f res = Z[t*l + n + 1];
	    Vec2f res1 = res;
	    Vec2f other = Z[t*l + n]; 
	    res1 -= other;
	    other *=2.f;
	    res -= other;
	    res += Z[t*l + n - 1];
	    res -= Z[(t+1)*l + n + 1];
	    res1 -= Z[(t+1)*l + n + 1];
	    other = Z[(t+1)*l+n]; 
	    res1 += other;
	    other *=2.f;
	    res += other;
	    res -= Z[(t+1)*l + n - 1];
	    thetas[t+1] += res.Len2(); 
	    //thetas1[t+1] += res1.Len2();
	  }
	}
	
	if (step==0) {
	  // Now we actually do real calculations
	  memset(S,0,20*sizeof(double));
	  // second derivative
	  if (n==0) {
	    S[0]=1.;  S[2]=-2.;  S[4]=1.;
	    S[11]=1.; S[13]=-2.; S[15]=1.;
	    s_offset=0;
	  }
	  else if (n==1) {
	    S[0]=-2.;  S[2]=5.;  S[4]=-4.;  S[6]=1.;
	    S[11]=-2.; S[13]=5.; S[15]=-4.; S[17]=1.;
	    s_offset=0;
	  }
	  else if (n==l-2) {
	    S[2]=1.;  S[4]=-4.;  S[6]=5.;  S[8]=-2.;
	    S[13]=1.; S[15]=-4.; S[17]=5.; S[19]=-2.;
	    s_offset=2*(l-5);
	  }
	  else if (n==l-1) {
	    S[4]=1.;  S[6]=-2.;  S[8]=1.;
	    S[15]=1.; S[17]=-2.; S[19]=1.;
	    s_offset=2*(l-5);
	  }
	  else {
	    S[0]=1.;  S[2]=-4.;  S[4]=6.;  S[6]=-4.;  S[8]=1.;
	    S[11]=1.; S[13]=-4.; S[15]=6.; S[17]=-4.; S[19]=1.;
	    s_offset=2*(n-2);
	  }
	  
	  if (!pL && t==0 && t==numFrames-1) {  // two frames, last not pinned; 
	    A_[0].addnbyk(2*n, s_offset, smoothCoeff2, S, 2, 10); 
	    rmmult(vp,S,x + s_offset, 2, 10, 1);
	    e[2*n]     += -smoothCoeff2*vp[0];
	    e[2*n + 1] += -smoothCoeff2*vp[1];
	    rmmult(vp,S,xplocs + s_offset, 2, 10, 1);
	    e[2*n]     += smoothCoeff2*vp[0];
	    e[2*n + 1] += smoothCoeff2*vp[1];
	  }
	  else if (pL && t==0 && t==numFrames-2) {
	    A_[0].addnbyk(2*n, s_offset, 2*smoothCoeff2, S, 2, 10);
	    rmmult(vp,S,x + s_offset, 2, 10, 1);
	    e[2*n]     += -2.*smoothCoeff2*vp[0];
	    e[2*n + 1] += -2.*smoothCoeff2*vp[1];
	    rmmult(vp,S,xplocs + s_offset, 2, 10, 1);
	    e[2*n]     += smoothCoeff2*vp[0];
	    e[2*n + 1] += smoothCoeff2*vp[1];
	    rmmult(vp,S,x + 2*l + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += smoothCoeff2*vp[1]; 
	  }
	  else if (t==0) {
	    A_[0].addnbyk(2*n, s_offset, 2.*smoothCoeff2, S, 2, 10);
	    rmmult(vp,S,x + s_offset, 2, 10, 1);
	    e[2*n]     += -2.*smoothCoeff2*vp[0];
	    e[2*n + 1] += -2.*smoothCoeff2*vp[1];
	    rmmult(vp,S,x + 2*l + s_offset, 2, 10, 1);
	    e[2*n]     += smoothCoeff2*vp[0];
	    e[2*n + 1] += smoothCoeff2*vp[1];
	    rmmult(vp,S,xplocs + s_offset, 2, 10, 1);
	    e[2*n]     += smoothCoeff2*vp[0];
	    e[2*n + 1] += smoothCoeff2*vp[1];
	  }
	  else if (!pL && t==numFrames-1) {
	    A_l[t-1].addnbyk(2*n, s_offset, -smoothCoeff2, S, 2, 10);
	    rmmult(vp,S,x + 2*l*(t-1) + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += smoothCoeff2*vp[1];
	    A_[t].addnbyk(2*n, s_offset, smoothCoeff2, S, 2, 10);
	    rmmult(vp,S,x + lt2 + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += -smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += -smoothCoeff2*vp[1];
	  }
	  else if (pL && t==numFrames-2) {
	    A_l[t-1].addnbyk(2*n, s_offset, -smoothCoeff2, S, 2, 10); 
	    rmmult(vp,S,x + 2*l*(t-1) + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += smoothCoeff2*vp[1];
	    A_[t].addnbyk(2*n, s_offset, 2.*smoothCoeff2, S, 2, 10);
	    rmmult(vp,S,x + lt2 + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += -2.*smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += -2.*smoothCoeff2*vp[1];
	    rmmult(vp,S,x + 2*l*(t+1) + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += smoothCoeff2*vp[1]; 
	  }
	  else {
	    A_l[t-1].addnbyk(2*n, s_offset, -smoothCoeff2, S, 2, 10);
	    rmmult(vp,S,x + 2*l*(t-1) + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += smoothCoeff2*vp[1];
	    A_[t].addnbyk(2*n, s_offset, 2.*smoothCoeff2, S, 2, 10);
	    rmmult(vp,S,x + lt2 + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += -2.*smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += -2.*smoothCoeff2*vp[1];
	    rmmult(vp,S,x + 2*l*(t+1) + s_offset, 2, 10, 1);
	    e[lt2 + 2*n]     += smoothCoeff2*vp[0];
	    e[lt2 + 2*n + 1] += smoothCoeff2*vp[1];
	  }
	}
      

      } // n
    
      if (step==0) {
	A_[t].pack(); 
	if (t>0)
	  A_l[t-1].pack();
      }
	
    } // t
  
    float thetaSum = 0;
    for (i=0; i<numFrames; ++i) 
      thetaSum += theta[i]*invn + thetas[i]*smoothCoeff2;
    fprintf(stdout,"%9f = ",thetaSum);    
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%f + ",theta[i]*invn);
    fprintf(stdout,"\n            ");
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%f + ",thetas[i]*smoothCoeff2);
    printf("\n");

    if (step==0) {
      fprintf(stdout,"solving\n");
      MyCongMat im(A_,A_l,numVec,l2);
      
      int steps = MAX(2000,2*numVar);
      memset(x,0,numFrames*l2*sizeof(double));
      double error = ConjGrad(numVar,&im,x,e,double(numVar)/4000.,&steps);
      double *sol = x;
      fprintf(stdout,"done solving in %d steps with error %f.\n",steps,error);
      
      
      
      if (error <= double(numVar)/500.) {      
	offset=0;
	int varI=0;
	for (t=0; t<numVec; t++) {
	  for (n=0; n<l; ++n, ++offset, varI+=2) 
	    Z[offset].Inc(sol[varI], sol[varI+1]);
	}
      }
      else
	printf("Problem solving, using should be\n");
    }
      
  } // end step
  
  INTERNAL_RESHAPE_DELETEALL;
  return ok;
}



//-----------------------------------------------------------------------------------------------
/*
int KLT_TrackingContext::longPatchTrack(const KLT_FullCPyramid** pyrms, int numFrames, float ax, float ay, 
					 int window_width, int window_height, double* param) {

  if ((pinLast && numFrames <2) || (!pinLast && numFrames <1)) {
    printf("internal too short span %d %d\n", pinLast, numFrames);
    return KLT_TRACKED;
  }

  checkWindow(window_width, window_height);
  double xloc, yloc;
  double fsubs = float(subsampling);
  int r, i, offset, ul = numFrames*6, ok = KLT_TRACKED;
  float fw = window_width, fh = window_height;
  int lw, lh;

  xloc = ax; yloc = ay;
  for (r = nPyramidLevels - 1 ; r >= 0 ; r--)  {
    xloc /= fsubs;  yloc /=  fsubs;
    fw /= fsubs; fh /=fsubs;
    for (i=0, offset=0; i<numFrames; i++, offset+=6) {
      param[offset+4] /=  fsubs;  param[offset+5] /=  fsubs;
    }
  }
  
  // Beginning with coarsest resolution, do ...
  for (r = nPyramidLevels - 1 ; r >= 0 ; r--)  {
    
    // Update features to current resolution
    xloc *= fsubs;  yloc *= fsubs;
    fw *= fsubs; fh *= fsubs;
    for (i=0, offset=0; i<numFrames; i++, offset+=6) {
      param[offset+4] *=  fsubs;  param[offset+5] *=  fsubs;
    }

    if (!constantWindow) {
      lw = (int) ceil(fw); lh = (int) ceil(fh);
      checkWindow(lw,lh);   
    }
    else {
      lw = window_width; lh = window_height;
    }

    // Track feature at current resolution 
    if (useNormalEq)
      ok = internalLongPatchTrack(pyrms,numFrames, xloc, yloc, lw, lh, param, r);
    else
      ok = internalLongPatchTrack2(pyrms,numFrames, xloc, yloc, lw, lh, param, r);
    //if (ok != KLT_TRACKED)
    //return ok;
    

    for (int i=0; i<ul; i++)
      assert(finite(param[i]));    
  }
  return ok;
  
}
*/
// DEBUG
//
//QRgb makeMyColor(Vec3f c) {
//  c.Clamp(0,255);
//  QColor cl(c.r(), c.g(), c.b());
//  return cl.rgb();
//}



/*
int KLT_TrackingContext::internalLongPatchTrack(const KLT_FullCPyramid** pyrms, int numFrames, float ax, float ay, 
						 int window_width, int window_height, 
						 double* param, const int level) {
  int hw = window_width/2;
  int hh = window_height/2;
  assert(pyrms[0] && pyrms[0]->img);
  int ccstatus, i, j,iteration=0, t, di;
  int numVar;
  if (!pinLast)
    numVar = numFrames*6; 
  else
    numVar = (numFrames-1)*6;
  int numVar2 = numVar*numVar;
  double* Z = new double[numVar2];
  double* e = new double[numVar];
  double* x = new double[numVar];
  double M_t[18], M_tplus1[18], M_tminus1[18], strange[36], vp[6];//, strangeLinear[16], strangeAffine[4];
  bool noMovement = false;
  float tx, ty;
  Vec3f imgSum, dummy;
  Vec3f G_t[2], G_tplus1[2], G_tminus1[2]; // [I_x I_y], thus 3x2 Jacobian matrix
  int ul = numFrames-1;
  int  ok = KLT_TRACKED;
  double cw = window_width*window_height;
  double theta[numFrames], thetas[numFrames]; // DEBUG
  Vec3f th1,th2,th3, d1,d2,d3, sub; // DEBUG
  double pLast[2];
  if (pinLast) {
    pLast[0] = param[numVar+4]; pLast[1] = param[numVar + 5];
  }
  QImage transformed[numFrames+1];
  if (dumpWindows) {
    for (t=0; t<=numFrames; t++) // DEBUG
      transformed[t].create(window_width, window_height,32);
  }

  // SPEED: could calculate this in longPatchTrack, divide squared terms / 4, linear by 2.
  if (1) {
    if (!useDiffScale)
      calculateStrange(strange,hh,hw); // if using limited patch, would need to pass this here
                  else {
      calculateStrangeLinear(strangeLinear,hh,hw);
      calculateStrangeAffine(strangeAffine,hh,hw);
                  }
  }
  //FILE* shit2 = fopen("strange","w");
  //fmatprt(shit2,strange,6,6,"%.2f ");
  //fclose(shit2);

  printf("initial param from %f %f:\n",ax,ay);
  for (i=0; i<numFrames*6; i++) {
    if (i%6==0) printf("), (");
    printf("%.2f, ",param[i]);
  }
  printf("\n");
  
  do {
    printf("\nstarting iteration %d, window %d x %d\n",iteration, window_width, window_height);
    // ZERO e,Z
    memset(Z,0,numVar2*sizeof(double));
    memset(e,0,numVar*sizeof(double));
    memset(theta,0,numFrames*sizeof(double));
    memset(thetas,0,numFrames*sizeof(double));
    
    // iterate over variables
    int numVec = (pinLast) ? numFrames-1 : numFrames; 
    for (t=0; t<numVec; t++) {
      // v_t is a vector of 6 variables, t=0 transforms a+1 image,
      // t=numFrames-1 transforms b image
      

      // iterate over window, calc 6 rows of system
      for (j = -hh ; j <= hh ; j++) { 
	for (i = -hw ; i <= hw ; i++)  { 
	  
	  // calculate image sum and gradients (equation 19)
	  if (t<ul) {
	    transform(i,j,&tx,&ty,param + 6*(t+1));
	    tx += ax; ty += ay;
	    di=pyrms[t+2]->img->color(tx,ty, level, imgSum);
	    ok = MIN( di, ok); 
	    pyrms[t+2]->gradx->color(tx,ty, level, G_tplus1[0]);	  // Calculate G_t
	    pyrms[t+2]->grady->color(tx,ty, level, G_tplus1[1]);
	  }
	  else 
	    imgSum.Set(0,0,0);
	  if (pinLast && t==numFrames-2) {// DEBUG
	    th3 = imgSum;
	    d3.Set(tx, ty,0);
	    if (dumpWindows)
	      transformed[numFrames].setPixel(i+hw, j+hh,makeMyColor(th3));
	  }

	  if (t>0)
	    transform(i,j,&tx,&ty,param + 6*(t-1));
	  else {  // frame a has identity transform
	    tx = i; ty = j;
	  }
	  tx += ax; ty += ay;
	  di=pyrms[t]->img->color(tx,ty, level, dummy);
	  ok = MIN( di, ok);
	  imgSum += dummy;
	  d1.Set(tx,ty,0);
	  th1 = dummy; // DEBUG
	  if (dumpWindows && t==0)
	    transformed[0].setPixel(i+hw, j+hh,makeMyColor(th1));
	  pyrms[t]->gradx->color(tx,ty, level, G_tminus1[0]);	  // Calculate G_t
	  pyrms[t]->grady->color(tx,ty, level, G_tminus1[1]);

	  transform(i,j,&tx,&ty,param + 6*t);
	  tx += ax; ty += ay;
	  di=pyrms[t+1]->img->color(tx,ty, level, dummy);
	  ok = MIN( di, ok);
	  if (ok!=0) {
	    delete[] Z; delete[] e; return ok;
	  }
	  imgSum -= dummy;
	  th2 = dummy; // DEBUG
	  d2.Set(tx,ty,0);
	  if (dumpWindows)
	    transformed[t+1].setPixel(i+hw, j+hh,makeMyColor(th2)); // DEBUG
	  if (t<ul)
	    imgSum -= dummy;
	  pyrms[t+1]->gradx->color(tx,ty, level, G_t[0]);	  // Calculate G_t
	  pyrms[t+1]->grady->color(tx,ty, level, G_t[1]);

	  // DEBUG
	  Vec3f_Sub(sub,th1,th2);
	  theta[t] += sub.Len2();
	  Vec3f_Sub(sub,d1,d2);
	  if (1)
	    thetas[t] += sub.Len2();
	  if (pinLast && t==numFrames-2) {
	    Vec3f_Sub(sub,th2,th3);
	    theta[t+1] += sub.Len2();
	    Vec3f_Sub(sub,d2,d3);
	    if (1)
	      thetas[t+1] += sub.Len2();
	  }
	  // SPEED: convert i,j to double or float only once
	  // Calculate M_t
	  calculateM_t(M_t, G_t, i, j);
	  if (t<ul)
	    calculateM_t(M_tplus1, G_tplus1, i, j);
	  if (t>0)
	    calculateM_t(M_tminus1, G_tminus1, i, j);

	  
	  // calculate e (right side of equation 19)
	  MTranV(e+6*t, M_t, imgSum);

	  
	  // calculate 3 terms for Z matrix
	  // SPEED: get rid of 6*t multiple multiplies
	  if (t>0) {
	    MTranM(Z + numVar*6*t + 6*(t-1), numVar, -1, M_t, M_tminus1);
	  }
	  if (t<ul) {
	    MTranM(Z + numVar*6*t + 6*t, numVar, 2, M_t, M_t); 	    //printf("mtranm %d %d\n",t,t);
	    if (t!=numFrames-2 || !pinLast) {
	      MTranM(Z + numVar*6*t + 6*(t+1), numVar, -1, M_t, M_tplus1);
	    }
	  }
	  else {
	    MTranM(Z + numVar*6*t + 6*t, numVar, 1, M_t, M_t);
	  }
	  
	  //for (di=0; di<numVar2; di++)
	  //assert(finite(Z[di]));
	} // i
      } // j


      if (1) {
	if (!useDiffScale) {
	  if (t==0) { // SPEED: t is 0, substitute it
	    add6by6(Z,  numVar*6*t + 6*(t+1), numVar, -smoothAlpha, strange);
	    vmul(vp, strange, param+6*(t+1), 6);
	    add6Vec(e+6*t, vp, smoothAlpha);
	    add6by6(Z,  numVar*6*t + 6*t, numVar, 2.*smoothAlpha, strange);
	    vmul(vp, strange, param+6*t, 6);
	    add6Vec(e+6*t, vp, -2.*smoothAlpha);
	  }
	  else if (!pinLast && t==ul) {
	    add6by6(Z,  numVar*6*t + 6*t, numVar, smoothAlpha, strange);
	    vmul(vp, strange, param+6*t, 6);
	    add6Vec(e+6*t, vp, -smoothAlpha);
	    add6by6(Z,  numVar*6*t + 6*(t-1), numVar, -smoothAlpha, strange);
	    vmul(vp, strange, param+6*(t-1), 6);
	    add6Vec(e+6*t, vp, smoothAlpha);
	  }
	  else if (pinLast && t==numFrames-2) { // also t==ul-1
	    add6by6(Z,  numVar*6*t + 6*(t-1), numVar, -smoothAlpha, strange);
	    vmul(vp, strange, param+6*(t-1), 6);
	    add6Vec(e+6*t, vp, smoothAlpha);
	    add6by6(Z,  numVar*6*t + 6*t, numVar, 2.*smoothAlpha, strange);
	    vmul(vp, strange, param+6*t, 6);
	    add6Vec(e+6*t, vp, -2.*smoothAlpha);
	    vmul(vp, strange, param+6*(t+1), 6);  // only for initials, not vars
	    add6Vec(e+6*t, vp, smoothAlpha);
	  }
	  else {
	    add6by6(Z,  numVar*6*t + 6*(t-1), numVar, -smoothAlpha, strange);
	    vmul(vp, strange, param+6*(t-1), 6);
	    add6Vec(e+6*t, vp, smoothAlpha);
	    add6by6(Z,  numVar*6*t + 6*(t+1), numVar, -smoothAlpha, strange);
	    vmul(vp, strange, param+6*(t+1), 6);
	    add6Vec(e+6*t, vp, smoothAlpha);
	    add6by6(Z,  numVar*6*t + 6*t, numVar, 2.*smoothAlpha, strange);
	    vmul(vp, strange, param+6*t, 6);
	    add6Vec(e+6*t, vp, -2.*smoothAlpha);
	  }
	}
	else {    // different scale case
	  if (t==0) { 
	    
	  }
	  else if (t==ul) {
	    
	  }
	  else {
	    
	  }
	} // end different scale case
      }
      //for (di=0; di<numVar2; di++)
      //assert(finite(Z[di]));
      
    }
    //giveArrays(Z,e,numVar); 

    // normalize system by window size // SPEED: not necessary once done
    for (i=0; i<numVar2; i++)
      Z[i] /= cw;
    for (i=0; i<numVar; i++)  
      e[i] /= cw;

    // DEBUG
    float thetaSum = 0;
    for (i=0; i<numFrames; i++) 
      thetaSum += theta[i] + thetas[i];
    printf("%9.3f = ",thetaSum/cw);
    
    for (i=0; i<numFrames; i++) 
      printf("%.3f + ",theta[i]/cw);
    printf("\n            ");
    if (1) {
      for (i=0; i<numFrames; i++) 
	printf("%.3f + ",thetas[i]/cw);
    }
    
    printf("\n");
    if (dumpWindows) {
      
      char named[200];
      for (i=0; i<numFrames; i++) {
	sprintf(named,"window%.3d.png",i);
	transformed[i].save(named,"PNG");
      }
      
      //MyMontage mont(transformed,numFrames+1,5,5);
      //sprintf(named,"windowMont.png");
      //mont.save(named);
    }
  
    
    printf("column norms: ");
    for (i=0; i<numVar; i++) {
      double thetaSum = 0;
      for (j=0; j<numVar; j++)
	thetaSum += Z[j*numVar + i] * Z[j*numVar + i];
      printf("%f, ", thetaSum / cw);
    }
    printf("\n");
    
    //printDoubleArray(stderr,Z,numVar, numVar); printf("\n");
    //printDoubleArray(stderr, e, numVar, 1);

    
    printf("solving...\n");
    FILE* shit = fopen("Z","w");
    fmatprt(shit,Z,numVar,numVar,"%5.2f ");
    fclose(shit);
    shit = fopen("e","w");
    fmatprt(shit,e,numVar,1,"%5.2f ");
    fclose(shit);
    

    // DEBUG        
    if (printSingulars) {
      double* sins = new double[numVar];
      double* copyZ = new double[numVar2];
      memcpy(copyZ,Z, numVar2*sizeof(double));
      int sinsStatus = svdval(sins,copyZ,numVar,numVar);
      assert(sinsStatus==0);
      delete[] copyZ;
      printf("singulars: ");
      for (int sini=0; sini<numVar; sini++)
	printf("%.3f ",sins[sini]);
      printf("\n");
      delete[] sins;
    }

    double *sol = NULL;
    if (!usePseudo) {
      ccstatus = solvps(Z,e,numVar); // e becomes delta v
      assert(ccstatus==0);
      sol = e;
      //printf("done solving.\n");
    }
    else {
      ccstatus = pseudoSolve(x, Z, e, numVar);
      if (ccstatus != numVar)
	printf("Deficient rank %d from %d\n",ccstatus, numVar);
      sol = x;
    }

    noMovement = true; // hypothesis to be tested
    for (i=0; i<numVar; i++) {
      param[i] += sol[i];
      if (fabs(sol[i]) > min_displacement)
	noMovement = false;
      assert(finite(param[i]));
    }
    
    for (i=0; i<numFrames*6; i++) { 
      if (i%6==0) printf("), (");
      printf("%.2f, ",param[i]);
    }
    printf("\n"); 
    if (pinLast) 
      assert(pLast[0] == param[numVar+4] && pLast[1] == param[numVar+5]);

  } while (!noMovement && iteration++ < max_iterations);

  delete[] Z; delete[] e; delete[] x;

  if (iteration >= max_iterations)
    ok = KLT_MAX_ITERATIONS;

  return ok;
}
*/

// SPEED: get rid of indexing multiplies
void KLT_TrackingContext::add6by6(double* Z, int offset, int span, double coeff, double* adder) {
  int i,j;
  for (j=0; j<6; j++) 
    for (i=0; i<6; i++)
      Z[offset + span*j + i] += coeff*adder[6*j+i];
}
/*
void KLT_TrackingContext::addnbyk(double* Z, const int span, const double coeff, 
				  const double* adder, const int n, const int k) {
  int i,j, shop, khop;
  for (j=0, shop=0, khop=0; j<n; ++j, shop+=span, khop+=k) 
    for (i=0; i<k; ++i)
      Z[shop + i] += coeff*adder[khop+i];
}

void KLT_TrackingContext::addnbyk(double* Z, const int span,  
				  const double* adder, const int n, const int k) {
  int i,j, shop, khop;
  for (j=0, shop=0, khop=0; j<n; j++, shop+=span, khop+=k) 
    for (i=0; i<k; ++i)
      Z[shop + i] += adder[khop+i];
      }
*/

void KLT_TrackingContext::add6Vec(double* out, const double* in, const double coeff) {
  for (int i=0; i<6; ++i)
    out[i] += coeff*in[i];
}
void KLT_TrackingContext::addVec3f(double* out, const Vec3f& v, const double coeff) {
  out[0] += coeff*v.r();
  out[1] += coeff*v.g();
  out[2] += coeff*v.b();
}

// Assumes M is 3x6, G is 2 vectors
void KLT_TrackingContext::calculateM_t(double* M, Vec3f* G, double x, double y) { 
  // row 1
  M[0] = G[0].r()*x; M[1] = G[1].r()*x; M[2] = G[0].r()*y; M[3] = G[1].r()*y;
  M[4] = G[0].r(); M[5] = G[1].r();

  // row 2
  M[6] = G[0].g()*x; M[7] = G[1].g()*x; M[8] = G[0].g()*y; M[9] = G[1].g()*y;
  M[10] = G[0].g(); M[11] = G[1].g();

  // row 3
  M[12] = G[0].b()*x; M[13] = G[1].b()*x; M[14] = G[0].b()*y; M[15] = G[1].b()*y;
  M[16] = G[0].b(); M[17] = G[1].b();

}
// Assumes M is 3x6, G is 2 vectors
void KLT_TrackingContext::addM_t(double* M, Vec3f* G, double x, double y) { 
  // row 1
  M[0] += G[0].r()*x; M[1] += G[1].r()*x; M[2] += G[0].r()*y; M[3] += G[1].r()*y;
  M[4] += G[0].r(); M[5] += G[1].r();

  // row 2
  M[6] += G[0].g()*x; M[7] += G[1].g()*x; M[8] += G[0].g()*y; M[9] += G[1].g()*y;
  M[10] += G[0].g(); M[11] += G[1].g();

  // row 3
  M[12] += G[0].b()*x; M[13] += G[1].b()*x; M[14] += G[0].b()*y; M[15] += G[1].b()*y;
  M[16] += G[0].b(); M[17] += G[1].b();

}

// Computes E += M_t^T V, where M_t is 3x6, V is 3x1
void KLT_TrackingContext::MTranV(double* E, const double* M_t, const Vec3f& V) { 
  int i,y, yhop;
  for (i=0; i<6; i++) {
    for (y=0,yhop=0; y<3; y++, yhop+=6)
      E[i] += M_t[yhop + i] * V.data()[y];
  }
}

void KLT_TrackingContext::calculateStrange(double* strange, const int hh, const int hw) {
  int t,i,j;
  double d;
  for (t=0; t<36; t++) // SPEED
    strange[t] = 0;
  for (j = -hh ; j <= hh ; ++j) 
    for (i = -hw ; i <= hw ; ++i)  {
      d=i*i; strange[0]+=d; strange[7]+=d;                                   // x^2
      d=i*j; strange[2]+=d; strange[9]+=d; strange[12]+=d; strange[19]+=d;  // xy
      strange[4]+=i; strange[11]+=i; strange[24]+=i; strange[31]+=i;         // x
      d=j*j; strange[14]+=d; strange[21]+=d;                                 // y^2
      strange[16]+=j; strange[23]+=j; strange[26]+=j; strange[33]+=j;        // y
      strange[28]+=1.;  strange[35]+=1;                                      // 1
    }
}
/*
void KLT_TrackingContext::calculateStrangeLinear(double* strange, const int hh, const int hw) { // 4x4 // STUB
  
}

void KLT_TrackingContext::calculateStrangeAffine(double* strange, const int hh, const int hw) { // 2x2 // STUB


}
*/
//-----------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------


float KLT_TrackingContext::minEigenvalue(float gxx, float gxy, float gyy) {
    return (gxx + gyy - sqrt((gxx - gyy)*(gxx - gyy) + 4*gxy*gxy))/2.0;
}


void KLT_TrackingContext::transform(const float x, const float y, float* xOut, float* yOut, const float* Z) {
  *xOut = (1.+Z[0])*x + Z[2]      *y + Z[4];
  *yOut = Z[1]      *x + (1.+Z[3])*y + Z[5];
}

void KLT_TrackingContext::transform(const float x, const float y, float* xOut, float* yOut, const double* Z) {
  *xOut = (1.+Z[0])*x + Z[2]      *y + Z[4];
  *yOut = Z[1]      *x + (1.+Z[3])*y + Z[5];
}

/*********************************************************************
 * _computeIntensityDifference
 *
 * Given two images and the window center in both images,
 * aligns the images wrt the window and computes the difference 
 * between the two overlaid images.
 * Aseem
 * Actually, given one window center and one transform Z, blah-blah-blah
 */

void KLT_TrackingContext::computeIntensityDifference(KLT_FloatImage* img1,   /* images */
						     KLT_FloatImage* img2,
						     const int window_width, const int window_height,
						     float x1, float y1,     /* center of window in 1st img */
						     const float* Z,
						     float* imgdiff, int* big_status) {

  int hw = window_width/2, hh = window_height/2;
  float g1, g2;
  int i, j;
  float tx, ty;
  int statusa,statusb;
  
  /* Compute values */
  for (j = -hh ; j <= hh ; j++)
    for (i = -hw ; i <= hw ; i++)  {
      g1 = img1->interpolate(x1+i, y1+j,&statusa);
      transform(i,j,&tx,&ty,Z);
      tx += x1; ty += y1;
      g2 = img2->interpolate(tx,ty, &statusb);
      *imgdiff++ = g1 - g2;
      if (!statusb || !statusa) {
	big_status = 0;
	return;
      }
    }
  *big_status = 1;
}


/*********************************************************************
 * _computeGradientSum
 *
 * Given two gradients and the window center in both images,
 * aligns the gradients wrt the window and computes the sum of the two 
 * overlaid gradients.

 * Aseem
 * Actually, given one window center and one transform Z, blah-blah-blah

 */

void KLT_TrackingContext::computeGradientSum(KLT_FloatImage* gradx1,  /* gradient images */
					     KLT_FloatImage* grady1,
					     KLT_FloatImage* gradx2,
					     KLT_FloatImage* grady2,
					     const int window_width, const int window_height,
					     float x1, float y1,      /* center of window in 1st img */
					     float* Z,
					     float* gradx,      /* output */
					     float* grady, int* big_status) {
  int hw = window_width/2, hh = window_height/2;
  float g1, g2;
  int i, j;
  float tx, ty;
  int statusa,statusb;

  /* Compute values */
  for (j = -hh ; j <= hh ; j++)
    for (i = -hw ; i <= hw ; i++)  {
      transform(i, j,&tx,&ty,Z);
      tx += x1; ty += y1;
      g1 = gradx1->interpolate(x1+i, y1+j,&statusa);
      g2 = gradx2->interpolate(tx, ty,&statusb);
      *gradx++ = g1 + g2;
      if (!statusb || !statusa) {
	big_status = 0;
	return;
      }
      g1 = grady1->interpolate(x1+i, y1+j,&statusa);
      g2 = grady2->interpolate(tx, ty,&statusa);
      *grady++ = g1 + g2;
    }
  *big_status=1;
}


/*********************************************************************
 * _compute6by6GradientMatrix
 *
 */

#define _TA(a,b) floatT[(a)*6 + (b)]  // row, column

void KLT_TrackingContext::compute6by6GradientMatrix(float* gradx, float* grady, const int window_width, 
						    const int window_height, double T[36]) {
  float gx, gy;
  int hw = window_width/2, hh=window_height/2;
  float floatT[36];
  int x,y;
  float gxgx, gygy, gxgy;
    
  /* Compute values */
  memset(floatT,0,36*sizeof(float));
  for (y= -hh; y<=hh; y++) {
    for (x = -hw ; x <= hw; x++)  {
      gx = *gradx++;
      gy = *grady++;
      gxgx = gx*gx; gygy = gy*gy; gxgy = gx*gy;

      // U
      _TA(0,3) += x*y*gxgy;
      _TA(0,0) += x*x*gxgx;
      _TA(1,1) += x*x*gygy;
      _TA(2,2) += y*y*gxgx;
      _TA(3,3) += y*y*gygy;
      _TA(0,1) += x*x*gxgy;
      _TA(0,2) += x*y*gxgx;
      _TA(1,3) += x*y*gygy;
      _TA(2,3) += y*y*gxgy;

      // V~
      _TA(4,0) += x*gxgx;
      _TA(4,1) += x*gxgy;
      _TA(4,2) += y*gxgx;
      _TA(4,3) += y*gxgy;
      _TA(5,1) += x*gygy;
      _TA(5,3) += y*gygy;

      // Z
      _TA(4,4) += gxgx;
      _TA(5,5) += gygy;
      _TA(5,4) += gxgy;
    }
  }

  // U
  _TA(3,0) = _TA(2,1) = _TA(1,2) = _TA(0,3);
  _TA(1,0) = _TA(0,1);
  _TA(2,0) = _TA(0,2);
  _TA(3,1) = _TA(1,3);
  _TA(3,2) = _TA(2,3);

  // V~
  _TA(5,0) = _TA(4,1);
  _TA(5,2) = _TA(4,3);

  // Z
  _TA(4,5) = _TA(5,4);

  // V
  for (y=4; y<=5; y++)
    for (x=0; x<=3; x++)
      _TA(x,y) = _TA(y,x);

  for (x=0; x<36; x++)
    T[x] = floatT[x];
  
}

	
/*********************************************************************
 * _compute2by1ErrorVector
 *
 */

void KLT_TrackingContext::compute6by1ErrorVector(float* imgdiff, float* gradx,float* grady, 
						 const int window_width, const int window_height, double* a) {
  float diff,dgx,dgy;
  int x,y;
  float floata[6];
  int hw = window_width/2, hh=window_height/2;
  
  floata[0] = floata[1] = floata[2] = floata[3] = floata[4] = floata[5] = 0;
  /* Compute values */
  for (y= -hh; y<=hh; y++) {
    for (x = -hw ; x <= hw; x++)  {
      diff = *imgdiff++;
      dgx = diff * (*gradx++);
      dgy = diff * (*grady++);
      floata[4] += dgx;
      floata[5] += dgy;
      floata[0] += x*dgx;
      floata[1] += x*dgy;
      floata[2] += y*dgx;
      floata[3] += y*dgy;
    }
  }

  for (x=0; x<6; x++)
    a[x]=floata[x];
}



// DEBUG
/*
void KLT_TrackingContext::giveArrays(double* A, double *b, int n) {
  if (_A == NULL) {
    _A = A;
    _b = b;
    _n = n;
  }
  else {
    int i,j;
    assert(_n == n);
    printf("A diffs: \n");
    for (j=0; j<n; j++) {
      for (i=0; i<n; i++)
	printf("%.2f ",_A[j*n+i] - A[j*n+i]);
      printf("\n");
    }
    
    printf("first A: \n");
    for (j=0; j<n; j++) {
      for (i=0; i<n; i++)
	printf("%.2f ",_A[j*n+i]);
      printf("\n");
    }

    printf("second A: \n");
    for (j=0; j<n; j++) {
      for (i=0; i<n; i++)
	printf("%.2f ",A[j*n+i]);
      printf("\n");
    }

    
    printf("b diffs: \n");
    for (i=0; i<n; i++)
      printf("%.2f ", b[i] - _b[i]);
    printf("\n");
    
    exit(0);
  }
}
*/


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



#ifndef _KLT_H_
#define _KLT_H_

using namespace std;
#include <stdio.h>
#include <qimage.h>
#include <qdatetime.h>
#include <float.h>
#include <qthread.h>
#include "error.h"
#include "../jl_vectors.h"
#include "mysparsemat.h"
#include "keeper.h"
#include "multiTrackData.h"
#include "multiSplineData.h"
#include "multiKeeper.h"
#include "splineKeeper.h"
#include "buildingSplineKeeper.h"
#include "obsCache.h"


#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "myassert.h"

typedef SplineKeeper CSplineKeeper;

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#define KLT_TRACKED           0
#define KLT_NOT_FOUND        -1
#define KLT_SMALL_DET        -2
#define KLT_MAX_ITERATIONS   -3
#define KLT_OOB              -4
#define KLT_LARGE_RESIDUE    -5

#include "pyramid.h"

enum KLT_ThreadTask {LCT2, SPLINE};

class KLT_TrackingContext;

class KLT_FullPyramid {  // greyscale version

 public:

  KLT_FullPyramid();
  void initMe(const QImage im, const KLT_TrackingContext* tc);
  void initMeFromEdges(const QImage im, const KLT_TrackingContext* tc);
  KLT_FullPyramid(const QImage im, const KLT_TrackingContext* tc);
  ~KLT_FullPyramid();

  void write(FILE* fp) const;
  bool load(FILE* fp, const KLT_TrackingContext* tc); // will return false if nlevels not same for tc
  void writeImages(char* imgname, char* gxname, char* gyname);

  KLT_Pyramid* img;
  KLT_Pyramid* gradx;
  KLT_Pyramid* grady;
  int nPyramidLevels;
};



class KLT_FullCPyramid {  // color version

 public:

  KLT_FullCPyramid();
  void initMe(const QImage im, const KLT_TrackingContext* tc);
  KLT_FullCPyramid(const QImage im, const KLT_TrackingContext* tc);
  ~KLT_FullCPyramid();

  void write(FILE* fp) const;
  bool load(FILE* fp, const KLT_TrackingContext* tc); // will return false if nlevels not same for tc
  void writeImages(char* imgname, char* gxname, char* gyname);

  KLT_ColorPyramid* img;
  KLT_ColorPyramid* gradx;
  KLT_ColorPyramid* grady;
  int nPyramidLevels;
};

void printDoubleArray(FILE* fp, const double* a, const int nrows, const int ncols);

class KLT_TrackingContext : public QThread {
  
 public:


  KLT_TrackingContext();
  #include "kltSpline.h"
  //KLT_TrackingContext(KLT_FullPyramid* pyr);

  //void calcAPyramids(QImage im);

  //int track(const KLT_FullPyramid *aPyrm,const KLT_FullPyramid *bPrym, 
  //    int window_width, int window_height,
  //    const float x, const float y, float *Z);
  // aPyrm: first frame pyramid
  // bPyrm: new frame pyramid
  // (x,y): location in current frame
  // w,h: window width & height
  // Z[4],Z[5]: estimate for new location in new frame (DELTA from currloc, not actual location)
  // add [1 0 x, 0 1 y] to get actual location
  // returns 1 for success, 0 out of bounds, 2 bad track
  //int colorTrack(const KLT_FullCPyramid *aPyrm,const KLT_FullCPyramid *bPyrm, 
  //	 int window_width, int window_height,
  //	 const float x, const float y, double *Z);
  // Z: should be 12, 0-5 is a's transform delta and 6-11 for b

  //int longPatchTrack(const KLT_FullCPyramid** pyrms, int numFrames, float ax, float ay, int window_width, 
  //	     int window_height, double* param);

  // pyrms: pyramids frames a->b inclusive
  // numFrames: b-a, so number of frames not including first
  // ax,ay: starting location in frame a
  // width, height: you know
  // param: float[numFrames*paramLength] of motion parameters.  param[0] refers to a+1 image
  // //paramLength: number of parameters, one image.  Now 6.
  void setupLongCurveTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, MultiTrackData* mt, bool redo);

  virtual void run();  // thread handler
  void runNoThread();  // thread handler

  //int curveTrack(const KLT_FullCPyramid *aPyrm,const KLT_FullCPyramid *bPyrm, 
  //	 const int numP, Vec2f* Plocs, int wwidth, int wlh, int whh, Vec2f* Z);
  // Plocs: numP vectors
  // Z: numP delta vectors
  int longCurveTrack(const KLT_FullCPyramid** pyrms, const int numFrames, const int numP,
		     const int wwidth, const int wlh, const int whh, Vec2f* Plocs, Vec2f* Z);
  int longCurveReshape(int numFrames, int numP, Vec2f* Plocs, 
		       Vec2f* Z, Vec2f* Sb, bool pL);


  void copySettings(const KLT_TrackingContext *o);

  ~KLT_TrackingContext();


  /* Available to user */


  /* Available, but hopefully can ignore */
  float min_displacement;	/* th for stopping tracking when pixel changes little */
  int max_iterations;		/* th for stopping tracking when too many iterations */
  float max_residue;		/* th for stopping tracking when residue is large */
  float grad_sigma;
  float smooth_sigma_fact;
  float pyramid_sigma_fact;
  int nPyramidLevels;
  int subsampling;	
  double smoothAlpha, edgeWeight;
  double smooth2Deriv, smooth1Deriv, smooth0Deriv, shape2Deriv;
  bool useImage;
  bool useEdges;
  bool useDiffScale;
  bool constantWindow;
  bool pinLast;
  bool useNormalEq;
  bool printSingulars;
  bool usePseudo;
  bool dumpWindows;
  int D2mode;
  bool _stateOk;

  KLT_ThreadTask _ttask;
  const KLT_FullCPyramid  **__pyrms; const KLT_FullPyramid **__pyrmsE;
  MultiTrackData* _mt; bool __redo;

  int steinum; //G!

 private:

  void longCurveTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, bool redo);
  double robust(const double r) const;
  //void initCrossDerivs(const KLT_FullCPyramid** pyrms, MultiTrackData* mt, const int level, const double invsubs);
  void initEdgeMins(const KLT_FullPyramid** pyrmsE, 
				       const int level, const double invsubs);
  double createMatrices2(const KLT_FullCPyramid** pyrms, const int numFrames, const int l,
			  const int wwidth, const int wlh, const int whh, const Vec2f* Plocs, Vec2f* Z,const int level, 
			  const double invsubs, Keeper* keep); 
  double createMatrices3(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, Vec2f* Z, 
			 const int level, const double invsubs, MultiKeeper* keep,MultiKeeper *imgKeep, MultiKeeper *edgeKeep);
  void calculateNewK(double* K, const double i, const double j, const double invsubs) const;
  void hnew(const int n, const float i, const float j, const float invc, const Vec2f* P, Vec2f& result) const;
  void hnew2(const int cp, const int c, const int t, const float i, const float j, const float invsubs, 
	     const Vec2f* Z, Vec2f& result) const;
  int steihaugSolver(GenKeeper* im, double* x, const double trustRadius, bool* boundaryHit, double* error2ptr = NULL); // returns # iterations
  double createTestSolution(const Vec2f* Z, const double* sol, const int numVec, const int l, Vec2f* Z2);
  void projectToTR(double *x, const double* d, const double* p, double trustRadius, int n);
  void projectToTR2(double *x, const double* d, const double* p, double trustRadius, int n);
  bool internalLongCurveTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, 
			      const int level, const double invsubs, const int curveSampling, const int hold, const int maxIterations);
  void safeTrack2(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE,
		  const int level, const double invsubs, const int curveSampling, const int hold, 
		  const int maxIterations);

  void checkWindow(int& window_width, int& window_height);
  void mutualInit();
  //int internalTrackFeature(float x1, float y1, float* Z, const int window_width, const int window_height,
  //		   const KLT_FullPyramid *aPyrm, const KLT_FullPyramid *bPyrm, const int level);
  float minEigenvalue(float gxx, float gxy, float gyy); 
  void transform(const float x, const float y, float* xOut, float* yOut, const float* Z);
  void transform(const float x, const float y, float* xOut, float* yOut, const double* Z);
  void computeIntensityDifference(KLT_FloatImage* img1, KLT_FloatImage* img2, const int window_width, 
				  const int window_height, float x1, float y1,
				  const float* Z,float* imgdiff, int* big_status);
  void computeGradientSum(KLT_FloatImage* gradx1, KLT_FloatImage* grady1, KLT_FloatImage* gradx2, 
			  KLT_FloatImage* grady2, const int window_width, const int window_height, 
			  float x1, float y1, float* Z, float* gradx, float* grady, int* big_status);
  void compute6by6GradientMatrix(float* gradx, float* grady, const int window_width, 
				 const int window_height, double T[36]);
  void compute6by1ErrorVector(float* imgdiff, float* gradx, float* grady, const int window_width, 
			      const int window_height, double* a);
  
  //int internalLongPatchTrack(const KLT_FullCPyramid** pyrms, int numFrames, float ax, float ay, 
  //int window_width, int window_height, double* param, int level);
  void MTranV(double* E, const double* M_t, const Vec3f& V); // E += M_t^T V
  //void MTranM(double* Z, const int span, const double coeff,double* MTran,double* M);
  // Z += coeff Mtran^T M
  // Z: array to put result in
  // span: actual width of Z (separates spans)
  // coeff: pre-multiplies all result values
  // MTran: array to transpose in M_T M
  // M: array not to transpose in M_T M
  // currently assumes Mtran, M are 3x6
  void calculateM_t(double* M, Vec3f* G, double x, double y);
  void add6by6(double* Z, int offset, int span, double coeff, double* adder);
  void add6Vec(double* out, const double* in, const double coeff);
  int internalColorTrackFeature(float x1, float y1, double* Z,
				const int window_width, const int window_height,
				const KLT_FullCPyramid *aPyrm,
				const KLT_FullCPyramid *bPyrm, 
				const int level);  
  //int internalLongPatchTrack2(const KLT_FullCPyramid** pyrms, int numFrames, float ax, float ay, 
  //		      int window_width, int window_height, 
  //		      double* param, const int level);

  void addnbyk(double* Z, const int span, const double coeff, const double* adder, const int n, const int k);
  void addnbyk(double* Z, const int span, const double* adder, const int n, const int k);
  void addVec3f(double* out, const Vec3f& v, const double coeff);
  void addM_t(double* M, Vec3f* G, double x, double y);
  //int pseudoSolve(double* x, double* Z, double* e, const int n); // returns number of vars below threshold

  void calculateStrange(double* strange, const int hh, const int hw); // 6x6
  //void calculateStrangeLinear(double* strange, const int hh, const int hw); // 4x4
  //void calculateStrangeAffine(double* strange, const int hh, const int hw); // 2x2

  //int internalCurveTrack(const KLT_FullCPyramid *aPyrm,const KLT_FullCPyramid *bPyrm, 
  //		 const int numP, const Vec2f* Plocs, int wwidth, int wlh, int whh, 
  //		 Vec2f* Z, int level);
  //int calculateMt(const int n, const int i, const int j, const int l, const double dinvc,
  //	  double *G, double *KS, double *MS, double *MtS);
  int calculateMtNorm(const int n, const int i, const int j, const int l, const Vec2f* P,
		      double *G, double *KS, double *MS); 
  void calculateK(double* K, const int n, const double i, const double j, const double c, const int l) const;
  void calculateKStrip(double* K, const double i, const double j, const double invc) const;
  void calculateKStrip0(double* K, const double i, const double j, const double invc) const;
  void calculateKStripl(double* K, const double i, const double j, const double invc) const;
  void calculateKStripNorm(double* K, const int n, const double i, const double j, const Vec2f* P) const;
  void calculateKStripNorm0(double* K, const double i, const double j, const Vec2f* P) const;
  void calculateKStripNorml(double* K, const int n, const double i, const double j, const Vec2f* P) const;
  void h(const int n, const float i, const float j, const int l, const float c,
	 const Vec2f* P, Vec2f& result) const;
  void hnorm(const int n, const float i, const float j, const int l, const Vec2f* P, Vec2f& result) const;
  void rmmultplus(double *rm,double *a,double *b,int n,int m,int l) const;
  void rmmult(double *mat,double *a,double *b,int m,int k,int n) const;
  int internalLongCurveTrack(const KLT_FullCPyramid** pyrms,const int numFrames, const int l,
			     const int wwidth, const int wlh, const int whh, const Vec2f* Plocs, Vec2f* Z, 
			     const int level, const double invsubs, const int curveSampling);
  int internalSetulbTrack(const KLT_FullCPyramid** pyrms,const int numFrames, const int l,
			     const int wwidth, const int wlh, const int whh, const Vec2f* Plocs, Vec2f* Z, 
			     const int level, const float spacing);
  double evaluateObj(double* objJac, const KLT_FullCPyramid** pyrms,const int numFrames, const int l,
					const int wwidth, const int wlh, const int whh, const Vec2f* Plocs, const double* Zd,
					const int level, const float spacing);
  void visMatrix(const double *A, const int w, const int h) const;
  void assertSym(const double *A, const int cols) const;
  // DEBUG
  //void giveArrays(double* A, double *b, int n);

  //double *_A, *_b;
  //int _n;

  double mtmtran[18], mtmresult[36]; // prevents PARALLEL
  mutable double _mstore[50];
};
    

inline void KLT_TrackingContext::rmmult(double *rm,double *a,double *b,int n,int m,int l) const { 
  double z,*q0,*p,*q; int i,j,k;
  q0 = _mstore;
  //q0=(double *)calloc(m,sizeof(double));
  for(i=0; i<l ;++i,++rm){
    for(k=0,p=b+i; k<m ;p+=l) q0[k++]= *p;
    for(j=0,p=a,q=rm; j<n ;++j,q+=l){
      for(k=0,z=0.; k<m ;) z+= *p++ * q0[k++];
      *q=z;
     }
   }
  //free(q0);
}


inline void KLT_TrackingContext::rmmultplus(double *rm,double *a,double *b,int n,int m,int l) const { 
  double z,*q0,*p, *q; int i,j,k;
  q0 = _mstore;
  //q0=(double *)calloc(m,sizeof(double));
  for(i=0; i<l ;++i,++rm){
    for(k=0,p=b+i; k<m ;p+=l) q0[k++]= *p;
    for(j=0,p=a,q=rm; j<n ;++j,q+=l){
      for(k=0,z=0.; k<m ;) z+= *p++ * q0[k++];
      *q+=z;
    }
  }
  //free(q0);
}


inline void KLT_TrackingContext::addnbyk(double* Z, const int span, const double coeff, 
				  const double* adder, const int n, const int k) {
  int i,j, shop, khop;
  for (j=0, shop=0, khop=0; j<n; ++j, shop+=span, khop+=k) 
    for (i=0; i<k; ++i)
      Z[shop + i] += coeff*adder[khop+i];
}

inline void KLT_TrackingContext::addnbyk(double* Z, const int span,  
				  const double* adder, const int n, const int k) {
  int i,j, shop, khop;
  for (j=0, shop=0, khop=0; j<n; ++j, shop+=span, khop+=k) 
    for (i=0; i<k; ++i)
      Z[shop + i] += adder[khop+i];
}


/*******************
 * Functions
 */


/* Aseem 
void KLTCreateFile(KLT_PixelType *img, int ncols, int nrows,
		   KLT_TrackingContext tc, FILE* fp1, FILE* fp2, FILE* fp3);
*/



/* Processing 


int KLTPreTrackAFeature(KLT_TrackingContext tc,  _KLT_Pyramid b, _KLT_Pyramid bgx,
			 _KLT_Pyramid bgy, float x, float y,
			 float *xOut, float *yOut);
*/




/* Writing/Reading 

void KLT_StorePyramids(KLT_TrackingContext tc, FILE* fp1, FILE* fp2, FILE* fp3);
*/

#endif







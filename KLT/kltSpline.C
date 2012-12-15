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



#include "klt.h"
#include "myassert.h"

// BIG STUB (S)
// 1. Edges
         // 2. HandleCOnstraints
// 3. Preconditioning
         // 4. Propagating uBlas vectors higher

//#define _TIME_

void printSparseVector(const Ubcv& J) {
  Ubcv::const_iterator i = J.begin();
  for ( ; i!=J.end(); ++i)
    printf("(%d, %.5e), ", i.index(), *i);
  printf("\n");
}

void KLT_TrackingContext::setupSplineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, 
					   MultiSplineData* mts, bool redo) {

   _ttask = SPLINE;
  __pyrms = pyrms; __pyrmsE = pyrmsE;
  _mts = mts; __redo = redo;


  _useD1 = _useD2 = _useD0 = true;
}


#define DELETE_STK 
void KLT_TrackingContext::splineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, 
					 bool redo) {

  if (_mts->_numFrames <2) {
    printf("Too short a span\n");
    DELETE_STK;
    return;
  }

  QTime totaltime;
  totaltime.start();
  tt1 = tt2 = tt3 = tt4 = tt5 = tt6 = tt7 = tt8 = tt9 = tt10 = tt11 = 0;
  //fpo = fopen("iterations.txt","w"); 

  if (redo) { 
    // create smooth initial position
    /*bool temp1 = useImage, temp2 = useEdges;
    useImage = false;  
    useEdges = false;
    innerSplineTrack(pyrms, pyrmsE, 0);
    useImage = temp1; useEdges = temp2;*/
    _mts->takeControls(&(_mts->_Z));
    _mts->discretizeAll(2, RESAMPLE_CONSISTENTLY); // comment out if RESAMPLE wanted    
    doSplineShapeInterp();
  }
  
  //innerSplineTrack(pyrms, pyrmsE, 0);  
  printf("Starting image term\n");
  int mylevels = nPyramidLevels;
  _mts->takeControls(&(_mts->_Z));
  _mts->discretizeAll(2, RESAMPLE_CONSISTENTLY); // comment out if RESAMPLE wanted     2 sampling not usual

  for (int r = mylevels - 1 ; r >= 0 ; r--)  {   
    if (useEdges) {
      initSplineEdgeMins(pyrmsE,r);  
    }
    printf("\nLevel %d\n",r);

    safeSplineTrack(pyrms, pyrmsE, r); 
    
  }
  //fclose(fpo); 

  _mts->takeControls(&(_mts->_Z));
  _mts->discretizeAll(2, RESAMPLE_CONSISTENTLY); // comment out if RESAMPLE wanted    

  printf("Total elapsed time %d\n",totaltime.elapsed());
  printf("time1: %d\n",tt1);
#ifdef _TIME_
  printf("time2: %d\n",tt2);
  printf("time3: %d\n",tt3);
  printf("time4: %d\n",tt4);
  printf("time5: %d\n",tt5);
  printf("time6: %d\n",tt6);
  printf("time7: %d\n",tt7);
  printf("time8: %d\n",tt8);
  printf("time9: %d\n",tt9);
  printf("time10: %d\n",tt10);
  printf("time11: %d\n",tt11);
  #endif

  DELETE_STK;
}


void KLT_TrackingContext::doOneStep(CSplineKeeper* keep, double* x) {
  keep->refresh();
  createSplineMatrices(NULL, NULL, &(_mts->_Z), 0, keep);
  //keep->outputMat("A.dat"); std::exit(0); 
  memset(x,0,sizeof(double)*keep->numVar());  
  int steps=0;
  double error2 = .0000001; //.00001; // .1

  bool bhit;
  steps = steihaugSolver(keep,x,DBL_MAX, &bhit, &error2);


  //error2 = ConjGrad(keep->numVar(), keep, x,keep->neg_g(), .000001, &steps);
  //printf("problem Solved by ConjGrad in %d\n with error %.4e\n",steps, error2);

  _mts->_z_mutex->lock(); 
  ZVec Z2 = _mts->_Z;
  _mts->createTestSol(&Z2, x, &(_mts->_Z));
  _mts->_z_mutex->unlock(); 

  
  /*FILE* fp = fopen("x.txt", "w"); 
  for (int i=0; i<keep->numVar(); ++i) {
    fprintf(fp, "%.10e\n",x[i]); 
  }
  fclose(fp);  std::exit(0);*/
}

void KLT_TrackingContext::doSplineShapeInterp() {

  bool temp1 = useImage, temp2 = useEdges;
  useImage = false;  
  useEdges = false;
  
  double* x = new double[_mts->numVars()];
  CSplineKeeper *keep = new CSplineKeeper(_mts);
  ZVec Z2;

  // first: D2 + D1
  _useD1 = false; 
  //_useD2 = false;  //G!
  doOneStep(keep, x); 

  //_mts->_z_mutex->lock(); _mts->createTestSol(&(_mts->_Z), x, &(_mts->_Z)); _mts->_z_mutex->unlock();  
  //  keep->outputMat("A0.dat");
  //keep->outputg("g0.txt");

  /*
  Z2 = _mts->_Z;
  _mts->_z_mutex->lock(); _mts->createTestSol(&(_mts->_Z), x, &(_mts->_Z)); _mts->_z_mutex->unlock();  
  printf("evaluate\n");
  keep->refresh(); 
  createSplineMatrices(NULL, NULL, &(_mts->_Z), 0, keep);  
  _mts->_Z = Z2;

  // next, D1 then D2+D1
  _useD2 = false;
  doOneStep(keep, x);
  _mts->_z_mutex->lock(); _mts->createTestSol(&(_mts->_Z), x, &(_mts->_Z)); _mts->_z_mutex->unlock();  
  _mts->outputControlVars("P1.txt");
  _useD2 = true;
  doOneStep(keep, x);
  _mts->_z_mutex->lock(); _mts->createTestSol(&(_mts->_Z), x, &(_mts->_Z)); _mts->_z_mutex->unlock();  
  keep->outputMat("A1.dat");
  keep->outputg("g1.txt");
  */
  //keep->shit(); 
  //keep->outputMat("A.dat"); 
    
  
  
  //_useD2 = true; 
  //doOneStep(keep, x);

  _useD2 = true; 
  _useD1 = true;   
  printf("Starting D1 shape interp\n");
  //safeSplineTrack(NULL, NULL, 0);     
  //keep->refresh();  
  //createSplineMatrices(NULL, NULL, &(_mts->_Z), 0, keep);   

  useImage = temp1; useEdges = temp2;


  delete[] x; 
  delete keep; 
}


void KLT_TrackingContext::safeSplineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, const int level) {
  bool res=true;
  while (1) {
    res = innerSplineTrack(pyrms, pyrmsE, level);

    if (res)
      break;
    else {
      _mts->_z_wait->wait();  // make sure dragging is done
      printf("Retry\n");
      _mts->takeControls(&(_mts->_Z));
      _mts->discretizeAll(2, RESAMPLE_CONSISTENTLY); // comment out if RESAMPLE wanted    
      _stateOk = true;
    }
  }

}


#define IST2_DELETE   delete[] x; delete keep1; delete keep2; //if (imgKeep) delete imgKeep; if (edgeKeep) delete edgeKeep;
bool KLT_TrackingContext::innerSplineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, 
					   const int level) {
    int iteration=0;
  
  bool toContinue = true;
  //int numVec = _mt->_numFrames-1;
  //int sizeZ = _mts->size_Z();

  //Vec2f* Z2 = new Vec2f[sizeZ];
  //memcpy(Z2, _mts->_Z, sizeZ*sizeof(Vec2f));
  ZVec Z2 = _mts->_Z;
  double trustRadius=10.;  // initial trust radius.   
  CSplineKeeper *keep1 = new CSplineKeeper(_mts), *keep2 = new CSplineKeeper(_mts);  
  /*CSplineKeeper *imgKeep = NULL, *edgeKeep = NULL;
  if (useImage)
    imgKeep = new CSplineKeeper(keep1);
  if (useEdges)
  edgeKeep = new CSplineKeeper(keep1);*/
  double *x = new double[keep1->numVar()];


  double prec;
  if (level==0)
    prec = .1;  // .25
  else
    prec = 0.5;  // 1
  double currTheta, newTheta;
  //keep1->printCurves(mt->_Z);
  currTheta = createSplineMatrices(pyrms, pyrmsE, &(_mts->_Z), level, keep1/*, imgKeep, edgeKeep*/);
  //keep1->outputMat("B.dat"); std::exit(0); 
  
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
  

  int maxIterations = 1000;
  do {
    fprintf(stdout,"starting iteration %d\n", iteration);
    
    double ro=0, maxStep = 10.; // 10 is just to force into loop initially
    
    bool boundaryHit;
    while (maxStep > .1  && trustRadius >= prec && iteration<=maxIterations) {
      
      /*memset(x,0,sizeof(double)*keep1->numVar());  
      int steps=0;
      ConjGrad(keep1->numVar(), keep1, x,keep1->g(), .000001, &steps);
      printf("Solved by COnjGrad in %d\n",steps);
      memset(x,0,sizeof(double)*keep1->numVar());*/

#ifdef _TIME_
      if (tt8 == 0) time8.start();
      else time8.restart();
#endif

      int numIter = steihaugSolver(keep1, x, trustRadius, &boundaryHit); // stei, remember to clear x

      //fprintf(fpo,"level: %d iter: %3d numIter: %4d\n",level,iteration,numIter); 
      //fflush(fpo);

#ifdef _TIME_
      tt8 += time8.elapsed();
#endif
      

      /*
      if (!res) {
	FILE* fp = fopen("/homes/gws/aseem/A.dat","w"); 
	keep1->printSparse(fp);
	keep1->outputg();
	fclose(fp);
	std::exit(0);
      }
      */

      /*
      FILE* fp = fopen("x.txt","w");
      for (int rt=0; rt<keep1->numVar(); rt++) 
	fprintf(fp,"%.10e\n",x[rt]);
      fclose(fp);
      std::exit(0);*/
      //exit(0);

      // eval solution
      if (_stateOk) {
	maxStep = _mts->createTestSol(&(_mts->_Z), x, &Z2); // write to Z2
	printf("max step %f\n",maxStep);
	assert(maxStep < trustRadius + .00001);
      }
      
      if (_stateOk)
	newTheta = createSplineMatrices(pyrms, pyrmsE, &Z2, level, keep2/*, imgKeep, edgeKeep*/);
      else
	newTheta = DBL_MAX;
      
      if (_stateOk) {
	ro = keep1->calculateRo(x, newTheta, currTheta);
	printf("ro: %f\n",ro);
	
	printf("old %f, new %f\n", currTheta, newTheta); 
	if (ro < .25) {// .25 otherwise , or 0
	  trustRadius *= .25;
	  printf("reducing trustRadius to %f\n",trustRadius);
	}
	else if (ro > .75 &&  boundaryHit) // add ro > .75 &&
	  trustRadius = MIN(2.*trustRadius, 10.);
      }

      _mts->_z_mutex->lock(); 
      if (_stateOk) {      
	if (ro <= 0) { // don't take step
	  printf("Not taking step\n");
	  keep2->refresh();
	}
	else { // step is fine
	  iteration++;
	  printf("Taking step, iteration %d\n\n", iteration);
	  CSplineKeeper* kswap = keep1;
	  keep1 = keep2; 
	  keep2 = kswap;
	  keep2->refresh();
	  currTheta = newTheta;
	  
	  _mts->_Z = Z2;  // only need to copy data that can change (not doing this currently)	
	  //memcpy(_mts->_Z,Z2,sizeZ*sizeof(Vec2f)); 
	}
	_mts->_z_mutex->unlock();
      }
      else {  // don't apply update, get outta here
	_mts->_z_mutex->unlock();
	IST2_DELETE;
	return false;
      } // end else, stateOk block

    } // while loop
    
    if (trustRadius < prec || maxStep < .1  || iteration>maxIterations)
      toContinue=false;
    else
      toContinue=true;
    printf("\n\n");
    
    
  } while (toContinue);

  if (iteration == maxIterations)
    printf("Hit max number of iterations\n");
  
  IST2_DELETE;
  return true;

}


bool fequal(const double a, const double b) {
  return (fabs(a-b) < .0001);
}

#define DELETE_CSM  delete[] theta; delete[] thetas; delete[] thetas1;  delete[] thetas0; delete[] thetaE; delete imgKeep; delete edgeKeep; delete D0Keep; delete D1Keep;

double KLT_TrackingContext::createSplineMatrices(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE,
						 ZVec* Z, const int level, CSplineKeeper* keep
						 /*,CSplineKeeper *imgKeep, CSplineKeeper *edgeKeep*/
						 ) {
  if (tt1 == 0) time1.start();
  else time1.restart();
#ifdef _TIME_
  if (tt9 == 0) time9.start();
  else time9.restart();
#endif

  int numFrames = _mts->_numFrames;
  int n,c,j, validPoint, di, t;
  Vec2f loc;
  Vec3f col1, col, dummy; 
  double K[16], K1[16], grad[24], grad1[24];
  double G_t[6],  G_t1[6]; // 3x2 Jacobian matrices
  int vars[4], vars1[4]; // full-2, but 1 per pair
  int ul = numFrames-1, un, un1;
  TrackSample ds,ds1;
  int ds1num;
  float invsubs = 1./float(pow2(level));
  float subs = float(pow2(level));
  double invcw, invew, invd0, invd1, invd2;
  //SamplePreComp snum, snum1;

  // for sample-based schemes
  int vars00[4], *vars01=NULL, vars02[4], vars10[4], *vars11=NULL, vars12[4];
  bool D2_parameter_based = false;
  bool D1_parameter_based = false;
  bool checkOk;

  double *theta = new double[numFrames], *thetaE = new double[numFrames], *thetas = new double[numFrames], 
    *thetas1 = new double[numFrames], *thetas0 = new double[numFrames]; 
  memset(theta,0,numFrames*sizeof(double));  // image term
  memset(thetas,0,numFrames*sizeof(double)); //2nd deriv term
  memset(thetas1,0,numFrames*sizeof(double));
  memset(thetas0,0,numFrames*sizeof(double));
  memset(thetaE,0,numFrames*sizeof(double));
  
  int imgCompCount=0, edgeCompCount=0, smooth0Count = 0, smooth1Count=0, smooth2Count=0;
  int stride = pow2(level);  // careful
  /*if (imgKeep)
    imgKeep->refresh(); 
  if (edgeKeep)
  edgeKeep->refresh();*/
  //printf("Here 1\n");
  CSplineKeeper *imgKeep = new CSplineKeeper(keep), *edgeKeep = new CSplineKeeper(keep),
    *D0Keep = new CSplineKeeper(keep), *D1Keep = new CSplineKeeper(keep); // D2Keep is just keep, which is already refreshed

  ObsCache imgCache(imgKeep), D0Cache(D0Keep), D1Cache(D1Keep), D2Cache(keep), ECache(edgeKeep);

  _mts->takeControls(Z);
  _mts->discretizeAll(2, REEVALUATE, true); // or RESAMPLE

#ifdef _TIME_
  tt9 += time9.elapsed();
#endif

  for (c=0; c < _mts->_nCurves && _stateOk; c++) { // iterate over curves

    
    for (t=0; t<_mts->_numFrames && _stateOk; ++t) {    // ONE,_numFrames-1,  -1, shouldn't be there

      //_mts->discretize(t, c, pow2(level));  

      un = _mts->getNumSamples(t,c)-1;
      un1 = _mts->getNumSamples(t+1, c) - 1;
      assert(!useEdges || !_mts->useEdges(c) || _mts->_edgeMins[c].size() == un);
      for (n=0; n<=un; ++n) {// iterate over samples of curve  

#ifdef _TIME_
	if (tt5 == 0) time5.start();
	else time5.restart();
#endif

	if (t==0)
	  _mts->getDiscreteSample(t,c,n,&ds);
	else
	  _mts->getDiscreteSample(t,c,n,&ds, vars);
	
	//snum.compute(ds);
	
	if (t==ul)
	  ds1num = _mts->getExistingCorrDiscreteSample(t,c,n,&ds1);
	else
	  ds1num = _mts->getExistingCorrDiscreteSample(t,c,n,&ds1, vars1);

	// get some common samples
	TrackSample ds00,ds02,ds12;
	
#ifdef _TIME_
	if (tt5 == 0) time5.start();
	else time5.restart();
#endif
	
	if (t>0) {
	  if (n!=un) _mts->getDiscreteSample(t,  c, n+1,      &ds02,  vars02);
	  vars01 = vars;
	}
	else
	  if (n!=un) _mts->getDiscreteSample(t,  c, n+1,      &ds02);
	TrackSample& ds01 = ds;  
	
	if (t<ul) {
	  //if (ds1num!=un1) _mts->getDiscreteSample(t+1,c, ds1num+1, &ds12, vars12);
	  if (ds1num!=un1 && n!=un) _mts->getExistingCorrDiscreteSample(t,c, n+1, &ds12, vars12);
	  vars11 = vars1;
	}
	else
	  //if (ds1num!=un1) _mts->getDiscreteSample(t+1,c, ds1num+1, &ds12);
	  if (ds1num!=un1 && n!=un) _mts->getExistingCorrDiscreteSample(t,c, n+1, &ds12);
	TrackSample& ds11 = ds1;  
	
#ifdef _TIME_
	tt5 += time5.elapsed();
#endif


#ifdef _TIME_
	tt5 += time5.elapsed();
#endif

	//snum1.compute(ds1);
	if (useImage && n%stride==0) {

	  // ONLY FOR NON_NORMALIZED METHOD!
	  if (n==un || ds1num == un1) 
	    continue; 

	  for (j = _mts->_trackWidths[c].x() ; j <= _mts->_trackWidths[c].y() ; ++j) { 

#ifdef _TIME_
	    if (tt7 == 0) time7.start();
	    else time7.restart();
#endif

	    validPoint = 0;
	    splineLoc(ds1, j, invsubs, &loc);
	    //printf("n: %d aloc: %f %f ",n,loc.x(), loc.y());
	    di = pyrms[t+1]->img->color(loc.x(),loc.y(),level, col1);
	    validPoint = MIN( di, validPoint); 
	    pyrms[t+1]->gradx->color(loc.x(),loc.y(), level, dummy); 
	    G_t1[0]=dummy.r(); G_t1[2]=dummy.g(); G_t1[4]=dummy.b();
	    pyrms[t+1]->grady->color(loc.x(),loc.y(), level, dummy);
	    G_t1[1]=dummy.r(); G_t1[3]=dummy.g(); G_t1[5]=dummy.b();

	    splineLoc(ds,j,invsubs,&loc);
	    //printf("bloc: %f %f\n",loc.x(), loc.y());
	    di = pyrms[t]->img->color(loc.x(),loc.y(),level, col);
	    validPoint = MIN( di, validPoint); 
	    pyrms[t]->gradx->color(loc.x(),loc.y(), level, dummy); 
	    G_t[0]=dummy.r(); G_t[2]=dummy.g(); G_t[4]=dummy.b();
	    pyrms[t]->grady->color(loc.x(),loc.y(), level, dummy);
	    G_t[1]=dummy.r(); G_t[3]=dummy.g(); G_t[5]=dummy.b();

	    if (validPoint == 0) {
	      noSubSplineLoc(ds,j,subs,&loc);
	      validPoint = _mts->occluded(ds._loc, t);
	      if (validPoint == 0) {
		noSubSplineLoc(ds1,j,subs,&loc);
		validPoint = _mts->occluded(ds1._loc, t+1);
	      }
	      if (validPoint == 1)
		printf("Occluded\n");
	    }
	  
	    // THis can happen if two adjacent samples fall on top of each other
	    if (validPoint==0 && (!finite(ds._k) || isnan(ds._k) ||
				  !finite(ds1._k) || isnan(ds1._k))) {
	      printf("bad _k, skipping");
	      validPoint = 1;
	    }

	    if (validPoint==0) {
	      imgCompCount++;
	      col -= col1;
	      theta[t] += col.Len2();

	      		
	      if (1) {  // NORMALIZE normal method
		
		if (t>0) {
		  calculateSplineK(ds, ds._k, j, invsubs, K);
		  rmmult(grad,G_t,K,3,2,8);      // G_t K_t
		}
		if (t<ul) {
		  calculateSplineK(ds1, ds1._k, j, invsubs, K1);  // K is 2x8
		  for (di=0; di<16; ++di)
		    K1[di] *=-1.;
		  rmmult(grad1, G_t1,K1,3,2,8);   // - G_(t+1) K_(t+1)
		}
		
#ifdef _TIME_
		tt7 += time7.elapsed();
		
		if (tt6 == 0) time6.start();
		else time6.restart();
		
#endif

		checkOk = true;
		if (t>0)
		  checkOk = checkOk & imgCache.checkVars(vars,0);
		if (t<ul)
		  checkOk = checkOk & imgCache.checkVars(vars1,1);
		
		if (checkOk) {
		  imgCache.newObs();
		  if (t>0)  imgCache.takeJac(grad,0); 
		  if (t<ul) imgCache.takeJac(grad1,1); 
		  imgCache.outer_prod(col.r());
		  imgCache.newObs();
		  if (t>0)  imgCache.takeJac(grad+8,0); 
		  if (t<ul) imgCache.takeJac(grad1+8,1); 
		  imgCache.outer_prod(col.g());
		  imgCache.newObs();
		  if (t>0)  imgCache.takeJac(grad+16,0); 
		  if (t<ul) imgCache.takeJac(grad1+16,1); 
		  imgCache.outer_prod(col.b());

		}

		else {
		  imgCache.recordMiss(3);

		  if (t==0) {
		    imgKeep->take8Key(vars1, grad1,    col.r(), col.g(), col.b()); 
		    //imgKeep->take8Key(vars1, grad1+8,  col.g());
		    //imgKeep->take8Key(vars1, grad1+16, col.b());
		  }
		  else if (t==ul) {
		    imgKeep->take8Key(vars, grad,    col.r(), col.g(), col.b());
		    //imgKeep->take8Key(vars, grad+8,  col.g());
		    //imgKeep->take8Key(vars, grad+16, col.b());		   
		  }
		  else {

		    imgKeep->take88(vars, vars1, grad,   grad1,     col.r());
		    imgKeep->take88(vars, vars1, grad+8, grad1+8,   col.g());
		    imgKeep->take88(vars, vars1, grad+16, grad1+16, col.b());		    
		  }
		} // end cache miss
	      } // end normalized-method
		
#ifdef _TIME_
	      tt6 += time6.elapsed();
#endif
		
	      /*
	      else {  // old/new non-normalized normal method
	    
		calculateNewK(K,0,j,invsubs); 
		if (t>0)
		  rmmult(grad,G_t,K,3,2,4);      // G_t K_t
		if (t<ul) {
		  for (di=0; di<8; di++)
		    K[di] *=-1.;
		  rmmult(grad1, G_t1,K,3,2,4);   // - G_(t+1) K_(t+1)
		}

		Ubcv J1 ( _mts->numVars() ), J2 ( _mts->numVars() ), J3( _mts->numVars() ); // for r,g,b observations
		
		//imgCache.newObs();
		if (t>0) { // SPEED: these J's share structure, obviously
		  putVarsOnJ(J1, ds01, vars01, grad[0], grad[1]);  // r
		  putVarsOnJ(J1, ds02, vars02, grad[2], grad[3]);
		  putVarsOnJ(J2, ds01, vars01, grad[4], grad[5]);  // g
		  putVarsOnJ(J2, ds02, vars02, grad[6], grad[7]);
		  putVarsOnJ(J3, ds01, vars01, grad[8], grad[9]);  // b
		  putVarsOnJ(J3, ds02, vars02, grad[10], grad[11]);
		  //imgCache.checkVars(vars01,0);
		  //imgCache.checkVars(vars02,0);
		}

		if (t<ul) {
		  putVarsOnJ(J1, ds11, vars11, grad1[0], grad1[1]);  // r
		  putVarsOnJ(J1, ds12, vars12, grad1[2], grad1[3]);
		  putVarsOnJ(J2, ds11, vars11, grad1[4], grad1[5]);  // g
		  putVarsOnJ(J2, ds12, vars12, grad1[6], grad1[7]);
		  putVarsOnJ(J3, ds11, vars11, grad1[8], grad1[9]);  // b
		  putVarsOnJ(J3, ds12, vars12, grad1[10], grad1[11]);
		  //imgCache.checkVars(vars11,1);
		  //imgCache.checkVars(vars12,1);
		}

		//printSparseVector(J1);
		//printSparseVector(J2);
		//printSparseVector(J3);
		
		//imgCache.obsFinished();
		imgKeep->takeJVector(J1, col.r());
		imgKeep->takeJVector(J2, col.g());
		imgKeep->takeJVector(J3, col.b());
	      } // end else
	      */

	      /*
	      Ubcv J1 ( _mts->numVars() ), J2 ( _mts->numVars() ), J3( _mts->numVars() ); // for r,g,b observations
	      if (t>0) 
		putGradOnJ(J1,J2,J3,grad,vars);
	      
	      if (t<ul) 
		putGradOnJ(J1,J2,J3,grad1,vars1);
	      
	      imgKeep->takeJVector(J1, col.r());
	      imgKeep->takeJVector(J2, col.g());
	      imgKeep->takeJVector(J3, col.b());
	      */
	    } // if validPoint

	    else {
	      imgCompCount++;
	      theta[t] += 5000;
	    }
	 
	  }// samples perp to curve (j)

	
	} // useImage

	// EDGE TERM 
	if (useEdges && n%stride==0 && n!=un && t>0 && t<ul &&_mts->useEdges(c) && _mts->_edgeMins[c][n] > 0) {
	  
	  splineLoc(ds,0,invsubs, &loc);
	  double eval = pyrmsE[t]->img->getFImage(level)->interpolate(loc.x(), loc.y(), &validPoint); 
	  G_t[0] = pyrmsE[t]->gradx->getFImage(level)->interpolate(loc.x(), loc.y(), &di);
	  G_t[1] = pyrmsE[t]->grady->getFImage(level)->interpolate(loc.x(), loc.y(), &di);
	  
	  if (validPoint == 1) {
	    noSubSplineLoc(ds,0,subs,&loc);
	    validPoint = !_mts->occluded(ds._loc, t);
	  }
	  
	  
	  // THis can happen if two adjacent samples fall on top of each other
	  if (validPoint==1 && (!finite(ds._k) || isnan(ds._k) ||
				!finite(ds1._k) || isnan(ds1._k))) {
	    printf("bad _k, skipping");
	    validPoint = 0;
	  }
	  
	  if (validPoint==1) {
	    edgeCompCount++;
	    const vector<double>& edgeMins = _mts->_edgeMins[c];
	    double tmp = eval / edgeMins[n];
	    thetaE[t] += tmp*tmp; //printf("%f\n",tmp*tmp);
	    
	    calculateSplineK(ds, ds._k, 0, invsubs, K);
	    rmmult(grad,G_t, K, 1,2,8);
	    // grad jac, eval residual, 1/(_edgeMins[n]*_edgeMins[n]) coefficient
	    
	    checkOk = ECache.checkVars(vars,0);
	    if (checkOk) {
	      ECache.newObs();
	      ECache.takeJac(grad,0);
	      ECache.outer_prod(eval,  1./(edgeMins[n]*edgeMins[n]) );
	    }
	    else {
	      ECache.recordMiss(1);
	      imgKeep->take8(vars,grad,eval, 1./(edgeMins[n]*edgeMins[n]));
	    }
	    
	  }
	  else {
	    edgeCompCount++;
	    thetaE[t] += 2;
	  }

	}
	

	
	// D0 TERM  ------------------------- 
	if (_useD0) {    

#ifdef _TIME_
	  if (tt10 == 0) time10.start();
	  else time10.restart();
#endif

	  ++smooth0Count;
	  Vec2f delta(ds._loc, ds1._loc);
	  thetas0[t] += delta.Len2();

	  /*if (t>0) {
	    memset(grad,0,9*sizeof(double));
	    grad[1] = ds._a;  grad[3] = ds._b;  grad[5] = ds._c;  grad[7] = ds._d;
	  }

	  if (t<ul) {
	    memset(grad1,0,9*sizeof(double));
	    grad1[1] = -ds1._a;  grad1[3] = -ds1._b;  grad1[5] = -ds1._c;  grad1[7] = -ds1._d;
	  }

	  if (tt6 == 0) time6.start();
	  else time6.restart();
	  if (t==0) {
	    D0Keep->take8Key(vars1, grad1+1, delta.x());
	    D0Keep->take8Key(vars1, grad1,   delta.y());
	  }
	  else if (t==ul) {
	    D0Keep->take8Key(vars, grad+1, delta.x());
	    D0Keep->take8Key(vars, grad,   delta.y());
	  }
	  else {
	    D0Keep->take88(vars, vars1, grad+1, grad1+1, delta.x());
	    D0Keep->take88(vars, vars1, grad, grad1, delta.y());
	    }*/

	  Ubcv J1( _mts->numVars() ), J2( _mts->numVars() );  // for x,y observations

#ifdef _TIME_
	  tt10 += time10.elapsed();
#endif
	  
	  checkOk = true;
	  if (t>0)
	    checkOk = checkOk & D0Cache.checkVars(vars01,0);
	  if (t<ul)
	    checkOk = checkOk & D0Cache.checkVars(vars11,1);

	  if (checkOk) {
	    D0Cache.newObs();
	    if (t>0)  D0Cache.takeJac(ds01,0, 1.,0);
	    if (t<ul) D0Cache.takeJac(ds11,1,-1.,0);
	    D0Cache.outer_prod(delta.x());
	    D0Cache.newObs();
	    if (t>0)  D0Cache.takeJac(ds01,0,0, 1.);
	    if (t<ul) D0Cache.takeJac(ds11,1,0,-1.);
	    D0Cache.outer_prod(delta.y());
	  }

	  else {
	    D0Cache.recordMiss(2);
	    if (t>0) {
	      putVarsOnJ(J1, J2, ds01, vars01, 1);
	    }
	    if (t<ul) {
	      putVarsOnJ(J1, J2, ds11, vars11, -1.);
	    }

	    D0Keep->takeJVector(J1, delta.x());
	    D0Keep->takeJVector(J2, delta.y());
	  }

#ifdef _TIME_
	  if (tt4 == 0) time4.start();
	  else time4.restart();
#endif

	  
	  /*
	  if (n==0 || n == un) {
	    printf("tcn: %d %d %d, loc: %.3f %.3f\n",  t,c,n,ds01._loc.x(), ds01._loc.y());
	    printSparseVector(J1);
	    printSparseVector(J2);
	    printf("sqrt(res): %.2f, %.2f\n\n",
		   delta.x(),delta.y());
		   }*/
	  //D0Cache.obsFinished();


#ifdef _TIME_
	  tt4 += time4.elapsed();
#endif
	  
	}
	// D1 TERM  -------------------------
	if (_useD1) {  
	  
	  
	  if (D1_parameter_based) {
	    /*
	    double tx, ty, r;	
	    r = ds._tangent.Len2() - ds1._tangent.Len2();
	    //assert(t==_mts->_numFrames-1 || fabs(r)<.01); 
	    thetas1[t] += r*r;

	    if (t>0) {
	      tx = ds._tangent.x();   ty = ds._tangent.y();
	      grad[0] = 2. * tx * ds._ap; grad[1] = 2. * ty * ds._ap;
	      grad[2] = 2. * tx * ds._bp; grad[3] = 2. * ty * ds._bp;
	      grad[4] = 2. * tx * ds._cp; grad[5] = 2. * ty * ds._cp;
	      grad[6] = 2. * tx * ds._dp; grad[7] = 2. * ty * ds._dp;
	    }

	    if (t<ul) {
	      tx = ds1._tangent.x();   ty = ds1._tangent.y();
	      grad1[0] = -2. * tx * ds1._ap; grad1[1] = -2. * ty * ds1._ap;
	      grad1[2] = -2. * tx * ds1._bp; grad1[3] = -2. * ty * ds1._bp;
	      grad1[4] = -2. * tx * ds1._cp; grad1[5] = -2. * ty * ds1._cp;
	      grad1[6] = -2. * tx * ds1._dp; grad1[7] = -2. * ty * ds1._dp;
	    }

	    if (t==0)
	      D1Keep->take8Key(vars1, grad1, r);
	    else if (t==ul)
	      D1Keep->take8Key(vars, grad, r);
	    else
	      D1Keep->take88(vars, vars1, grad, grad1, r);
	    ++smooth1Count;
	    */
	  }

	  else {  // sample-based D1
  
	    if (n==un || ds1num==un1) continue; // can't do last samples

	    //if (t==ul)
	    //printf("shit\n"); 
	    //TrackSample ds12;

	    /*
#ifdef _TIME_	    
	    if (tt5 == 0) time5.start();
	    else time5.restart();
#endif

            // get next and corresponding samples
	    if (t<ul)
	      _mts->getExistingCorrDiscreteSample(t,c,n+1,&ds12,vars12);
	    else
	      _mts->getExistingCorrDiscreteSample(t,c,n+1,&ds12);
#ifdef _TIME_
	    tt5 += time5.elapsed();
#endif
	    */

#ifdef _TIME_
	    if (tt10 == 0) time10.start();
	    else time10.restart();	  
#endif

	    Vec2f delta1, delta2;
	    Vec2f_Sub(delta1,ds01._loc, ds02._loc);
	    Vec2f_Sub(delta2, ds11._loc, ds12._loc);
	    double resid = delta1.Len2() - delta2.Len2();
	    thetas1[t] += resid*resid;

#ifdef _TIME_
	    tt10 += time10.elapsed();
#endif

	    Ubcv J ( _mts->numVars() );
	    checkOk = true;
	    if (t>0)
	      checkOk = checkOk & D1Cache.checkVars(vars01,0) & D1Cache.checkVars(vars02,0);
	    if (t<ul)
	      checkOk = checkOk & D1Cache.checkVars(vars11,1) & D1Cache.checkVars(vars12,1);
	    
	    if (checkOk) {
	      D1Cache.newObs();
	      if (t>0) {
		D1Cache.takeJac(ds01,0, 2*delta1.x(),  2*delta1.y());
		D1Cache.takeJac(ds02,0,-2*delta1.x(), -2*delta1.y());
	      }
	      if (t<ul) {
		D1Cache.takeJac(ds11,1,-2*delta2.x(),  -2*delta2.y());
		D1Cache.takeJac(ds12,1,2*delta2.x(), 2*delta2.y());
	      }
	      D1Cache.outer_prod(resid);
	      //D1Cache.obsFinished();
	    }

	    else {
	      D1Cache.recordMiss(1);
	      if (t>0) {
		putVarsOnJ(J, ds01, vars01, 2*delta1.x(),  2*delta1.y());
		putVarsOnJ(J, ds02, vars02, -2*delta1.x(), -2*delta1.y());
		
	      }
	      if (t<ul) {
		putVarsOnJ(J, ds11, vars11, -2*delta2.x(),  -2*delta2.y());
		putVarsOnJ(J, ds12, vars12, 2*delta2.x(), 2*delta2.y());
	      }
	      D1Keep->takeJVector(J, resid);
	    }

#ifdef _TIME_
	    if (tt4 == 0) time4.start();
	    else time4.restart();
#endif


#ifdef _TIME_
	    tt4 += time4.elapsed();
#endif
	    ++smooth1Count;
	  }
	}



	// D2 TERM  -------------------------
	if (_useD2) {  
	  
	  if (D2_parameter_based) {
	    /*
	    float speed2 = ds._tangent.x()*ds._tangent.x() + ds._tangent.y()*ds._tangent.y();
	    float speed12 = ds1._tangent.x()*ds1._tangent.x() + ds1._tangent.y()*ds1._tangent.y();
	    Vec2f curv1 = ds._curvature;  curv1 /= speed2;
	    Vec2f curv2 = ds1._curvature; curv2 /= speed12;
	    Vec2f delta(curv1, curv2);
	    //Vec2f delta(ds._curvature/speed2, ds1._curvature/speed12);
	    thetas[t] += delta.Len2();
	    //printf("c %d t %d n %d, %f\n",c,t,n,delta.Len2());
	  
	    if (t>0) {
	      memset(grad,0,9*sizeof(double));
	      grad[1] = ds._app/speed2;  grad[3] = ds._bpp/speed2;  
	      grad[5] = ds._cpp/speed2;  grad[7] = ds._dpp/speed2;
	    }
	  
	    if (t<ul) {
	      memset(grad1,0,9*sizeof(double));
	      grad1[1] = -ds1._app/speed12;  grad1[3] = -ds1._bpp/speed12; 
	      grad1[5] = -ds1._cpp/speed12;  grad1[7] = -ds1._dpp/speed12;
	    }
	  
	    if (t==0) {
	      keep->take8Key(vars1, grad1+1, delta.x());
	      keep->take8Key(vars1, grad1,   delta.y());
	    }
	    else if (t==ul) {
	      keep->take8Key(vars, grad+1, delta.x());
	      keep->take8Key(vars, grad,   delta.y());
	    }
	    else {
	      keep->take88(vars, vars1, grad+1, grad1+1, delta.x());
	      keep->take88(vars, vars1, grad, grad1, delta.y());
	    }
	    ++smooth2Count;
	    */
	  }
	  //--------------
	  
	  else {  // sample-based term for D2
	    
	    if (n==0 || n==un || ds1num==0 || ds1num ==un1) continue; // no sample terms for endpoints
	    //printf("time %d: %d,%.3f corr to %d,%.3f: %.2f %.2f %.2f %.2f, %.2f %.2f %.2f %.2f\n",t,n,ds._t+ds._baset,
	    //   ds1num, ds1._t + ds1._baset,
	    //   ds._a, ds._b, ds._c, ds._d, ds1._a, ds1._b, ds1._c, ds1._d);
	    
	    
	    // get relevant samples for smoothness terms
	    // SPEED: can keep these between iterations, reducing calls by 2/3
#ifdef _TIME_
	    if (tt5 == 0) time5.start();
	    else time5.restart();
#endif

	    TrackSample ds00,ds10;
	    if (t>0) 
	      _mts->getDiscreteSample(t,  c, n-1,      &ds00, vars00);       
	    else 
	      _mts->getDiscreteSample(t,  c, n-1,      &ds00);
	    	    
	    if (t<ul) 
	      //_mts->getDiscreteSample(t+1,c, ds1num-1, &ds10, vars10);
	      _mts->getExistingCorrDiscreteSample(t,c, n-1, &ds10, vars10);
	    else 
	      //_mts->getDiscreteSample(t+1,c, ds1num-1, &ds10);
	    _mts->getExistingCorrDiscreteSample(t,c, n-1, &ds10);

#ifdef _TIME_
	    tt5 += time5.elapsed();
#endif

	    // calculate residual
#ifdef _TIME_
	    if (tt10 == 0) time10.start();
	    else time10.restart();	  
#endif
	    Vec2f delta = ds01._loc; delta *=-2;
	    delta += ds00._loc; delta += ds02._loc;
	    Vec2f_AddScale(delta, delta, ds11._loc, 2.);
	    delta -= ds10._loc; delta -= ds12._loc;
	    thetas[t] += delta.Len2();
	    
	    Ubcv J1 ( _mts->numVars() ), J2 ( _mts->numVars() ); // for x,y observations
#ifdef _TIME_
	    tt10 += time10.elapsed();
#endif

	    checkOk = true;
	    if (t>0)
	      checkOk = checkOk & D2Cache.checkVars(vars00,0) & D2Cache.checkVars(vars01,0) & D2Cache.checkVars(vars02,0);
	    if (t<ul)
	      checkOk = checkOk & D2Cache.checkVars(vars10,1) & D2Cache.checkVars(vars11,1) & D2Cache.checkVars(vars12,1);
	    
	    if (checkOk) {
	      D2Cache.newObs();
	      if (t>0)  {
		D2Cache.takeJac(ds00,0, 1.,0);
		D2Cache.takeJac(ds01,0, -2.,0);
		D2Cache.takeJac(ds02,0, 1.,0);
	      }
	      if (t<ul) {
		D2Cache.takeJac(ds10,1, -1.,0);
		D2Cache.takeJac(ds11,1, 2.,0);
		D2Cache.takeJac(ds12,1, -1.,0);
	      }
	      D2Cache.outer_prod(delta.x());
	      D2Cache.newObs();
	      if (t>0)  {
		D2Cache.takeJac(ds00,0,0, 1.);
		D2Cache.takeJac(ds01,0,0, -2.);
		D2Cache.takeJac(ds02,0,0, 1.);
	      }
	      if (t<ul) {
		D2Cache.takeJac(ds10,1,0, -1.);
		D2Cache.takeJac(ds11,1,0, 2.);
		D2Cache.takeJac(ds12,1,0, -1.);
	      }
	      D2Cache.outer_prod(delta.y());
	    }

	    else {
	      D2Cache.recordMiss(2);
	      if (t>0) {
		putVarsOnJ(J1, J2, ds00, vars00, 1.);
		putVarsOnJ(J1, J2, ds01, vars01, -2.);
		putVarsOnJ(J1, J2, ds02, vars02, 1.);
	      }
	      if (t<ul) {
		putVarsOnJ(J1, J2, ds10, vars10, -1.);
		putVarsOnJ(J1, J2, ds11, vars11, 2.);
		putVarsOnJ(J1, J2, ds12, vars12, -1.);
	      }
	      keep->takeJVector(J1, delta.x());
	      keep->takeJVector(J2, delta.y());
	    }

	    //if (ds._t > .99 && t>0)
	    //printf("shit\n");
#ifdef _TIME_
	    if (tt4 == 0) time4.start();
	    else time4.restart();
#endif

	    //D2Cache.obsFinished();


#ifdef _TIME_
	    tt4 += time4.elapsed();
#endif

	    ++smooth2Count;

	  } // end sample-based D2

	} // if useD2
	
      } //samples on curve (n)

    } // time (t)
    
  } // curves

  imgCache.writeBack(); D0Cache.writeBack();
  D1Cache.writeBack(); D2Cache.writeBack();  ECache.writeBack();
  printf("Image hit rate %.5f\n",imgCache.hitRate());
  printf("D0 hit rate %.5f\n",D0Cache.hitRate());
  printf("D1 hit rate %.5f\n",D1Cache.hitRate());
  printf("D2 hit rate %.5f\n",D2Cache.hitRate());
  printf("E hit rate %.5f\n",ECache.hitRate());
#ifdef _TIME_
  if (tt2 == 0) time2.start();
  else time2.restart();
#endif  

  if (smooth2Count>0) {
    invd2 = _mts->_numFrames * smooth2Deriv/ double(smooth2Count);  
    printf("invd2 %f\n",invd2);
    keep->scalarMult(invd2);   
  }
  else invd2 = 0;
  
  if (smooth0Count>0) {
    invd0 = _mts->_numFrames * smooth0Deriv/ double(smooth0Count);   
    D0Keep->scalarMult(invd0);
    keep->addMat(D0Keep);
  }
  else invd0 = 0;

  if (smooth1Count>0) {
    invd1 = _mts->_numFrames * smooth1Deriv/ double(smooth1Count);   
    D1Keep->scalarMult(invd1);
    keep->addMat(D1Keep);
  }
  else invd1 = 0;
  

  if (imgCompCount>0) {
    invcw = _mts->_numFrames * 1. / double(imgCompCount);   
    imgKeep->scalarMult(invcw);
    keep->addMat(imgKeep);
  }
  else {
    printf("No image count\n");
    invcw = 0;
  }
  if (edgeCompCount>0) {
    invew = edgeWeight / double(edgeCompCount);
    edgeKeep->scalarMult(invew);
    keep->addMat(edgeKeep);
  }
  else {
    printf("No edge count\n");
    invew = 0;
  }

#ifdef _TIME_
  tt2 += time2.elapsed();
#endif

  //keep->outputMat("B.dat");
  //keep->shit();    
  //std::exit(0);

#ifdef _TIME_
  if (tt11 == 0) time11.start();
  else time11.restart();
#endif

  double thetaSum = 0;
  int i;
  for (i=0; i<numFrames; ++i) 
    thetaSum += theta[i]*invcw + thetas1[i]*invd1 + 
      thetas[i]*invd2 + thetas0[i]*invd0 + thetaE[i]*invew;
  
  if (_stateOk) {
    fprintf(stdout,"%9f = ",thetaSum);
    
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",theta[i]*invcw);   // image
    fprintf(stdout,"\nD2          ");
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",thetas[i]*invd2);  // D2 
    fprintf(stdout,"\nD1          ");
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",thetas1[i]*invd1); // D1
    fprintf(stdout,"\nD0          ");
    for (i=0; i<numFrames; ++i) 
      fprintf(stdout,"%.3f + ",thetas0[i]*invd0); // D0
    fprintf(stdout,"\nEdge        ");
    for (i=0; i<numFrames-1; ++i) 
      fprintf(stdout,"%.3f + ",thetaE[i]*invew);  // Edge
    printf("\n");
  }  


  //std::exit(0);
#ifdef _TIME_
  tt11 += time11.elapsed();
#endif

  tt1 += time1.elapsed();

  DELETE_CSM;
  return thetaSum;
}

void KLT_TrackingContext::putGradOnJ(Ubcv& J1, Ubcv& J2, Ubcv& J3, const double grad[24], const int vars[4]) {
  if (tt6 == 0) time6.start();
  else time6.restart();
  J1[vars[0]] += grad[0];   J1[vars[0]+1] += grad[1];
  J1[vars[1]] += grad[2];   J1[vars[1]+1] += grad[3];
  J1[vars[2]] += grad[4];   J1[vars[2]+1] += grad[5];
  J1[vars[3]] += grad[6];   J1[vars[3]+1] += grad[7];

  J2[vars[0]] += grad[8];   J2[vars[0]+1] += grad[9];
  J2[vars[1]] += grad[10];   J2[vars[1]+1] += grad[11];
  J2[vars[2]] += grad[12];   J2[vars[2]+1] += grad[13];
  J2[vars[3]] += grad[14];   J2[vars[3]+1] += grad[15];
  
  J3[vars[0]] += grad[16];   J3[vars[0]+1] += grad[17];
  J3[vars[1]] += grad[18];   J3[vars[1]+1] += grad[19];
  J3[vars[2]] += grad[20];   J3[vars[2]+1] += grad[21];
  J3[vars[3]] += grad[22];   J3[vars[3]+1] += grad[23];
  tt6 += time6.elapsed();
}

void KLT_TrackingContext::putVarsOnJ(Ubcv& J1, Ubcv& J2, const TrackSample& ds, const int vars[4], const double coef) {
  if (tt3 == 0) time3.start();
  else time3.restart();
  J1[vars[0]]   += coef * ds._a;
  J2[vars[0]+1] += coef * ds._a;
  J1[vars[1]]   += coef * ds._b;
  J2[vars[1]+1] += coef * ds._b;
  J1[vars[2]]   += coef * ds._c;
  J2[vars[2]+1] += coef * ds._c;
  J1[vars[3]]   += coef * ds._d;
  J2[vars[3]+1] += coef * ds._d;
  tt3 += time3.elapsed();
}

void KLT_TrackingContext::putVarsOnJ(Ubcv& J, const TrackSample& ds, const int vars[4], const double c1, const double c2) {
  if (tt3 == 0) time3.start();
  else time3.restart();
  J[vars[0]]   += c1 * ds._a;
  J[vars[0]+1] += c2 * ds._a;
  J[vars[1]]   += c1 * ds._b;
  J[vars[1]+1] += c2 * ds._b;
  J[vars[2]]   += c1 * ds._c;
  J[vars[2]+1] += c2 * ds._c;
  J[vars[3]]   += c1 * ds._d;
  J[vars[3]+1] += c2 * ds._d;
  tt3 += time3.elapsed();
}


// Given a TrackSample with loc & normal, calculate location j off
void KLT_TrackingContext::noSubSplineLoc(const TrackSample& ds, int j, float subs, Vec2f* loc) { 
  *loc = ds._loc;
  //*loc *= invsubs;  // downsample before adding tangent term
  if (j!=0) 
    Vec2f_AddScale(*loc, *loc, ds._nnormal, float(j) * subs);  
}


// Given a TrackSample with loc & normal, calculate location j off
void KLT_TrackingContext::splineLoc(const TrackSample& ds, int j, float invsubs, Vec2f* loc) { 
  *loc = ds._loc;
  *loc *= invsubs;  // downsample before adding tangent term
  if (j!=0) 
    Vec2f_AddScale(*loc, *loc, ds._nnormal, float(j));  
}

// Calculate 2x8 derivative K matrix for sample, j off tangent
void KLT_TrackingContext::calculateSplineK(const TrackSample& s, double k,
					   int j, double invsubs, double K[16]) {  
  assert(finite(k)  && !isnan(k));
  //double t = s._t;
  //double mt = 1.-t;
  memset(K,0,16*sizeof(double));
  assert(invsubs > 0);

  // from sample itself
  K[0] = invsubs * s._a;    K[9] =  K[0];
  K[2] = invsubs * s._b;    K[11] = K[2];
  K[4] = invsubs * s._c;    K[13] = K[4];
  K[6] = invsubs * s._d;    K[15] = K[6];
  

  // from march off sample
  if (j==0) return;
  double jk = double(j) * k;
  assert(k > 0);
  K[8]  = jk * s._ap;   K[1] = -K[8];  //a'
  K[10] = jk * s._bp;   K[3] = -K[10]; //b'
  K[12] = jk * s._cp;   K[5] = -K[12]; //c'
  K[14] = jk * s._dp;   K[7] = -K[14]; //d'

}



void KLT_TrackingContext::initSplineEdgeMins(const KLT_FullPyramid** pyrmsE, const int level) {
  printf("initing Edge mins\n");
  int c, l=0, hw, n, ut = _mts->_numFrames-1, validPoint, di;
  double eval, eval2;
  Vec2f loc;
  float finvsubs = 1./float(pow2(level));
  TrackSample ds0,dst;
  //float finvsubs = (float) invsubs;

  for (c=0; c < _mts->_nCurves; c++) { // iterate over curves

    //CurveRecord* cr = _mt->_crs + c;
    //if (!cr->_useEdges) continue;
    if (!_mts->useEdges(c)) continue;
    //l = cr->_nPoints;
    l = _mts->getNumSamples(0,c);  assert(l>0);
    vector<double>& edgeMins = *(_mts->refreshEdgeMins(c,l-1));

    //if (cr->_edgeMin) 
    // delete[] cr->_edgeMin;
    //cr->_edgeMin = new double[2*l-1];

    for (n=0; n<l-1; n++) { // iterate over curve points

      if (n==l-2)
	hw=2;
      else
	hw=1;

      //for (i = 0 ; i <= hw ; ++i, ++ict)  {  

      validPoint = 1;
      //assert(cr->_edgeMin);

      _mts->getDiscreteSample(0,c,n,&ds0);
      _mts->getExistingCorrDiscreteSample(ut-1,c,n,&dst);  // second to last -> last
      splineLoc(ds0,0,finvsubs,&loc);
      //hnew2(n,c,0,i,0,finvsubs, _mt->_Z, loc);
      eval = pyrmsE[0]->img->getFImage(level)->interpolate(loc.x(), loc.y(), &di);
      validPoint = MIN( di, validPoint);
      //hnew2(n,c,ut,i,0,finvsubs,_mt->_Z, loc);
      splineLoc(dst,0,finvsubs, &loc);
      eval2 = pyrmsE[ut]->img->getFImage(level)->interpolate(loc.x(), loc.y(), &di);
      //eval2 = 0;  // ONE
      validPoint = MIN( di, validPoint); // ONE
      if (validPoint==1)
	edgeMins.push_back(MAX(eval, eval2));
      else
	edgeMins.push_back(-1);
      //assert(cr->_edgeMin[ict] ==-1);
	

    } // end iterate over curve points
    assert((int) edgeMins.size() == l-1);
  } // end iterate over curves
  printf("finished initing Edge mins with %d samples\n",l-1);  
}

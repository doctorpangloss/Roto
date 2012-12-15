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



#include "multiKeeper.h"

MultiKeeper::MultiKeeper(MultiKeeper* other) {
  _mt = other->_mt;
  _precond = false;
  int c;

  _b_oelems = new OffElems*[_mt->_nCurves];  
  _e_oelems = new OffElems*[_mt->_nCurves];
  memset(_b_oelems, 0, _mt->_nCurves*sizeof(OffElems*));
  memset(_e_oelems, 0, _mt->_nCurves*sizeof(OffElems*));


  _g = NULL;
  _gDirty = true;

  _keepers = new Keeper[_mt->_nCurves];
  for (c=0; c<_mt->_nCurves; c++)
    _keepers[c].init(&(other->_keepers[c]));
  _subsNvar = other->_subsNvar;
  //_sub_offset = other->_sub_offset;
  _sub_offset = new int[_mt->_nCurves];
  memcpy(_sub_offset, other->_sub_offset, _mt->_nCurves*sizeof(int));
  
  _sweeper = NULL;
}

MultiKeeper::MultiKeeper(const MultiTrackData* mt, int curveSampling, int hold, bool preconde) : _mt(mt), _precond(preconde) {

  //assert(curveSampling==1); 
  _b_oelems = new OffElems*[mt->_nCurves];  
  _e_oelems = new OffElems*[mt->_nCurves];
  memset(_b_oelems, 0, mt->_nCurves*sizeof(OffElems*));
  memset(_e_oelems, 0, mt->_nCurves*sizeof(OffElems*));
  _g = NULL;
  _gDirty = true;

  _keepers = new Keeper[mt->_nCurves];
  _subsNvar=0;
  _sub_offset = new int[mt->_nCurves];
  for (int c=0; c<mt->_nCurves; c++) {
    _keepers[c].init(mt->_numFrames-1, mt->_crs[c]._nVars*2, curveSampling, hold);

    if (c==0)
      _sub_offset[c] = 0;
    else
      _sub_offset[c] = _sub_offset[c-1] + _keepers[c-1].numVar();

    _subsNvar += _keepers[c].numVar() / 2;
  }

  if (precond()) {
    int index, i, j;
    vector<int>* fixers = new vector<int>[mt->_nCurves];
    int* numvars = new int[mt->_nCurves];
    for (i=0; i<_mt->_nCurves; ++i) {
      const FixedV& flocs = _mt->_crs[i]._fixedLocs;
      for (j=0; j<flocs.size(); ++j) { //SPEED, not good
	index = internal_cp_to_subvar(i, flocs[j]._i); // will only work for endpoints (because of curve subdv)
	index += flocs[j]._frame * internal_cp_to_keeper(i,flocs[j]._i)->numVarOneBlock();
	assert (index+1 < _subsNvar*2);
	fixers[i].push_back(index);
      }      
      numvars[i] = _keepers[i].numVarOneBlock()/2;
    }
    _sweeper = new HB_Sweep(mt->_nCurves, numvars, mt->_numFrames-1, fixers);
    delete[] fixers;
    delete[] numvars;
  }
  else
    _sweeper=NULL;

}

MultiKeeper::~MultiKeeper(){
  delete[] _keepers;
  delete[] _sub_offset;

  for (int i=0; i<_mt->_nCurves; i++) {
    if (_b_oelems[i])
      delete _b_oelems[i];
    if (_e_oelems[i])
      delete _e_oelems[i];
  }
  delete[] _b_oelems;
  delete[] _e_oelems;

  if (_g) delete[] _g;
  if (_sweeper) delete _sweeper;
}

void MultiKeeper::add(MultiKeeper* other, double c) {
  int i;
  for (i=0; i<_mt->_nCurves; ++i) {
    _keepers[i].add(&(other->_keepers[i]), c);
    if (_b_oelems[i])
      _b_oelems[i]->add(other->_b_oelems[i], c);
    if (_e_oelems[i])
      _e_oelems[i]->add(other->_e_oelems[i], c);
  }
    
}


void OffElems::add(OffElems* other, double c) {
  assert(other);  
  assert(other->_nFrames == _nFrames);
  int ml = _nFrames*8, ol = (_nFrames-1)*8;
  int i;
  double *optr1, *ptr1;
  double *optr2, *ptr2;

  optr1 = other->_main; ptr1 = _main;
  for (i=0; i<ml; ++i, ++optr1, ++ptr1)
    *ptr1 += c * *optr1;
  
  optr1 = other->_left;  ptr1 = _left;
  optr2 = other->_right; ptr2 =_right;
  for (i=0; i<ol; ++i, ++optr1, ++ptr1,++optr2, ++ptr2) {
    *ptr1 += c * *optr1;
    *ptr2 += c * *optr2;
  }
  
  
}



// r = B*x
void MultiKeeper::matVecMult(const double x[], double r[]) const {
  memset(r,0,_subsNvar*2*sizeof(double));
  int i;

  const double *mx=x;
  double *mr=r;
  for (i=0; i<_mt->_nCurves; ++i) {
    _keepers[i].matVecMult(mx, mr);
    mx += _keepers[i].numVar();
    mr += _keepers[i].numVar();
    
    if (_b_oelems[i]) 
      _b_oelems[i]->matVecMult(x,r);

    if (_e_oelems[i]) 
      _e_oelems[i]->matVecMult(x,r);
  }

  /*
  FILE* fp;  
  fp = fopen("/homes/gws/aseem/r.txt","w");
  for (i=0; i<_subsNvar*2; ++i)
    fprintf(fp,"%.10f\n",r[i]);
  fclose(fp);

  fp = fopen("/homes/gws/aseem/x.txt","w");
  for (i=0; i<_subsNvar*2; ++i)
    fprintf(fp,"%.10f\n",x[i]);
  fclose(fp);
  exit(0);
  */
}

void OffElems::matVecMult(const double x[], double r[]) const {
  int nf1 = _nFrames-1;
  int t, j;
  const double *ms;
  int mc, mr;

  //FILE *fp = fopen("/homes/gws/aseem/shit.dat", "w");

  // main strip
  mc = _k; mr = _l; ms = _main;
  for (t=0; t<_nFrames; ++t, mr += _Oh, mc += _Ow) { 
    for (j=0; j<4; ++j) {
      *(r+mr+j) += *ms * *(x+mc);      //fprintf(fp,"%d %d 1\n",mr+j+1, mc+1);
      *(r+mc) += *ms * *(x+mr+j);       //fprintf(fp,"%d %d 1\n", mc+1,mr+j+1);
      ++ms; 
      *(r+mr+j) += *ms * *(x+mc+1);    //fprintf(fp,"%d %d 1\n",mr+j+1, mc+1+1);
      *(r+mc+1) += *ms * *(x+mr+j);       //fprintf(fp,"%d %d 1\n", mc+1+1,mr+j+1);
      ++ms; 
    }
  }

  // right strip
  mc = _k+_Ow; mr = _l; ms = _right;
  for (t=0; t<nf1; ++t, mr += _Oh, mc += _Ow) { 
    for (j=0; j<4; ++j) {
      *(r+mr+j) += *ms * *(x+mc);      //fprintf(fp,"%d %d 1\n",mr+j+1, mc+1);
      *(r+mc) += *ms * *(x+mr+j);       //fprintf(fp,"%d %d 1\n", mc+1,mr+j+1);
      ++ms; 
      *(r+mr+j) += *ms * *(x+mc+1);    //fprintf(fp,"%d %d 1\n",mr+j+1, mc+1+1);
      *(r+mc+1) += *ms * *(x+mr+j);       //fprintf(fp,"%d %d 1\n", mc+1+1,mr+j+1);
      ++ms; 
    }
  }

  // left strip
  mc = _k; mr = _l+_Oh; ms = _left;
  for (t=0; t<nf1; ++t, mr += _Oh, mc += _Ow) { 
    for (j=0; j<4; ++j) {
      *(r+mr+j) += *ms * *(x+mc);        //fprintf(fp,"%d %d 1\n",mr+j+1, mc+1);
      *(r+mc) += *ms * *(x+mr+j);       //fprintf(fp,"%d %d 1\n", mc+1,mr+j+1);
      ++ms; 
      *(r+mr+j) += *ms * *(x+mc+1);      //fprintf(fp,"%d %d 1\n",mr+j+1, mc+1+1);
      *(r+mc+1) += *ms * *(x+mr+j);       //fprintf(fp,"%d %d 1\n", mc+1+1,mr+j+1);
      ++ms; 
    }
  }
  //fclose(fp); 
  //exit(0);
}  
   

void MultiKeeper::printSparse(FILE* fp) {
  int i, mx=0;
  for (i=0; i<_mt->_nCurves; ++i) {
    _keepers[i].printSparse(mx, fp);    
    mx += _keepers[i].numVar();

    if (_b_oelems[i])
      _b_oelems[i]->printSparse(fp);
    if (_e_oelems[i])
      _e_oelems[i]->printSparse(fp);
  }
}

void OffElems::printSparse(FILE* fp) const {
  int nf1 = _nFrames-1;
  int t, j;
  const double *ms;
  int mc, mr;
  
   // main strip
  mc = _k; mr = _l; ms = _main;
  for (t=0; t<_nFrames; ++t, mr += _Oh, mc += _Ow) { 
    for (j=0; j<4; ++j) {
      fprintf(fp,"%d %d %.10f\n",mr+j+1, mc+1, *ms);
      fprintf(fp,"%d %d %.10f\n", mc+1,mr+j+1, *ms);
      ++ms; 
      fprintf(fp,"%d %d %.10f\n",mr+j+1, mc+1+1, *ms);
      fprintf(fp,"%d %d %.10f\n", mc+1+1,mr+j+1, *ms);
      ++ms; 
    }
  }

  // right strip
  mc = _k+_Ow; mr = _l; ms = _right;
  for (t=0; t<nf1; ++t, mr += _Oh, mc += _Ow) { 
    for (j=0; j<4; ++j) {
      fprintf(fp,"%d %d %.10f\n",mr+j+1, mc+1, *ms);
      fprintf(fp,"%d %d %.10f\n", mc+1,mr+j+1, *ms);
      ++ms; 
      fprintf(fp,"%d %d %.10f\n",mr+j+1, mc+1+1, *ms);
      fprintf(fp,"%d %d %.10f\n", mc+1+1,mr+j+1, *ms);
      ++ms; 
    }
  }

  // left strip
  mc = _k; mr = _l+_Oh; ms = _left;
  for (t=0; t<nf1; ++t, mr += _Oh, mc += _Ow) { 
    for (j=0; j<4; ++j) {
      fprintf(fp,"%d %d %.10f\n",mr+j+1, mc+1, *ms);
      fprintf(fp,"%d %d %.10f\n", mc+1,mr+j+1, *ms);
      ++ms; 
      fprintf(fp,"%d %d %.10f\n",mr+j+1, mc+1+1, *ms);
      fprintf(fp,"%d %d %.10f\n", mc+1+1,mr+j+1, *ms);
      ++ms; 
    }
  }
}

void MultiKeeper::refresh(){  // set B, g to zero
  
  int c;
  for (c=0; c<_mt->_nCurves; c++) {
    if (_b_oelems[c])
      _b_oelems[c]->refresh();
    if (_e_oelems[c])
      _e_oelems[c]->refresh();
  }
  
  _gDirty = true;
  
  for (c=0; c<_mt->_nCurves; c++) 
    _keepers[c].refresh();
} 
  

double MultiKeeper::calculateRo(const double* p, const double newTheta, double const oldTheta) const {
  int nv = _subsNvar*2;
  double num = oldTheta - newTheta;
  double denom = vecDot(nv, g(), p);
  double *Bp = new double[nv];
  matVecMult(p, Bp);
  denom += .5 * vecDot(nv, p, Bp);
  denom *= -1.;
  if (denom == 0) {
    printf("Calculate Ro: num %.5f, dem 0!\n",num);
    delete[] Bp;
    return 1;
  }
  printf("Calculate Ro: num %.5f, dem %.5f\n", num, denom);
  if (denom < 0) {
    printf("Bug: denom < 0\n");
    delete[] Bp; return 0;
  }
  double result = num/denom;
  delete[] Bp;
  return result; 
}
    
void MultiKeeper::handleConstraints(double* r) { 
  int index, i, j;
  for (i=0; i<_mt->_nCurves; ++i) {

    if (!precond()) {
      const FixedV& flocs = _mt->_crs[i]._fixedLocs;
      for (j=0; j<flocs.size(); ++j) { //SPEED, not good
	index = internal_cp_to_subvar(i, flocs[j]._i); // will only work for endpoints (because of curve subdv)
	index += flocs[j]._frame * internal_cp_to_keeper(i,flocs[j]._i)->numVarOneBlock();
	assert (index+1 < _subsNvar*2);
	r[index] = r[index+1] = 0;
	//printf("Doing set\n");
      } 
    }

    else
      _sweeper->zeroFixed(r);

  }
}

double MultiKeeper::createTestSol(const Vec2f* Z, const double* sol, Vec2f* Z2) const{ 
  int zup=0;
  double maxStep = -1, tmp;
  const double* msol = sol;

  for (int c=0; c<_mt->_nCurves; ++c) {
    zup += _mt->_crs[c]._nVars;                   // skip key 0 data

    tmp = _keepers[c].createTestSol(Z+zup, msol, Z2+zup);

    if (tmp > maxStep) {
      maxStep=tmp;
      assert(maxStep < 1000);
    }

    msol += _keepers[c].numVar();               // Up solution pointer by what keeper thinks it should (subsampling!)
    zup += _mt->_numFrames * _mt->_crs[c]._nVars; // add rest data, plus skip key 1
  }

  return maxStep;
}


const double* MultiKeeper::g() const { 
  if (_g && !_gDirty) return _g;
  
  if (!_g)
    _g = new double[_subsNvar*2];
  
  if (_gDirty) {
    double* gptr=_g;
    for (int i=0; i<_mt->_nCurves; i++) {
      const double* o = _keepers[i].g();
      memcpy(gptr, o, _keepers[i].numVar()*sizeof(double));
      gptr +=  _keepers[i].numVar();
    }
    _gDirty = false;
  }

  return _g;
}
  
//-----------------------------------------------------------------------------------------------------------
/*
void MultiKeeper::debug_jac(int t, int c, int cp, int size, const double* g, double res, double cons) {
  DebugJacobian dj;
  for (int i=0; i<size/2; ++i) {
    int var = internal_cp_to_var(c, cp+i);
    int curve = _mt->_crs[c].varToCurve(cp+i);
    dj.add(var + t*_mt->_crs[curve]._nVars*2, g[2*i]);
    dj.add(var + t*_mt->_crs[curve]._nVars*2+1, g[2*i+1]);
    dj.cons = cons;
    dj.res = res;
  }
  djacs.push_back(dj);
}

void MultiKeeper::debug_jac(int t, int c, int cp, int size, const double* g1, const double* g2, double res, double cons) {
  DebugJacobian dj;
  for (int i=0; i<size/2; i++) {
    int var = internal_cp_to_var(c, cp+i);
    int curve = _mt->_crs[c].varToCurve(cp+i);
    dj.add(var + t*_mt->_crs[curve]._nVars*2, g1[2*i]);
    dj.add(var + t*_mt->_crs[curve]._nVars*2+1, g1[2*i+1]);
    dj.add(var + (t+1)*_mt->_crs[curve]._nVars*2, g2[2*i]);
    dj.add(var + (t+1)*_mt->_crs[curve]._nVars*2+1, g2[2*i+1]);
    dj.cons = cons;
    dj.res = res;
  }
  djacs.push_back(dj);
}
*/
void MultiKeeper::outputg() {
  FILE* fp = fopen("/homes/gws/aseem/g.txt","w");
  const double* ge = g();
  for (int i=0; i<_subsNvar*2; i++)
    fprintf(fp,"%.5f\n",ge[i]);
  fclose(fp);
}

void MultiKeeper::outputSol(double* x) {
  FILE* fp = fopen("/homes/gws/aseem/sol.txt","w");  
  for (int i=0; i<_subsNvar*2; i++)
    fprintf(fp,"%.5f\n",x[i]);
  fclose(fp);
}

//-----------------------------------------------------------------------------------------------------------

void MultiKeeper::take2Key0(int curve, int ib, const double* g2, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  //int nib, maxj;
  //maxj = _keepers[curve].convert1(2*ib, 1, nib, g2, _G1);
  //_keepers[curve]._takejKey0(maxj, nib, _G1, res, c);
  internal_cp_to_keeper(curve, ib)->take2Key0(internal_cp_to_keeper_var(curve, ib),g2,res,c);
  //debug_jac(0,curve,ib,2,g2,res,c);
}

void MultiKeeper::take2(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c) { 
  memset(_G1,0,6*sizeof(double));
  memset(_G2,0,6*sizeof(double));
  //int nib, maxj;
  //maxj = _keepers[curve].convert2(2*ib, 1, nib, g1, g2, _G1, _G2);
  //_keepers[curve]._takej(maxj, b, nib, _G1, _G2, res, c);
  internal_cp_to_keeper(curve, ib)->take2(b, internal_cp_to_keeper_var(curve, ib), g1,g2, res, c);
  //debug_jac(b,curve,ib,2,g1,g2,res,c);
}

void MultiKeeper::take2Key1(int curve, int ib, const double* g1, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  //int nib, maxj;
  //maxj = _keepers[curve].convert1(2*ib, 1, nib, g1, _G1);
  //_keepers[curve]._takejKey1(maxj, nib, _G1, res, c);
  internal_cp_to_keeper(curve, ib)->take2Key1(internal_cp_to_keeper_var(curve, ib),g1,res,c);
  //debug_jac(_mt->_numFrames-2,curve,ib,2,g1,res,c);
}
  
void MultiKeeper::take4Key0(int curve, int ib, const double* g2, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  if (ib==0 && !_mt->begWithRest(curve)) { 
     int nib1, nib2, maxj1, maxj2;
     convert1( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g2); assert(maxj1==2);
     _keepers[curve]._takejKey0(maxj2, nib2, _G1+2, res, c);
     internal_cp_to_keeper(curve, 0)->_takejKey0(maxj1,nib1, _G1, res,c);
     processOElemsKey0(true, curve,1,0, maxj1+maxj2, _G1, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-2 && !_mt->endWithRest(curve)) {
     int nib1, nib2, maxj1, maxj2;
     convert1( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g2); assert(maxj2==2);
     _keepers[curve]._takejKey0(maxj1, nib1, _G1, res, c);
     internal_cp_to_keeper(curve, ib+1)->_takejKey0(maxj2,nib2, _G1+maxj1, res,c);
     processOElemsKey0(false, curve, ib-1, ib+1, maxj1+maxj2, _G1, c);
  }
  else {
    int nib, maxj;
    int lib = internal_cp_to_keeper_var(curve, ib);
    maxj = _keepers[curve].convert1(lib, 2, nib, g2, _G1);
    _keepers[curve]._takejKey0(maxj, nib, _G1, res, c);
  }
  //debug_jac(0,curve,ib,4,g2,res,c);
}

  /*
void MultiKeeper::take4Key0(int curve, int ib, const double* g2, double res, const double c) {
  assert(ib < _mt->_crs[curve]._nPoints-1);
  
  if (ib==0 && !_mt->begWithRest(curve)) { 
    _keepers[curve].take2Key0(internal_cp_to_keeper_var(curve, 1), g2+2, res, c);
    internal_cp_to_keeper(curve, 0)->take2Key0(internal_cp_to_keeper_var(curve, 0), g2, res, c);
    processOElemsKey0(true, curve, 1, 0, 4, g2, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-2 && !_mt->endWithRest(curve)) {
    _keepers[curve].take2Key0(internal_cp_to_keeper_var(curve, ib), g2, res, c);
    internal_cp_to_keeper(curve, ib+1)->take2Key0(internal_cp_to_keeper_var(curve, ib+1), g2+2, res, c);
    processOElemsKey0(false, curve, ib-1, ib+1, 4, g2, c);
  }
  else {
    _keepers[curve].take4Key0(internal_cp_to_keeper_var(curve, ib), g2,res,c);
    assert(_keepers+curve == internal_cp_to_keeper(curve, ib));
  }
  
  //debug_jac(0,curve,ib,4,g2,res,c);
}
  */

void MultiKeeper::take4(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  memset(_G2,0,6*sizeof(double));
  if (ib==0 && !_mt->begWithRest(curve)) { 
     int nib1, nib2, maxj1, maxj2;
    convert2( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g1,g2);  assert(maxj1==2);
    _keepers[curve]._takej(maxj2, b, nib2, _G1+2, _G2+2, res, c);
    internal_cp_to_keeper(curve, 0)->_takej(maxj1,b, nib1, _G1, _G2, res,c);
    processOElems(true, b, curve, maxj1+maxj2, _G1, _G2, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-2 && !_mt->endWithRest(curve)) {
     int nib1, nib2, maxj1, maxj2;
    convert2( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g1,g2);  assert(maxj2==2);
    _keepers[curve]._takej(maxj1, b, nib1, _G1, _G2, res, c);
    internal_cp_to_keeper(curve, ib+1)->_takej(maxj2,b, nib2, _G1+maxj1, _G2+maxj1, res,c);
    processOElems(false, b, curve, maxj1+maxj2, _G1, _G2, c);
  }
  else {
    int nib, maxj;
    int lib = internal_cp_to_keeper_var(curve, ib);
    maxj = _keepers[curve].convert2(lib, 2, nib, g1, g2, _G1, _G2);
    _keepers[curve]._takej(maxj, b, nib, _G1, _G2, res, c);
  }
  //debug_jac(b,curve,ib,4,g1,g2,res,c);
}

  /*
void MultiKeeper::take4(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c) {
  // SPEED: many repeat calls to private macros
  assert(ib < _mt->_crs[curve]._nPoints-1);

  if (ib==0 && !_mt->begWithRest(curve)) { 
    _keepers[curve].take2(b, internal_cp_to_keeper_var(curve, 1), g1+2, g2+2, res, c);
    internal_cp_to_keeper(curve, 0)->take2(b, internal_cp_to_keeper_var(curve, 0), g1, g2, res, c);
    //processOElems(true, b, curve, internal_cp_to_var(curve, 1), 
    //	  internal_cp_to_var(curve, 0), _mt->_crs[curve]._nVars, 4, g1, g2);
    processOElems(true, b, curve, 4, g1, g2, c);
  }

  else if (ib==_mt->_crs[curve]._nPoints-2 && !_mt->endWithRest(curve)) {
    _keepers[curve].take2(b, internal_cp_to_keeper_var(curve, ib), g1, g2, res, c);
    internal_cp_to_keeper(curve, ib+1)->take2(b, internal_cp_to_keeper_var(curve, ib+1), g1+2, g2+2, res, c);
    //processOElems(false, b, curve, internal_cp_to_var(curve, ib), 
    //	  internal_cp_to_var(curve, ib+1), _mt->_crs[curve]._nVars, 4, g1, g2);
    processOElems(false, b, curve, 4, g1, g2, c);
  }

  else {
    _keepers[curve].take4(b, internal_cp_to_keeper_var(curve, ib)  ,g1,g2,res,c);
    assert(_keepers+curve == internal_cp_to_keeper(curve, ib));
  }

  //debug_jac(b,curve,ib,4,g1,g2,res,c);
  }*/

void MultiKeeper::take4Key1(int curve, int ib, const double* g1, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  if (ib==0 && !_mt->begWithRest(curve)) { 
     int nib1, nib2, maxj1, maxj2;
     convert1( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g1);
     _keepers[curve]._takejKey1(maxj2, nib2, _G1+2, res, c);
     internal_cp_to_keeper(curve, 0)->_takejKey1(maxj1,nib1, _G1, res,c);
     processOElemsKey1(true, curve, maxj1+maxj2, _G1, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-2 && !_mt->endWithRest(curve)) {
     int nib1, nib2, maxj1, maxj2;
     convert1( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g1);
     _keepers[curve]._takejKey1(maxj1, nib1, _G1+2, res, c);
     internal_cp_to_keeper(curve, ib+1)->_takejKey1(maxj2,nib2, _G1+maxj1, res,c);
     processOElemsKey1(false, curve, maxj1+maxj2, _G1, c);
  }
  else {
    int nib, maxj;
    int lib = internal_cp_to_keeper_var(curve, ib);
    maxj = _keepers[curve].convert1(lib, 2, nib, g1, _G1);
    _keepers[curve]._takejKey1(maxj, nib, _G1, res, c);
  }
  //debug_jac(_mt->_numFrames-2,curve,ib,4,g1,res,c);  
}


/*
void MultiKeeper::take4Key1(int curve, int ib, const double* g1, double res, const double c) {
  assert(ib < _mt->_crs[curve]._nPoints-1);
  
  if (ib==0 && !_mt->begWithRest(curve)) { 
    _keepers[curve].take2Key1(internal_cp_to_keeper_var(curve, 1), g1+2, res, c);
    internal_cp_to_keeper(curve, 0)->take2Key1(internal_cp_to_keeper_var(curve, 0), g1, res, c);
    processOElemsKey1(true, curve, 4, g1, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-2 && !_mt->endWithRest(curve)) {
    _keepers[curve].take2Key1(internal_cp_to_keeper_var(curve, ib), g1, res, c);
    internal_cp_to_keeper(curve, ib+1)->take2Key1(internal_cp_to_keeper_var(curve, ib+1), g1+2, res, c);
    processOElemsKey1(false, curve, 4, g1, c);
  }
  else {
    _keepers[curve].take4Key1(internal_cp_to_keeper_var(curve, ib), g1,res,c);
    assert(_keepers+curve == internal_cp_to_keeper(curve, ib));
  }
  
  //debug_jac(_mt->_numFrames-2,curve,ib,4,g1,res,c);  
}


*/

void MultiKeeper::take6Key0(int curve, int ib, const double* g2, const double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  if (ib==0 && !_mt->begWithRest(curve)) { 
     int nib1, nib2, maxj1, maxj2;
     convert1( 1,2, ib, curve, nib1,nib2, maxj1, maxj2, g2); assert(maxj1==2);
     _keepers[curve]._takejKey0(maxj2, nib2, _G1+2, res, c);
     internal_cp_to_keeper(curve, 0)->_takejKey0(maxj1,nib1, _G1, res,c);
      processOElemsKey0(true, curve,1,0, maxj1+maxj2, _G1, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-3 && !_mt->endWithRest(curve)) {
     int nib1, nib2, maxj1, maxj2;
     convert1( 2,1, ib, curve, nib1,nib2, maxj1, maxj2, g2);  
     _keepers[curve]._takejKey0(maxj1, nib1, _G1, res, c);
     internal_cp_to_keeper(curve, ib+2)->_takejKey0(maxj2,nib2, _G1+maxj1, res,c);
     processOElemsKey0(false, curve, ib, ib+2, maxj1+maxj2, _G1, c);
  }
  else {
    int nib, maxj;
    int lib = internal_cp_to_keeper_var(curve, ib);
    maxj = _keepers[curve].convert1(lib, 3, nib, g2, _G1);
    _keepers[curve]._takejKey0(maxj, nib, _G1, res, c);
  }
  //debug_jac(0,curve,ib,6,g2,res,c);    
}
     
/*
void MultiKeeper::take6Key0(int curve, int ib, const double* g2, const double res, const double c) {
  if (ib==0 && !_mt->begWithRest(curve)) { 
    _keepers[curve].take4Key0(internal_cp_to_keeper_var(curve, 1), g2+2, res, c);
    internal_cp_to_keeper(curve, 0)->take2Key0(internal_cp_to_keeper_var(curve, 0), g2, res, c);
    processOElemsKey0(true, curve,1,0, 6, g2, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-3 && !_mt->endWithRest(curve)) {
    _keepers[curve].take4Key0(internal_cp_to_keeper_var(curve, ib), g2, res, c);
    internal_cp_to_keeper(curve, ib+2)->take2Key0(internal_cp_to_keeper_var(curve, ib+2), g2+4, res, c);
    processOElemsKey0(false, curve, ib, ib+2, 6, g2, c);
  }
  else {
    _keepers[curve].take6Key0(internal_cp_to_keeper_var(curve, ib), g2,res,c);
    //assert(curve == _mt->_crs[curve].varToCurve(ib)); 
  }

  //debug_jac(0,curve,ib,6,g2,res,c);    
  }*/



void MultiKeeper::take6(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  memset(_G2,0,6*sizeof(double));
  if (ib==0 && !_mt->begWithRest(curve)) { 
    int nib1, nib2, maxj1, maxj2;
    convert2( 1,2, ib, curve, nib1,nib2, maxj1, maxj2, g1,g2); assert(maxj1==2);
    _keepers[curve]._takej(maxj2, b, nib2, _G1+2, _G2+2, res, c);
    internal_cp_to_keeper(curve, 0)->_takej(maxj1,b, nib1, _G1, _G2, res,c);
    processOElems(true, b, curve, maxj1+maxj2, _G1, _G2, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-3 && !_mt->endWithRest(curve)) {
    int nib1, nib2, maxj1, maxj2;
    convert2( 2,1, ib, curve, nib1,nib2, maxj1, maxj2, g1,g2);
    _keepers[curve]._takej(maxj1, b, nib1, _G1, _G2, res, c);
    internal_cp_to_keeper(curve, ib+2)->_takej(maxj2,b, nib2, _G1+maxj1, _G2+maxj1, res,c);
    processOElems(false, b, curve, maxj1+maxj2, _G1, _G2, c);
  }
  else {
    int nib, maxj;
    int lib = internal_cp_to_keeper_var(curve, ib);
    maxj = _keepers[curve].convert2(lib, 3, nib, g1, g2, _G1, _G2);
    _keepers[curve]._takej(maxj, b, nib, _G1, _G2, res, c);
  }
  //debug_jac(b,curve,ib,6,g1,g2,res,c);
}

/*
void MultiKeeper::take6(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c) {
  // SPEED: many repeat calls to private macros
  assert(ib < _mt->_crs[curve]._nPoints-2);

  if (ib==0 && !_mt->begWithRest(curve)) { 
    _keepers[curve].take4(b, internal_cp_to_keeper_var(curve, 1), g1+2, g2+2, res, c);
    internal_cp_to_keeper(curve, 0)->take2(b, internal_cp_to_keeper_var(curve, 0), g1, g2, res, c);
    //processOElems(true, b, curve, internal_cp_to_var(curve, 1), 
    //	  internal_cp_to_var(curve, 0), _mt->_crs[curve]._nVars, 6, g1, g2);
    processOElems(true, b, curve, 6, g1, g2, c);
  }

  else if (ib==_mt->_crs[curve]._nPoints-3 && !_mt->endWithRest(curve)) {
    _keepers[curve].take4(b, internal_cp_to_keeper_var(curve, ib), g1, g2, res, c);
    internal_cp_to_keeper(curve, ib+2)->take2(b, internal_cp_to_keeper_var(curve, ib+2), g1+4, g2+4, res, c);
    //processOElems(false, b, curve, internal_cp_to_var(curve, ib), 
    //	  internal_cp_to_var(curve, ib+2), _mt->_crs[curve]._nVars, 6, g1, g2);
    processOElems(false, b, curve, 6, g1, g2, c);
  }

  else {
    _keepers[curve].take6(b, internal_cp_to_keeper_var(curve, ib)  ,g1,g2,res,c);
    //assert(curve == _mt->_crs[curve].varToCurve(ib)); 
  }

  //debug_jac(b,curve,ib,6,g1,g2,res,c);
}
*/

void MultiKeeper::take6Key1(int curve, int ib, const double* g1, double res, const double c) {
    memset(_G1,0,6*sizeof(double));
    if (ib==0 && !_mt->begWithRest(curve)) { 
      int nib1, nib2, maxj1, maxj2;
      convert1( 1,2, ib, curve, nib1,nib2, maxj1, maxj2, g1);  assert(maxj1==2);
      _keepers[curve]._takejKey1(maxj2, nib2, _G1+2, res, c);
      internal_cp_to_keeper(curve, 0)->_takejKey1(maxj1,nib1, _G1, res,c);
      processOElemsKey1(true, curve, maxj1+maxj2, _G1, c);
    }
    else if (ib==_mt->_crs[curve]._nPoints-3 && !_mt->endWithRest(curve)) {
      int nib1, nib2, maxj1, maxj2;
      convert1( 2,1, ib, curve, nib1,nib2,  maxj1, maxj2, g1);
      _keepers[curve]._takejKey1(maxj1, nib1, _G1, res, c);
      internal_cp_to_keeper(curve, ib+2)->_takejKey1(maxj2,nib2, _G1+maxj1, res,c);
      processOElemsKey1(false, curve, maxj1+maxj2, _G1, c);
    }
    else {
      int nib, maxj;
      int lib = internal_cp_to_keeper_var(curve, ib);
      maxj = _keepers[curve].convert1(lib, 3, nib, g1, _G1);
      _keepers[curve]._takejKey1(maxj, nib, _G1, res, c);
    }
    //debug_jac(_mt->_numFrames-2,curve,ib,6,g1,res,c);
}

/*
void MultiKeeper::take6Key1(int curve, int ib, const double* g1, double res, const double c) {
  if (ib==0 && !_mt->begWithRest(curve)) { 
    _keepers[curve].take4Key1(internal_cp_to_keeper_var(curve, 1), g1+2, res, c);
    internal_cp_to_keeper(curve, 0)->take2Key1(internal_cp_to_keeper_var(curve, 0), g1, res, c);
    processOElemsKey1(true, curve, 6, g1, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-3 && !_mt->endWithRest(curve)) {
    _keepers[curve].take4Key1(internal_cp_to_keeper_var(curve, ib), g1, res, c);
    internal_cp_to_keeper(curve, ib+2)->take2Key1(internal_cp_to_keeper_var(curve, ib+2), g1+4, res, c);
    processOElemsKey1(false, curve, 6, g1, c);
  }
  else {
    _keepers[curve].take6Key1(internal_cp_to_keeper_var(curve, ib), g1, res, c);
    //assert(curve == _mt->_crs[curve].varToCurve(ib)); 
  }

  //debug_jac(_mt->_numFrames-2,curve,ib,6,g1,res,c);
}
*/

void MultiKeeper::take4Special(int b, int curve, int ib, const double* g1, double res, const double c) {

  /*
  // Hessian is  c * (- a^2 (a^2 - 3r^2) / (a^2 + r^2)^3) * gg' (where g is column vector)
  // Gradient  c * (- a^2*r / (a^2 + r^2)^2) * g
  // a^2 = (1/16, or .0625) min_edge_at_keyframes

  double ch,cg, r2= res * res;
  double a2r2 =a2+r2;
  ch = - c * a2*(a2 - 3.*r2) / (a2r2*a2r2*a2r2);
  cg = - c * a2*res / (a2r2*a2r2);
  assert(!isnan(ch) && !isnan(cg));
  printf("r: %.5f r*r %.5f, p: %.5f, second: %.5f, first: %.5f\n", res, r2, -(r2/a2) / (1. + r2/a2),ch, cg);
  */

  memset(_G1,0,6*sizeof(double));
  if (ib==0 && !_mt->begWithRest(curve)) { 
    int nib1, nib2, maxj1, maxj2;
     convert1( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g1);
     _keepers[curve]._takejKeyb(b, maxj2, nib2, _G1+2, res, c);
     internal_cp_to_keeper(curve, 0)->_takejKeyb(b, maxj1,nib1, _G1, res, c);
     processOElemsKeyb(true, b, curve, maxj1+maxj2, _G1, c);
  }
  else if (ib==_mt->_crs[curve]._nPoints-2 && !_mt->endWithRest(curve)) {
     int nib1, nib2, maxj1, maxj2;
     convert1( 1,1, ib, curve, nib1,nib2, maxj1, maxj2, g1);
     _keepers[curve]._takejKeyb(b, maxj1, nib1, _G1+2, res, c);
     internal_cp_to_keeper(curve, ib+1)->_takejKeyb(b, maxj2,nib2, _G1+maxj1, res, c);
     processOElemsKeyb(false, b, curve, maxj1+maxj2, _G1, c);
  }
  else {
    int nib, maxj;
    int lib = internal_cp_to_keeper_var(curve, ib);
    maxj = _keepers[curve].convert1(lib, 2, nib, g1, _G1);
    _keepers[curve]._takejKeyb(b, maxj, nib, _G1, res, c);
  }
  


}


//-----------------------------------------------------------------------------------------------------------

OffElems* MultiKeeper::allocateElems(bool beg, int curve, int l, int k, int Ow, int Oh) {
  OffElems* oe = NULL;
  if (beg) {
    if (!_b_oelems[curve]) {
      oe = _b_oelems[curve] = new OffElems(_mt->_numFrames-1);
      oe->_l = l; oe->_k = k;  oe->_Ow = Ow; oe->_Oh = Oh;
    }
    else
      oe = _b_oelems[curve];
  }
  else {
    if (!_e_oelems[curve]) {
      oe = _e_oelems[curve] = new OffElems(_mt->_numFrames-1);
      oe->_l = l; oe->_k = k;  oe->_Ow = Ow; oe->_Oh = Oh;
    }
    else
      oe = _e_oelems[curve];
  }
  assert(oe);
  return oe;
}


// beg: beginning or end constraint
// size: 4 or 6 
// l: indice of first local variable (not div by 2 space)
// k: indice of non-local variable
//void MultiKeeper::processOElemsKey0(bool beg, int curve, int l, int k, int Ow, int Oh, int size, const double* g2, const double c) {
void MultiKeeper::processOElemsKey0(bool beg, int curve, int lcp, int kcp, int size, const double* g2, const double c) {
  int l = internal_cp_to_subvar(curve, lcp);                     //internal_cp_to_var(curve, lcp);
  int k = internal_cp_to_subvar(curve, kcp);                     //internal_cp_to_var(curve, kcp);
  int Ow = internal_cp_to_keeper(curve, kcp)->numVarOneBlock();   //2*_mt->cpToCurveNvars(curve, kcp);
  int Oh = internal_cp_to_keeper(curve, lcp)->numVarOneBlock();    //2*_mt->cpToCurveNvars(curve, lcp);


  OffElems* oe = allocateElems(beg, curve, l, k, Ow, Oh);

  if (size==6) {
    if (beg)
      outerProductPlus(g2+2, g2, 4, oe->_main, c);
    else
      outerProductPlus(g2, g2+4, 4, oe->_main, c);
  }
  else {
    assert(size==4);
    if (beg)
      outerProductPlus(g2+2, g2, 2, oe->_main, c);
    else
      outerProductPlus(g2, g2+2, 2, oe->_main + 4, c);
  }
  
}

// assume off elements already created
void MultiKeeper::processOElemsKey1(bool beg, int curve, int size, const double* g1, const double c) {
  OffElems* oe = NULL;
  if (beg)
    oe = _b_oelems[curve];
  else
    oe = _e_oelems[curve];
  assert(oe);
  
  int b = _mt->_numFrames-2;
  if (size==6) {
    if (beg)
      outerProductPlus(g1+2, g1, 4, oe->_main + 8*b, c);
    else
      outerProductPlus(g1, g1+4, 4, oe->_main + 8*b, c);
  }
  else {
    assert(size==4);
    if (beg)
      outerProductPlus(g1+2, g1, 2, oe->_main + 8*b, c);
    else
      outerProductPlus(g1, g1+2, 2, oe->_main + 8*b + 4, c);
  }

}

void MultiKeeper::processOElemsKeyb(bool beg, int b, int curve, int size, const double* g1, const double c) {
  OffElems* oe = NULL;
  oe = allocateElems(beg, curve, -1, -1, -1, -1);  // these numbers shouldn't matter for edgeKeep
  assert(oe);
  
  if (size==6) {
    if (beg)
      outerProductPlus(g1+2, g1, 4, oe->_main + 8*b, c);
    else
      outerProductPlus(g1, g1+4, 4, oe->_main + 8*b, c);
  }
  else {
    assert(size==4);
    if (beg)
      outerProductPlus(g1+2, g1, 2, oe->_main + 8*b, c);
    else
      outerProductPlus(g1, g1+2, 2, oe->_main + 8*b + 4, c);
  }

}


// assume off elements already created
void MultiKeeper::processOElems(bool beg, int b, int curve, int size, const double* g1, const double* g2, const double c) {
  OffElems* oe = NULL;
  if (beg)
    oe = _b_oelems[curve];
  else
    oe = _e_oelems[curve];
  assert(oe);


  if (size==6) {
    if (beg) {
      outerProductPlus(g1+2, g1, 4, oe->_main + 8*b, c);
      outerProductPlus(g1+2, g2, 4, oe->_right + 8*b, c);
      outerProductPlus(g2+2, g1, 4, oe->_left + 8*b, c);
      outerProductPlus(g2+2, g2, 4, oe->_main + 8*(b+1), c);
    }
    else {
      outerProductPlus(g1, g1+4, 4, oe->_main + 8*b, c);
      outerProductPlus(g1, g2+4, 4, oe->_right + 8*b, c);
      outerProductPlus(g2, g1+4, 4, oe->_left + 8*b, c);
      outerProductPlus(g2, g2+4, 4, oe->_main + 8*(b+1), c);
    }
  }
  else {
    assert(size==4);
    if (beg) {
      outerProductPlus(g1+2, g1, 2, oe->_main + 8*b, c);
      outerProductPlus(g1+2, g2, 2, oe->_right + 8*b, c);
      outerProductPlus(g2+2, g1, 2, oe->_left + 8*b, c);
      outerProductPlus(g2+2, g2, 2, oe->_main + 8*(b+1), c);
    }
    else {    // need + 4 since first 2x2 values are 0, then 2x2 numbers put in next 4
      outerProductPlus(g1, g1+2, 2, oe->_main + 8*b + 4, c);  
      outerProductPlus(g1, g2+2, 2, oe->_right + 8*b + 4, c);
      outerProductPlus(g2, g1+2, 2, oe->_left + 8*b + 4, c);
      outerProductPlus(g2, g2+2, 2, oe->_main + 8*(b+1) + 4, c);
    }
  }

  
}


// multiples dimensions An x 1  x  1 x 2, assume B has length 2, A has length An
// result is An x 2
void MultiKeeper::outerProductPlus(const double* A, const double* B, int An, double* putHere, const double co) {
  for (int j=0; j<An; ++j, putHere+=2, ++A) {
    (*putHere) += co * *A * B[0];
    (*(putHere+1)) += co * *A * B[1];
  }   
}


//-----------------------------------------------------------------------------------------------------------

// anum,bnum, ib,nib1,nib2 in /2 space
// sets maxj1, maxj2 in full 2 space
void MultiKeeper::convert2(int anum, int bnum, int ib, int curve, int& nib1, int& nib2, 
			  int& maxj1, int& maxj2, const double* g1, const double* g2) {
  int lib = internal_cp_to_keeper_var(curve, ib);
  maxj1 = internal_cp_to_keeper(curve, ib)->convert2(lib,anum,nib1,g1,g2,_G1,_G2);
  assert(maxj1<=anum*2);
  lib = internal_cp_to_keeper_var(curve, ib+anum);
  maxj2 = internal_cp_to_keeper(curve, ib+anum)->convert2(lib,bnum,nib2,g1 + 2*anum, g2 + 2*anum,
							  _G1 + maxj1,_G2 + maxj1);
}

void MultiKeeper::convert1(int anum, int bnum, int ib, int curve, int& nib1, int& nib2, 
			  int& maxj1, int& maxj2, const double* g1) {
  int lib = internal_cp_to_keeper_var(curve, ib);
  maxj1 = internal_cp_to_keeper(curve, ib)->convert1(lib,anum,nib1,g1,_G1);
  assert(maxj1<=anum*2);
  lib = internal_cp_to_keeper_var(curve, ib+anum);
  maxj2 = internal_cp_to_keeper(curve, ib+anum)->convert1(lib,bnum,nib2,g1 + 2*anum,_G1 + maxj1);
}

//-----------------------------------------------------------------------------------------------------------

void MultiKeeper::fprint(FILE* fp) {
  int tnvar = _subsNvar*2, i, j,c;
  
  double *mat = new double[tnvar*tnvar], *Z=mat;
  memset(mat, 0, tnvar*tnvar*sizeof(double));
  

  for (c=0; c<_mt->_nCurves; ++c) {
    // on-diagonal elements
    _keepers[c].fillMatrix(Z, tnvar);
    Z += (tnvar+1)*_keepers[c].numVar();
    
    // off    
    if (_b_oelems[c])
      _b_oelems[c]->fillMatrix(mat, tnvar, _mt->_numFrames-1);
    if (_e_oelems[c])
      _e_oelems[c]->fillMatrix(mat, tnvar, _mt->_numFrames-1);
    
  }    

  for (j=0,c=0; j<tnvar; ++j) {
    for (i=0; i<tnvar; ++i,++c)
      fprintf(fp,"%.10f ", mat[c]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  delete[] mat;
}


void OffElems::fillMatrix(double* mat, int stride, int nf) const {
  int t,a,j,i,r, c;

  // main strip
  for (t=0, a=0; t<nf; ++t) {
    for (j=0; j<4; ++j) {
      r = _Oh*t + _l + j;
      for (i=0; i<2; ++i, ++a) {
	c = t*_Ow + _k + i;
	mat[r*stride+c] = _main[a];
	mat[c*stride+r] = _main[a];
      }
    }
  }

  //right strip
  for (t=0, a=0; t<nf-1; ++t) {
    for (j=0; j<4; ++j) {
      r = _Oh*t + _l + j;
      for (i=0; i<2; ++i, ++a) {
	c = (t+1)*_Ow + _k + i;
	mat[r*stride+c] = _right[a];
	mat[c*stride+r] = _right[a];
      }
    }
  }

  //left strip
  for (t=0, a=0; t<nf-1; ++t) {
    for (j=0; j<4; ++j) {
      r = _Oh*(1+t) + _l + j;
      for (i=0; i<2; ++i, ++a) {
	c = t*_Ow + _k + i;
	mat[r*stride+c] = _left[a];
	mat[c*stride+r] = _left[a];
      }
    }
  }
  
}




void MultiKeeper::printCurves(const Vec2f* Z) {
  for (int c=0; c<_mt->_nCurves; ++c) {
    printf("CURVE %d\n",c);
    for (int t=0; t<=_mt->_numFrames; ++t) {
      printf("  TIME %d\n",t);
      for (int n=0; n<_mt->_crs[c]._nPoints; ++n) {
	int l = _mt->packed_index(c,n,t);
	printf("    %.3f %.3f\n",Z[l].x(), Z[l].y());
      }
    }
  }
  
}


void MultiKeeper::printSubs(FILE* fp) {
  int w=0, h=0;
  for (int c=0; c<_mt->_nCurves; ++c)
    for (int t=0; t<_mt->_numFrames-1; ++t, h += _keepers[c].numVarOneBlock(), w += _keepers[c].origNVOneBlock())
      _keepers[c].printSubs2(fp, w, h);
      
}


bool MultiKeeper::precond() const {
  return _precond;
}

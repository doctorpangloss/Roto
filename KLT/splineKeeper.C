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



//#include <iostream>
#include "splineKeeper.h"
#include "linearSolver.h"  // not forever SPEED 1
/*
#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
*/

// SPEED 1: overarching, both kltSpline and here, replace double* work with ublas vectors
#include "myassert.h"


SplineKeeper::SplineKeeper(MultiSplineData* mts) : _mts(mts) {  
  //  _fp = fopen("J.dat","w"); row=1;
  mutualInit();

#ifdef _CHECKJ_
  _jcheck = new BuildingSplineKeeper(_nv);
#endif


  _precond = true;   
  if (precond()) {
    vector<int>* fixers = new vector<int>[mts->_nCurves];
    int* numvars = new int[mts->_nCurves];
    
    const FixedControlV& flocs = _mts->getFixedLocs();
    for (FixedControlV::const_iterator j = flocs.begin(); j!= flocs.end(); ++j) {
      fixers[j->_c].push_back(2*_mts->tcn_var(j->_t, j->_c, j->_n));
    }
    
    for (int i=0; i<_mts->_nCurves; ++i) 
      numvars[i] = _mts->c_numVarsOneFrame(i)/2;
    
    _sweeper = new HB_Sweep(mts->_nCurves, numvars, mts->_numFrames-1, fixers);
    delete[] fixers;
    delete[] numvars;
    
  }
  else 
    _sweeper = NULL;

}




SplineKeeper::SplineKeeper(SplineKeeper* other) {  
  //  _fp = fopen("J.dat","w"); row=1;
  _mts = other->_mts;
  mutualInit();
  _sweeper = NULL;
#ifdef _CHECKJ_
  _jcheck = new BuildingSplineKeeper(other->_jcheck);
#endif
}

void SplineKeeper::mutualInit() {
  _nv = _mts->numVars();

  constructMultiDiagMatrix();

  //_mat = new MyMat(_nv, _nv);  // MAT

  //int maxWidth = 2*(_mts->tcn_var(2,0,1) - 1) + 14;  
  //_mat = new MyMat(_nv, maxWidth);

  _g = new double[_nv]; //Ubv(_nv);  
  //*_g =Zero_Ubv(_nv);
  memset(_g, 0, _nv*sizeof(double));
  //_crapg = new double[_nv]; 
 //_crapgDirty = true;
  _precond = false;
}

// BUG!
void SplineKeeper::constructMultiDiagMatrix() {
  int numCurves = _mts->_nCurves;
  std::vector<int> maxWidths, dims;
  maxWidths.reserve(numCurves); dims.reserve(numCurves);
  
  for (uint i = 0; i<(uint) numCurves; ++i) {
    dims.push_back(_mts->c_numVars(i));
    if (_mts->_numFrames!=2)
      maxWidths.push_back( 2*( _mts->tcn_var(2,i,1) - _mts->tcn_var(1,i,1)  ) + 14 );
    else
      maxWidths.push_back( 100000); // G!
  }
  _mat = new MultiDiagMatrix(numCurves, maxWidths, dims);  // MAT
}

double SplineKeeper::calculateRo(const double* p, const double newTheta, double const oldTheta) const { 
  double num = oldTheta - newTheta;
  double denom = vecDot(_nv, g(), p);
  double *Bp = new double[_nv];
  matVecMult(p, Bp);
  denom += .5 * vecDot(_nv, p, Bp);
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


void SplineKeeper::take8(int vars[4], double* J, double r, double c) {  
  int i, j;

  //_crapgDirty = true;
  MyMat& mat = *_mat;  //Ubv& g = *_g;
  double f;
  int iv, jv;
  for (i=0,iv=0; i<8; i+=2,++iv) {
    assert(vars[iv]+1<_nv);

    mat.takeUnrdF(vars[iv], vars[iv],     c*J[i]*J[i] );
    mat.takeUnrdF(vars[iv]+1, vars[iv]+1, c*J[i+1]*J[i+1] );

    f = c*J[i]*J[i+1];
    mat.takeUnrdF(vars[iv], vars[iv]+1, f);
    //mat(vars[iv]+1, vars[iv])   += f;  // MAT

    _g[vars[iv]]   += c*J[i]*r;
    _g[vars[iv]+1] += c*J[i+1]*r;
    assert(!isnan(_g[vars[iv]]));
    assert(!isnan(_g[vars[iv]+1]));

    for (j=i+2, jv=iv+1; j<8; j+=2, ++jv) {
      assert(vars[jv]+1<_nv);

      f = c*J[i] * J[j];
      mat.takeUnrdF(vars[iv], vars[jv], f);
      //mat(vars[jv], vars[iv])     += f;  // MAT

      f = c*J[i] * J[j+1];
      mat.takeUnrdF(vars[iv], vars[jv]+1, f);
      //mat(vars[jv]+1, vars[iv])   += f;  // MAT

      f = c*J[i+1] * J[j];
      mat.takeUnrdF(vars[iv]+1, vars[jv], f);
      //mat(vars[jv], vars[iv]+1)   += f;  // MAT

      f = c*J[i+1] * J[j+1];
      mat.takeUnrdF(vars[iv]+1, vars[jv]+1, f);
      //mat(vars[jv]+1, vars[iv]+1) += f;  // MAT

    }
  }
}

void SplineKeeper::take8Key(int vars[4], double* J, double r1, double r2, double r3) {  
  int i, j;

  //_crapgDirty = true;
  MyMat& mat = *_mat;  //Ubv& g = *_g;
  double f;
  int iv, jv;
  for (i=0,iv=0; i<8; i+=2,++iv) {
    assert(vars[iv]+1<_nv);

    mat.takeUnrdF(vars[iv], vars[iv],     J[i]*J[i] + J[8+i]*J[8+i]  + J[16+i]*J[16+i]  );
    mat.takeUnrdF(vars[iv]+1, vars[iv]+1, J[i+1]*J[i+1] + J[9+i]*J[9+i]  + J[17+i]*J[17+i]  );

    f = J[i]*J[i+1] + J[8+i]*J[9+i] + J[16+i]*J[17+i];
    mat.takeUnrdF(vars[iv], vars[iv]+1, f);
    //mat(vars[iv]+1, vars[iv])   += f;  // MAT

    _g[vars[iv]]   += J[i]*r1 + J[8+i]*r2 + J[16+i]*r3;
    _g[vars[iv]+1] += J[i+1]*r1 + J[9+i]*r2 + J[17+i]*r3;
    assert(!isnan(_g[vars[iv]]));
    assert(!isnan(_g[vars[iv]+1]));

    for (j=i+2, jv=iv+1; j<8; j+=2, ++jv) {
      assert(vars[jv]+1<_nv);

      f = J[i] * J[j] + J[8+i] * J[8+j] + J[16+i] * J[16+j];
      mat.takeUnrdF(vars[iv], vars[jv], f);
      //mat(vars[jv], vars[iv])     += f;  // MAT

      f = J[i] * J[j+1] + J[8+i] * J[9+j] + J[16+i] * J[17+j];
      mat.takeUnrdF(vars[iv], vars[jv]+1, f);
      //mat(vars[jv]+1, vars[iv])   += f;  // MAT

      f = J[i+1] * J[j] + J[9+i] * J[8+j] + J[17+i] * J[16+j];
      mat.takeUnrdF(vars[iv]+1, vars[jv], f);
      //mat(vars[jv], vars[iv]+1)   += f;  // MAT

      f = J[i+1] * J[j+1] + J[i+9] * J[j+9] + J[i+17] * J[j+17];
      mat.takeUnrdF(vars[iv]+1, vars[jv]+1, f);
      //mat(vars[jv]+1, vars[iv]+1) += f;  // MAT

    }
  }

}

// same, but two sets of J, var, for 4 blocks
// (J = concat J1, J2);
void SplineKeeper::take88(int vars1[4], int vars2[4], double* J1, double* J2, double r) {  
  int i, j;

  //_crapgDirty = true;
  memcpy(_vars, vars1, 4*sizeof(int));    memcpy(_vars+4, vars2, 4*sizeof(int));
  memcpy(_J, J1, 8*sizeof(double));       memcpy(_J+8, J2, 8*sizeof(double));
  MyMat& mat = *_mat;  //Ubv& g = *_g;
  double f;
  int iv, jv;
  for (i=0,iv=0; i<16; i+=2,++iv) {
    assert(_vars[iv]+1<_nv);

    mat.takeUnrdF(_vars[iv], _vars[iv],     _J[i]*_J[i]);
    mat.takeUnrdF(_vars[iv]+1, _vars[iv]+1, _J[i+1]*_J[i+1]);

    f = _J[i]*_J[i+1];
    mat.takeUnrdF(_vars[iv], _vars[iv]+1, f);
    //mat(_vars[iv]+1, _vars[iv])   += f;  // MAT

    _g[_vars[iv]] += _J[i]*r;
    _g[_vars[iv]+1] += _J[i+1]*r;
    assert(!isnan(_g[_vars[iv]]));
    assert(!isnan(_g[_vars[iv]+1]));

    for (j=i+2, jv=iv+1; j<16; j+=2, ++jv) {
      assert(_vars[jv]+1<_nv);

      f = _J[i] * _J[j];
      mat.takeUnrdF(_vars[iv], _vars[jv], f);
      //mat(_vars[jv], _vars[iv])     += f;  // MAT

      f = _J[i] * _J[j+1];
      mat.takeUnrdF(_vars[iv], _vars[jv]+1, f);
      //mat(_vars[jv]+1, _vars[iv])   += f;  // MAT

      f = _J[i+1] * _J[j];
      mat.takeUnrdF(_vars[iv]+1, _vars[jv], f);
      //mat(_vars[jv], _vars[iv]+1)   += f; // MAT


      f = _J[i+1] * _J[j+1];
      mat.takeUnrdF(_vars[iv]+1, _vars[jv]+1, f);
      //mat(_vars[jv]+1, _vars[iv]+1) += f;  // MAT
    }
  } 
}


// J must be of length numVar, J'J added to A and J'r added to b
void SplineKeeper::takeJVector(const Ubcv& J, double r) {  
  //_crapgDirty = true;
  Ubcv::const_iterator i1, i2, jend;
  int num = 0;
  double f;
  MyMat& mat = *_mat;  //Ubv& g = *_g;

  jend = J.end();
  for (i1=J.begin(); i1!=jend; ++i1, ++num) {
    _g[i1.index()] += *i1 * r;
    assert(!isnan(_g[i1.index()]));
    //fprintf(_fp, "%d %d %.10f\n",row+1, i1.index()+1, *i1);  
    for (i2=J.begin(); i2!=jend; ++i2) {

      if (i2.index()<i1.index()) continue; // below diagonal

      f = *i2 * *i1;
      if (f!=0) {  // SPEED: should time effect of this
	
	//mat(i1.index(), i2.index()) += f; 
	mat.takeF(i1.index(), i2.index(), f);
      }
      //else
      //printf("%d %d %.10e\n",i1.index(), i2.index(), f);
      
      //if (i1.index() != i2.index())  // MAT
      //mat(i2.index(), i1.index()) += f;

    }
  }
  //assert(num < 33);
  //++row;
  //if (num>8)
  //printf("Took %d jacobian, %f residual\n",num,r);

#ifdef _CHECKJ_
  _jcheck->takeJVector(J,r);
#endif
}





// SPEED: eventually propagate uBLAS types up farther for efficiency
//r = Ax;
void SplineKeeper::matVecMult(const double x[], double r[]) const { 
  //int i;
  //Ubv xvec(_nv);  // MAT
  //Ubv rvec(_nv);
  //for (i=0; i<_nv; ++i)  // MAT
  //xvec(i) = x[i];

  //axpy_prod(*_mat, xvec, rvec);  // MAT
  _mat->matVecMult(x,r);

  //for (i=0; i<_nv; ++i) // MAT
  //r[i] = rvec(i);
}


void SplineKeeper::refresh() { 
  _mat->clear();
  //*_g =Zero_Ubv(_nv);
  memset(_g, 0, _nv*sizeof(double));
  //_crapgDirty = true;

#ifdef _CHECKJ_
  _jcheck->refresh();
#endif
}


void SplineKeeper::shit() {
  //if (_fp) {fclose(_fp); _fp = NULL;}

  /*
  // y = Ax;
  double* x = new double[_nv];
  double* y = new double[_nv];
  int i;
  for (i=0; i<_nv; ++i)
    x[i] = 1.;
  matVecMult(x, y); 
  FILE* fp = fopen("y.dat", "w");
  for (i=0; i<_nv; ++i)
    fprintf(fp, "%.10e\n", y[i]);
  fclose(fp);
    delete[] x; delete[] y;
  */

  outputg("g.txt");
  outputMat("B.dat"); 
#ifdef _CHECKJ_
  _jcheck->outputJ("BJ.dat");
#endif
}

void SplineKeeper::outputg(char* name) const {
  
  FILE* fp = fopen(name, "w");
  for (int i=0; i<_nv; ++i)
    fprintf(fp, "%.10e\n", _g[i]);
  fclose(fp);
}

void SplineKeeper::outputMat(char* name) {
  assert(_mat);
  
  FILE* fp = fopen(name, "w");
  _mat->outputMat(fp);  // MAT
  
  /*
  for (MyMat::iterator1 it1 = _mat->begin1(); it1 !=_mat->end1(); ++it1)
    for (MyMat::iterator2 it2 = it1.begin(); it2!=it1.end(); ++it2) {
      double e = *it2;
      fprintf(fp, "%d %d %.10f\n",it2.index1()+1, it2.index2()+1, e);
    }
  */
  fclose(fp);

}

void SplineKeeper::scalarMult(const double m) {

  //(*_g) *= m;
  for (int i=0; i<_nv; ++i) {
    assert(!isnan(_g[i]));
    _g[i] *= m;
    assert(!isnan(_g[i]));
  }

  (*_mat) *= m;

  //_crapgDirty = true;

#ifdef _CHECKJ_
  _jcheck->scalarMult(m);
#endif
}

void SplineKeeper::addMat(const SplineKeeper* sp) {
  //(*_g) += *(sp->_g);
  for (int i=0; i<_nv; ++i) {
    _g[i] += sp->_g[i];
    assert(!isnan(_g[i]));
  }

  (*_mat) += *(sp->_mat);
  //_crapgDirty = true;

#ifdef _CHECKJ_
  _jcheck->addMat(sp->_jcheck);
#endif
}


/*
const double* SplineKeeper::neg_g() {
  if (_crapgDirty)
    g();

  for (int i=0; i<_nv; ++i)
    _crapg[i] = -_crapg[i];

  _crapgDirty = true;
  return _crapg;
  }*/

const double* SplineKeeper::g() const {
  /*if (!_crapgDirty) return _crapg;

  for (int i=0; i<_nv; ++i) {
    _crapg[i] = (*_g)[i];
    assert(!isnan(_crapg[i]));
  }

  return _crapg;*/
  for (int i=0; i<_nv; ++i) // G!
    assert(!isnan(_g[i]));

  return _g;
}

void SplineKeeper::handleConstraints(double* r) {
  int index;
  if (!precond()) {
    const FixedControlV& V = _mts->getFixedLocs();
    for (FixedControlV::const_iterator i = V.begin(); i!=V.end(); ++i) {
      index = 2*_mts->tcn_var(i->_t, i->_c, i->_n);
      r[index] = r[index+1] = 0;
    }
  }
  else
    _sweeper->zeroFixed(r);
}

bool SplineKeeper::precond() const {
  return _precond;
}


SplineKeeper::~SplineKeeper() {
  if (_mat)
    delete _mat;
  delete[] _g;
  //delete[] _crapg;
  if (_sweeper) delete _sweeper;

#ifdef _CHECKJ_
  delete _jcheck;
#endif
}



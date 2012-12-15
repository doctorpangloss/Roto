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



#ifndef SPLINEKEEPER_H
#define SPLINEKEEPER_H

#include <iostream>
#include <stdio.h>
#include "genKeeper.h"
#include "multiSplineData.h"
#include "HB_Sweep.h"


#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/ublas/vector_of_vector.hpp>
//#include <boost/numeric/ublas/symmetric.hpp>


#include "diagMatrix.h"
#include "multiDiagMatrix.h"


//#define _CHECKJ_
#ifdef _CHECKJ_
#include "buildingSplineKeeper.h"
#endif

typedef boost::numeric::ublas::compressed_matrix<double> Ubm;
//typedef boost::numeric::ublas::generalized_vector_of_vector< double, boost::numeric::ublas::row_major, boost::numeric::ublas::vector<boost::numeric::ublas::coordinate_vector<double> >  > Ubm;
typedef boost::numeric::ublas::vector<double> Ubv;
typedef boost::numeric::ublas::zero_vector<double> Zero_Ubv;
typedef boost::numeric::ublas::compressed_vector<double> Ubcv;
//typedef boost::numeric::ublas::symmetric_adaptor<Ubm, boost::numeric::ublas::lower> Ubsm;

typedef MultiDiagMatrix MyMat;  // MAT 

class SplineKeeper : public GenKeeper {


 public:

  
  SplineKeeper(MultiSplineData* mts);

  SplineKeeper(SplineKeeper* other);

  void mutualInit();

  ~SplineKeeper();

  int numVar() const {return _nv; }  // full-2

  double calculateRo(const double* p, const double newTheta, double const oldTheta) const;



  void refresh(); 

  void handleConstraints(double* r);
  
  bool precond() const;
  void precondVec(double* v) {
    _sweeper->sweepUp(v); _sweeper->sweepDown(v);
  }
  

  // Variable numbers in vars (full-2, but each represents a pair), 3 jacobians in J[24], 3 residual r's,
  // adds J'J, J'r to system
  void take8Key(int vars[4], double* J, double r1, double r2, double r3);
  // takes one residual, also one c so that we add cJ'J and cJ'r;
  void take8(int vars[4], double* J, double r, double c); 

  // same, but two sets of J, var, for 4 blocks
  // (J = concat J1, J2);
  void take88(int vars1[4], int vars2[4], double* J1, double* J2, double r);

  // J must be of length numVar, J'J added to A and J'r added to b
  void takeJVector(const Ubcv& J, double r);

  void takeF(const int i, const int j, const double f);
  void takeg(const int i, const double f);

  void matVecMult(const double x[], double r[]) const;

  void scalarMult(const double m);
  void addMat(const SplineKeeper* sp);

  void outputMat(char* name);

  const double* g() const;
  void outputg(char* name) const;
  //const double* neg_g();

  void shit();

 private:

  void constructMultiDiagMatrix();


  MultiSplineData* _mts;
  int _nv; // full-2
  FILE* _fp;
  int row;
  MyMat* _mat;
  //Ubsm* _mats;
  double* _g; // or r
  //mutable double* _crapg; 
  //mutable bool _crapgDirty;
  int _vars[8];   // temp
  double _J[16];  // temp
  bool _precond;
  HB_Sweep* _sweeper;

#ifdef _CHECKJ_
  BuildingSplineKeeper* _jcheck;
#endif
};

inline void SplineKeeper::takeF(const int i, const int j, const double f) {
  assert(_mat);
  _mat->takeUnrdF(i,j,f);
}

inline void SplineKeeper::takeg(const int i, const double f) {
  assert(_g);
  _g[i] += f;
  assert(!isnan(_g[i]));
}


#endif

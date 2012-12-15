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



#ifndef MULTIKEEPER_H
#define MULTIKEEPER_H


#include "keeper.h"
#include "multiTrackData.h"
#include "HB_Sweep.h"
//#include "debugJac.h"


class OffElems {
 public:

  OffElems(int nFrames) {
    _main = new double[8*nFrames]; 
    memset(_main, 0, 8*nFrames*sizeof(double));
    _left = new double[8*(nFrames-1)];
    memset(_left, 0, 8*(nFrames-1)*sizeof(double));
    _right = new double[8*(nFrames-1)];
    memset(_right, 0, 8*(nFrames-1)*sizeof(double));
    _nFrames = nFrames;
  }

  OffElems(OffElems* other) {
    _nFrames = other->_nFrames;
    _l = other->_l; _k = other->_k;
    _Ow = other->_Ow; _Oh = other->_Oh;
    _main = new double[8*_nFrames]; 
    memset(_main, 0, 8*_nFrames*sizeof(double));
    _left = new double[8*(_nFrames-1)];
    memset(_left, 0, 8*(_nFrames-1)*sizeof(double));
    _right = new double[8*(_nFrames-1)];
    memset(_right, 0, 8*(_nFrames-1)*sizeof(double));
  }

  void refresh() {
    memset(_main, 0, 8*_nFrames*sizeof(double));
    memset(_left, 0, 8*(_nFrames-1)*sizeof(double));
    memset(_right, 0, 8*(_nFrames-1)*sizeof(double));
  }

  void fillMatrix(double* mat, int stride, int nf) const;

  void matVecMult(const double x[], double r[]) const;

  void printSparse(FILE* cp) const;

  void add(OffElems* other, double c);

  ~OffElems() { delete[] _main; delete[] _left; delete[] _right; }
  
  double *_main, *_left, *_right;
  
  int _l, _k, _Ow, _Oh; // _l and _k are for first frame,  to traverse frames
  int _nFrames;
  
 private:
  OffElems();
};

//typedef vector<OffElems*> OffElemsVec;

class MultiKeeper : public GenKeeper {

 public:

  MultiKeeper(const MultiTrackData* mt, int curveSampling, int hold, bool precond = TRUE);
  MultiKeeper(MultiKeeper* other); // careful!  Needs other to live while this lives

  void matVecMult(const double x[], double r[]) const ;  // r = B*x
  
  int numVar() const { return _subsNvar*2; } // full2 space
  
  void refresh(); // set B, g to zero
  
  const double* g() const;
  
  double calculateRo(const double* p, const double newTheta, double const oldTheta) const;
    
  double createTestSol(const Vec2f* Z, const double* sol, Vec2f* Z2) const;

  void handleConstraints(double* r); // set components of r to 0 if they can't move

  bool precond() const;
  void precondVec(double* v) {
    _sweeper->sweepUp(v); _sweeper->sweepDown(v);
  }

  void add( MultiKeeper* other, double c );
  
  ~MultiKeeper();
  
  // Note that this ib is in /2 space!
  void take2Key0(int curve, int ib, const double* g2, double res, const double c);
  void take2(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c);
  void take2Key1(int curve, int ib, const double* g1, double res, const double c);
  
  void take4Key0(int curve, int ib, const double* g2, double res, const double c);
  void take4(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c);
  void take4Key1(int curve, int ib, const double* g1, double res, const double c);

  void take6Key0(int curve, int ib, const double* g2, const double res, const double c);
  void take6(int b, int curve, int ib, const double* g1, const double* g2, double res, const double c);
  void take6Key1(int curve, int ib, const double* g1, double res, const double c);

  void take4Special(int b, int curve, int ib, const double* g1, double res, const double c);

  void fprint(FILE* fp);
  void printSparse(FILE* fp);
  void outputg();
  void outputSol(double* x);
  void printCurves(const Vec2f* Z);
  void printSubs(FILE* fp);

  //DebugJacobians djacs; 

 private:
  // full 2 output, /2 input
  int internal_cp_to_keeper_var(const int curve, const int p) {
    return (2 * _mt->_crs[curve].varToCurveVar(p));
    //return (2 * _cp_varInKeeper[_mt->_curve_map[curve] + p]); } // SPEED: use shift
  }
  Keeper* internal_cp_to_keeper(const int curve, const int p) {
    return (_keepers +  _mt->_crs[curve].varToCurve(p));
    //return (_cp_keeper[_mt->_curve_map[curve] + p]);
  }
  int internal_cp_to_var(const int curve, const int p) {
    return (2 * _mt->_crs[curve].varToGlobalVar(p));
    //return (2 * _mt->_cp_varMap[_mt->_curve_map[curve] + p]);
  }

  //result is global subvar, +1 could also be affected
  int internal_cp_to_subvar(const int curve, const int p) {
    int keeperNum = _mt->_crs[curve].varToCurve(p);
    int keeperVar = _mt->_crs[curve].varToCurveVar(p);
    return _sub_offset[keeperNum] + 2*_keepers[keeperNum].subVar(keeperVar);
  }

  //void assignVarsKeepers(const int curveSampling);
  OffElems* allocateElems(bool beg, int curve, int l, int k, int Ow, int Oh);
  void processOElemsKey0(bool beg, int curve, int lcp, int kcp, int size, const double* g2, const double c);
  void processOElemsKey1(bool beg, int curve, int size, const double* g1, const double c);
  void processOElemsKeyb(bool beg, int b, int curve, int size, const double* g1, const double c);
  void processOElems(bool beg, int b, int curve, int size, const double* g1, const double* g2, const double c);
  void outerProductPlus(const double* A, const double* B, int An, double* putHere, const double co);
  
  void convert2(int anum, int bnum, int ib, int curve, int& nib1, int& nib2, 
		int& maxj1, int& maxj2, const double* g1, const double* g2);
  void convert1(int anum, int bnum, int ib, int curve, int& nib1, int& nib2, 
		int& maxj1, int& maxj2, const double* g1);

  //void debug_jac(int t, int c, int cp, int size, const double* g, double res, double cons); 
  //void debug_jac(int t, int c, int cp, int size, const double* g1, const double* g2, double res, double cons); 

  Keeper* _keepers;
  OffElems **_b_oelems, **_e_oelems;
  const MultiTrackData* _mt;
  HB_Sweep* _sweeper;

  int _subsNvar; // number of subsampled variables, total, /2 space
  int* _sub_offset;  // for each curve/keeper, starting global sub variable, full 2 space

  mutable double* _g;
  mutable bool _gDirty;

  bool _precond;

  double _G1[6], _G2[6];




};







#endif

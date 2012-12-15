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



#ifndef KEEPER_H
#define KEEPER_H

#include <assert.h>
#include <stdio.h>
#include "linearSolver.h"
#include "../jl_vectors.h"
#include "genKeeper.h"

#define _WIDTH 13
#define _DIAG 6
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#ifdef __APPLE__
#define _D _D
#define _P _P
#endif




struct SubRcd {

  void set(const double k, const int P) { _k=k; _P=P; }

  double _k;  // _k * _P + (1-k) * _P_(+1)
  int _P;

};

class Keeper : public GenKeeper{

 public:

  Keeper() {_subs = NULL;}
  Keeper(int num_blocks, int block_size, int curveSampling);
  void init(int num_blocks, int block_size, int curveSampling, int hold);
  void init(Keeper* other);

  void matVecMult(const double x[], double r[]) const ;  // r = B*x
  
  int numVar() const { return _numVar; }  // this is subsampled, full 2, coordinate system
  int numVarOneBlock() const { return _block_size; } // subsampled, full 2, one frame
  int origNVOneBlock() const { return _o_bs;} // full 2, one frame

  void refresh(); // set B, g to zero

  const double* g() const { return _g; }

  double calculateRo(const double* p, const double newTheta, double const oldTheta) const;

  // i is a /2 space curve var.  Returns FIRST sub variable (/2 space) this curve var influences.
  int subVar(const int i) const {
    return _subs[i]._P;
  }

  double createTestSol(const Vec2f* Z, const double* sol, Vec2f* Z2) const;
  
  void print();
  void printSparse(int off, FILE* fp);
  void printFull(FILE* fp);
  
  void printSubs(FILE* fp);
  void printSubs2(FILE* fp, int w, int h); // sparse
  void printE(FILE* fp) const;
  void fillMatrix(double* mat, int stride) const;

  void add(Keeper* other, double c); 
  
  ~Keeper();

  // these are in original coordinate system
  // /2 space
  //void fill_e(double* e); 
  
  void take2Key0(int ib, const double* g2, double res, const double c);
  void take2(int b, int ib, const double* g1, const double* g2, double res, const double c);
  void take2Key1(int ib, const double* g1, double res, const double c);
  
  void take4Key0(int ib, const double* g2, double res, const double c);
  void take4(int b, int ib, const double* g1, const double* g2, double res, const double c);
  void take4Key1(int ib, const double* g1, double res, const double c);

  void take6Key0(int ib, const double* g2, const double res, const double c);
  void take6(int b, int ib, const double* g1, const double* g2, double res, const double c);
  void take6Key1(int ib, const double* g1, double res, const double c);

  //void take4Keyb(int b, int ib, const double* g2, const double ch, const double cg);

  // these are in subsampled coordinate system
  // SPEED: get rid of g1,g2, use _G1, _G2
  // /2 space
  
  void _take2Key0(const int ib, const double* g2, double res, const double c);
  void _take2(const int b, const int ib, const double* g1, const double* g2, double res, const double c);
  void _take2Key1(const int ib, const double* g1, double res, const double c);
  
  void _take4Key0(const int ib, const double* g2, double res, const double c);
  void _take4(const int b, const int ib, const double* g1, const double* g2, double res, const double c);
  void _take4Key1(const int ib, const double* g1, double res, const double c);

  void _take6Key0(const int ib, const double* g2, const double res, const double c);
  void _take6(const int b, const int ib, const double* g1, const double* g2, double res, const double c);
  void _take6Key1(const int ib, const double* g1, double res, const double c);

  void _takejKey0(int j, const int ib, const double* g2, double res, const double c);
  void _takej(int j, const int b, const int ib, const double* g1, const double* g2, double res, const double c);
  void _takejKey1(int j, const int ib, const double* g1, double res, const double c);

  void _takejKeyb(const int b, int j, const int ib, const double* g1, double res, const double c);

  void _take2Keyb(const int b, const int ib, const double* g1, double res, const double c);
  void _take4Keyb(const int b, const int ib, const double* g1, double res, const double c);

  int convert1(int ib, const int hm, int& nib, const double* g1, double* G1);
  int convert2(int ib, const int hm, int& nib, const double* g1, const double* g2, double* G1, double* G2);

 protected:
  int convert1(int ib, const int hm, int& nib, const double* g1);
  int convert2(int ib, const int hm, int& nib, const double* g1, const double* g2);
  void make2Block(const double* a, const double* b, const double c);
  void add2Block(double* strip, const int row);
  void make4Block(const double* a, const double* b, const double c);
  void add4Block(double* strip, const int row);
  void make6Block(const double* a, const double* b, const double c);
  void add6Block(double* strip, const int row);

  // data
  // these are in full2 coordinates
  int _numVar, _offNumVar, _cs;
  int _num_blocks, _block_size, _o_bs;
  double *_main, *_left, *_right, *_g;
  SubRcd* _subs;  // for each curve point, which sub point and how much
  // currently, _g not negated, remember this when passing to SteihaugSolver
  double _D[36], _G1[6], _G2[6];
  
};

#endif

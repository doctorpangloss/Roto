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
 * pyramid.h
 *********************************************************************/

#ifndef _PYRAMID_H_
#define _PYRAMID_H_


#include "klt_util.h"
#include "../jl_vectors.h"

class KLT_Pyramid  {

 public:

  //KLT_Pyramid();

  //void setParams(int basecols, int baserows, int subsampling, int nlevel);

  KLT_Pyramid(int basecols, int baserows, int subsampling, int nlevel, bool initMe = true);

  void computePyramid(KLT_FloatImage* floatimg, float sigma_fact);

  void dumpPyramid();

  void initMe();

  bool inited() const { return _inited; }

  ~KLT_Pyramid();

  KLT_Pyramid(FILE* fp);
  void write(FILE* fp);
  void writeImages(char* name) const;
  void writeDerivImages(char* name) const;
  
  KLT_FloatImage* getFImage(const int r) {
    assert(r>=0 && r<nLevels);
    return img+r;
  }
  const KLT_FloatImage* getFImage(const int r) const {
    assert(r>=0 && r<nLevels);
    return img+r;
  }
  int getNCols(const int r) const {
    return ncols[r];
  }

  int getNRows(const int r) const {
    return nrows[r];
  }

  int getNLevels() const { return nLevels; }


 private:
  int subsampling;
  int nLevels;
  KLT_FloatImage *img;
  int *ncols, *nrows;
  bool _inited;
  int _basencols, _basenrows;

 private:

  void setupStructure(int basencols, int basenrows);
};

class KLT_ColorPyramid {
 public:

  KLT_ColorPyramid(int basecols, int baserows, int subsampling, int nlevel);
  KLT_ColorPyramid(FILE* fp);
  ~KLT_ColorPyramid();

  void smoothAndComputePyramid(const QImage im, const Kernels* kern, float sigma_fact);

  int color(const float x, const float y, const int level, Vec3f& putHere) const;
  // returns 0 if ok, -1 if out of bounds

  KLT_Pyramid* r() { return _r;}
  KLT_Pyramid* g() { return _g;}
  KLT_Pyramid* b() { return _b;}

  int getNLevels() const;

  void write(FILE* fp) const;
  void writeImages(char* name) const;
  void writeDerivImages(char* name) const;

 private:
  KLT_Pyramid *_r, *_g, *_b;
};



#endif

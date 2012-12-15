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



#ifndef GENKEEPER_H
#define GENKEEPER_H

#include "linearSolver.h"

class GenKeeper : public implicitMatrix {
 public:
  virtual ~GenKeeper() {};
  virtual void matVecMult(const double x[], double r[]) const  = 0;
  virtual int numVar() const = 0;
  virtual const double* g() const = 0;  
  virtual void handleConstraints(double*) {}
  virtual bool precond() const {return 0; }
  virtual void precondVec(double*) {}
  virtual void outputMat(char*) {}
  virtual void outputg(char*) const {}
};

#endif

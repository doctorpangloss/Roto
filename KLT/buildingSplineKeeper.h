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



#ifndef BUILDINGSPLINEKEEPER_H
#define BUILDINGSPLINEKEEPER_H

#include <stdio.h>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
typedef boost::numeric::ublas::compressed_vector<double> Ubcv;

//#include "splineKeeper.h"

class BuildingSplineKeeper {


 public:

  BuildingSplineKeeper(int nv);
  BuildingSplineKeeper(BuildingSplineKeeper* other);
  
  ~BuildingSplineKeeper();

  void refresh();

    // J must be of length numVar, J'J added to A and J'r added to b
  void takeJVector(const Ubcv& J, double r);
  
  //void matVecMult(const double x[], double r[]) const;

  void scalarMult(const double m);
  void addMat(const BuildingSplineKeeper* sp);

  void outputJ(char* name);

 private:

  std::vector< std::vector<Ubcv> > _Js;
  std::vector<double> _constants;

  int _nv;
};

#endif

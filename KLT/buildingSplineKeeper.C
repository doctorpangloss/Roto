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



#include "buildingSplineKeeper.h"
#include "myassert.h"

BuildingSplineKeeper::BuildingSplineKeeper(int nv)  : _nv(nv) {
  _constants.push_back(0);
  std::vector < Ubcv > crap;
  _Js.push_back(crap);
}

BuildingSplineKeeper::BuildingSplineKeeper(BuildingSplineKeeper* other) {
  _constants.push_back(0);
  std::vector < Ubcv > crap;
  _Js.push_back(crap);
  _nv = other->_nv;
}
  
BuildingSplineKeeper::~BuildingSplineKeeper() {

}

void BuildingSplineKeeper::refresh() {
  _Js.clear();
  _constants.clear();
  _constants.push_back(0);
  std::vector < Ubcv > crap;
  _Js.push_back(crap);
  //*_g =Zero_Ubv(_nv);
  //_crapgDirty = true;
}


// J must be of length numVar, J'J added to A and J'r added to b
void BuildingSplineKeeper::takeJVector(const Ubcv& J, double r) {
  _Js[0].push_back(J);
  //_crapgDirty = true;
  //Ubcv::const_iterator i1;
  //Ubv& g = *_g;
  //for (i1=J.begin(); i1!=J.end(); ++i1)
  //g(i1.index()) += *i1 * r;
}

/*  
void BuildingSplineKeeper::matVecMult(const double x[], double r[]) const {
  std::vector< std::vector<Ubcv>  >::const_iterator i1;
  std::vector<Ubcv>::const_iterator i2;
  int i;

  Ubv xvec(_nv); 
  Ubv rvec(_nv);
  for (i=0; i<_nv; ++i) 
    xvec(i) = x[i];

  double c;
  i=0;
  for (i1=_Js.begin(); i1 != _Js.end(); ++i1, ++i) {  // iterate over groups
    c = _constants[i];
    for (i2 = i1->begin(); i2 != i1->end(); ++i2) {
      //double con = c * inner_prod(xvec, *i2);
      //Ubcv u = con * *i2; 
      //rvec += u;
	rvec += *i2 * (  c * inner_prod(xvec, *i2) );
    }
  }

  for (i=0; i<_nv; ++i) 
    r[i] = rvec(i);
  
}
*/

void BuildingSplineKeeper::scalarMult(const double m) {
  _constants[0] = m;
  //(*_g) *= m;
  //_crapgDirty = true;
}


void BuildingSplineKeeper::addMat(const BuildingSplineKeeper* sp) {
  assert(sp->_Js.size() == 1);
  _Js.push_back(sp->_Js[0]);
  _constants.push_back(sp->_constants[0]);
  //(*_g) += *(sp->_g);
}


void BuildingSplineKeeper::outputJ(char* name) {
  FILE* fp =fopen(name,"w");
  FILE* fpc = fopen("c.txt","w");

  std::vector< std::vector<Ubcv>  >::const_iterator i1;
  std::vector<Ubcv>::const_iterator i2;
  Ubcv::const_iterator i3;
  double c;
  int i=0, row=0;
  for (i1=_Js.begin(); i1 != _Js.end(); ++i1, ++i) {  // iterate over groups
    c = _constants[i];
    printf("%d observations of constant %.10e\n",i1->size(), c);
    for (i2 = i1->begin(); i2 != i1->end(); ++i2, ++row) { // iterate over vectors
      for (i3 = i2->begin(); i3 != i2->end(); ++i3)  // one vector, list it
	fprintf(fp, "%d %d %.10e\n",row+1, i3.index()+1, *i3);
      fprintf(fpc,"%.10e\n",c);
    }
  }
  fclose(fp);
  fclose(fpc);
}

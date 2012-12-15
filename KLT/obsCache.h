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



#ifndef OBSCACHE_H
#define OBSCACHE_H

#include "splineKeeper.h"
#include "samples.h"

class ObsCache {
  
 public:

  ObsCache(SplineKeeper* sk);

  void newObs();
  
  //void obsFinished();

  bool checkVars(const int v[4], const int i);

  void takeJac(const double J[8], const int ind);

  void takeJac(const TrackSample& ds, const int ind, const double xcoeff, const double ycoeff);

  void outer_prod(const double r);
  void outer_prod(const double r, const double c);
  
  void writeBack();

  double hitRate() const { return double(_hits) / double(_totalObs); }
  void recordMiss(int i) { _totalObs+=i;}

 private:

  int _vars0[4], _vars1[4];
  bool _ok;
  int _totalObs, _hits, _obsWaiting;
  SplineKeeper* _sk;
  double _mat[256], _g[16], _J[16];
};


#endif

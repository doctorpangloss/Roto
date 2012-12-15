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



#ifndef TRACKGRAPH_H
#define TRACKGRAPH_H

#ifdef __APPLE__
using namespace std;
#endif
#include <vector>
#include <algorithm>
#include <stdio.h>
#include "rotoPath.h"
#include "KLT/multiTrackData.h"
#include "KLT/multiSplineData.h"

class RotoPath;
typedef vector<RotoPath*> PathV;


struct JointRecord {
  
  JointRecord() {
    which[0] = which[1] = side[0] = side[1] = -1; 
    //fixed[0] = fixed[1] = false;
    fixer[0] = fixer[1] = NULL; 
  }

  // order is b, e
  int which[2];   // which curve, -1 for resident
  int side[2]; // side, -1 for resident
  //bool fixed[2];  // whether location is fixed, default false
  const RotoPath* fixer[2]; // constraining rotoPath, NULL by default
  int fixerSide[2]; // if fixer, side of fixer
};

typedef vector<JointRecord> JointRV;

class TrackGraph {

 public:

  TrackGraph() {}

  PathV* paths() { return &_paths; }

  bool pathAlreadyThere(const RotoPath* rp);

  void addPath(RotoPath* rp);
  
  void processConstraint(const RotoPath* constrained, const RotoPath* constrainee, 
		       const int edside, const int eeside);

  void traceVarLoc(PathV::const_iterator pc, JointRecord* jr, int pcside, const int side);

  void processJointRecords(const RotoPath* rp);

  void printStructure();

  //void startMulti(MultiTrackData* mt);
  
  //void finishStartingMulti(MultiTrackData* mt);

  void buildMulti(MultiSplineData* mts);
  
  void clear();

  void buildBackReconcileJoints(int bFrame, int aFrame);

  //void copyLocs(MultiTrackData* mt);
  void copyLocs(MultiSplineData* mt);

  const PathV& getKey0Paths() const { return _key0paths; }

 private:

  void getKey0Paths(PathV* putHere, const int numFrames); // can only get called once, not for public-consumption
  
  PathV _paths, _key0paths;
  JointRV _joints;

};


#endif

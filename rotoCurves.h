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



#ifndef ROTOCURVES_H
#define ROTOCURVES_H

#include <fstream>
#include <qdatastream.h>
//#include "llist.h"
#include <list>
#include <deque>
#include "rotoPath.h"
#include "rotoRegion.h"

class RotoRegion;
class RotoPath;
class Joint;
typedef std::vector<Joint> Joints;
typedef std::list<RotoPath*> RotoPathList;
typedef std::vector<RotoPath*> PathV;

class RotoCurves {

  // One image of roto curves

 public:

  RotoCurves();
  void setFrame(int f) { _frame = f; }

  int getNumCurves() const { return _paths.size(); }

  /*
  void startCurvesIterator() const {
    _paths.SetIterationHead();
  }

  RotoPath* getNextCurve();

  const RotoPath* getNextCurve() const;*/

  RotoPathList::const_iterator begin() const { return _paths.begin(); }
  RotoPathList::const_iterator end() const { return _paths.end(); }
  

  RotoPath* getCurve(const int i);

  void addPath(RotoPath* newPath);

  RotoPath* getTail() { return *(_paths.end()); }
  void deleteTail();

  void save(FILE* fp, int version) const;
  void saveqt(QDataStream* qd, int version);
  void load(FILE* fp, int version);
  void loadqt(QDataStream *qd, int version);
  void documentSelf();

  void saveMatte(int w, int h) const;

  void saveRotoXML(std::ofstream& fp, int frame) const;
  void clearXMLLabels();
  void createXMLLabels(int& labelNum);

  void setupJointsOneSide(RotoPath* rp, int side);
  // careful, these procedures may modify number of points to eliminate duplicates at end
  void setupJoints(RotoPath* rp);   
  void setupBegJoints(RotoPath* rp);
  void setupEndJoints(RotoPath* rp);
  void makeJoint(RotoPath* rp1, RotoPath* rp2, int w1, int w2);
  void checkAllJoints();
  //void reconcileJoints();

  //void setActive(RotoPath* a) {
  //_active = a;
  //}
  //RotoPath* getActive() {
  //return _active;
  //}  

  RotoPath* cycle(RotoPath* p);

  void deletePath(RotoPath* dp);

  void addNPaths(const int i);

  RotoPath* pickPrimitive(const int x, const int y);
  RotoPath* pickFreePrimitive(const int x, const int y);

  RotoPath* pickSplineSegment(const Vec2f loc, float& t);

  void renderAllPaths(bool showTrackPoints, bool visEffort=false);
  void renderFreePaths(bool showTrackPoints);
  void renderNextCorrLines();
  void renderRegions();

  void assertNoDraws();

  void printJoints();

  void splitCurve(RotoPath* rp, const int ctrl);
  
  void calcLerp();

  RotoPath* distanceToEnds2(const Vec2f& loc, int& whichEnd);
  RotoPath* distanceToCtrls2(const Vec2f& loc, int* whichControl);

  Vec2i calcStat();

  void addRotoRegion(RotoRegion* r);
  RotoRegion* pickRegion(const float x, const float y);
  int numRegions() const { return _regions.size(); }

  // makes mask by finding lowest region in deck that these paths are
  // part of, and then AND'ing together the masks of regions above
  // it which would occlude them.
  unsigned char* makeCummMask(const PathV& paths);


 private:


  //void reconcileJoint(Joints* js, Vec2f* rploc);

  RotoPathList _paths;
  //RotoPath* _active;
  int _frame;
  //ImgSequence* _is;
  std::deque<RotoRegion*> _regions; // from bottom to top
};


#endif

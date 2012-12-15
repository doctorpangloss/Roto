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



#ifndef DRAWCURVES_H
#define DRAWCURVES_H

#include <qdatastream.h>
#include "llist.h"
#include "drawPath.h"
//#include "drawPatch.h"

class DrawPath;
class DrawPatch;

class DrawCurves {

 public:

  DrawCurves();
  void setFrame(int f) { _frame = f; }

  void addPath(DrawPath* dp);
  //void addPatch(DrawPatch* dp) {
  // _patches.AddToTail(dp);
  //}

  int getNumCurves() const { return _paths.getSize(); }
  //int getNumPatches() const { return _patches.getSize(); }

  void renderAllPaths(const int renderMode);
  //void renderAllPatches();
  
  void renderAllCorr();

  void startDrawPathIterator() {
    _paths.SetIterationHead();
  }
  //void startDrawPatchIterator() {
  //_patches.SetIterationHead();
  //}

  void deletePath(DrawPath* dp);
  void deleteAll();

  //void deletePatch(DrawPatch *dch);

  //DrawPath* pickCurve(const int x, const int y);
  //DrawPatch* pickPatch(const int x, const int y);
  void* pickPrimitive(const int x, const int y, int* which);

  DrawPath* IterateNext();
  //DrawPatch* IterateNextPatch() {
  //return _patches.IterateNext(); }

  DrawPath* getCurve(const int i);
  //DrawPatch* getPatch(const int i);

  void setJustDrawn(const int i);

  void save(FILE* fp, int version) const;
  void saveqt(QDataStream* qd, int version) const;
  void load(FILE* fp, int version);
  void loadqt(QDataStream *qd, int version);
  void documentSelf();

  void addNPaths(const int i);
  //void addNPatches(const int i);

  void MoveToTop(DrawPath* dp);
  void MoveToBottom(DrawPath* dp);

  void MoveOneUp(DrawPath* dp);
  void MoveOneDown(DrawPath* dp);

  DrawPath* cycle(DrawPath* p);

  Vec2i calcStat();

  //void setActivePatch(DrawPatch* a) { _activePatch = a; }
  //DrawPatch *getActivePatch() { return _activePatch; }

 private:

  LList<DrawPath> _paths;
  //LList<DrawPatch> _patches;

  int _frame;
  //DrawPatch *_activePatch;
};



#endif

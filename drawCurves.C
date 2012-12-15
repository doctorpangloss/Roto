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



#include <OpenGl/gl.h>
#include <float.h>
#include "drawCurves.h"

DrawCurves::DrawCurves() {
  //_activePatch = NULL;
}

void DrawCurves::renderAllPaths(const int renderMode) {
  DrawPath* curr;
  _paths.SetIterationHead();
  int name = 0;
  while ((curr=_paths.IterateNext()) != NULL) {
    glPushName(name++);
    curr->render(renderMode);
    glPopName();
  }
}
/*
void DrawCurves::renderAllPatches() {
  DrawPatch* curr;
  _patches.SetIterationHead();
  int name = 1000;
  while ((curr=_patches.IterateNext()) != NULL) { 
    glPushName(name++);
    curr->render();
    glPopName();
  }
}
*/
void DrawCurves::renderAllCorr() {
  DrawPath* curr;
  _paths.SetIterationHead();
  while ((curr=_paths.IterateNext()) != NULL) 
    curr->renderCorr();
}

void DrawCurves::deletePath(DrawPath* dp) {
  if (dp->prevC())
    dp->prevC()->setNextC(NULL);
  if (dp->nextC())
    dp->nextC()->setPrevC(NULL);
  //if (dp->getCorrespondRoto())
  //dp->getCorrespondRoto()->notifyDrawDelete(dp);
  dp->vacateRotoCorr();

  _paths.RemoveNode(dp,1);  
}

/*
void DrawCurves::deletePatch(DrawPatch *dch) {
  if (dch->prevP())
    dch->prevP()->setNextP(NULL);
  if (dch->nextP())
    dch->nextP()->setNextP(NULL);
  _patches.RemoveNode(dch,1);
}
*/
void DrawCurves::deleteAll() {
  startDrawPathIterator();
  DrawPath *curr;
  while ((curr=IterateNext()) != NULL) {
    deletePath(curr);
  }
  assert(_paths.getSize() == 0);
}
/*
DrawPath* DrawCurves::pickCurve(const int x, const int y) {
  _paths.SetIterationHead();
  DrawPath *curr, *minPath;
  double minDist=DBL_MAX, d;
  Vec2f loc(float(x)+.5, float(y)+.5);
  while ((curr=_paths.IterateNext()) != NULL) {
    d = curr->distanceTo2(loc);
    if (d < minDist) {
      minDist = d;
      minPath = curr;
    }
  }
  if (minDist<225)
    return minPath;
  else
    return NULL;
}

DrawPatch* DrawCurves::pickPatch(const int x, const int y) {
  _patches.SetIterationTail();  // better for respecting user expectation of layering
  DrawPatch* curr;
  while ((curr=_patches.IteratePrev()) != NULL) {
    if (curr->onPatch(x,y)) {
      printf("selected patch\n");
      return curr;
    }
  }
  return NULL;
}
*/

void* DrawCurves::pickPrimitive(const int x, const int y, int* which) { 
  printf("picking\n");
  GLint viewport[4];
  GLdouble projection[16];
  GLuint selectBuf[100];
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glSelectBuffer(100, selectBuf);
  glRenderMode(GL_SELECT); 
  glInitNames();
  
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPickMatrix(x, AbstractPath::_globalh-y, 5., 5., viewport);
  glMultMatrixd(projection);
  glMatrixMode(GL_MODELVIEW);

  renderAllPaths(0);
  //renderAllPatches();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glFlush();

  int hits = glRenderMode(GL_RENDER);
  printf("%d hits\n",hits);
  if (hits==0) return NULL;
  uint maxHit = selectBuf[3];
  for (int i=1; i<hits; i++) {
    if (selectBuf[3+4*i] > maxHit) {
      maxHit = selectBuf[3+4*i];
    }
  }
  
  if (maxHit >= 1000) { // patches
    //_patches.SetIterationHead();
    //*which = 0;
    //printf("selected patch\n");
      //return _patches.IterateForward(maxHit - 1000);
  }
  else {
    _paths.SetIterationHead();
    *which = 1;
    printf("selected path\n");
    return _paths.IterateForward(maxHit);
  }
  return NULL;
}


void DrawCurves::setJustDrawn(const int i) {
  DrawPath* curr;
  _paths.SetIterationHead();
  while ((curr=_paths.IterateNext()) != NULL) {
    curr->setJustDrawn(i);
  }
  //_patches.SetIterationHead();
  //DrawPatch *currP;
  //while ((currP=_patches.IterateNext()) != NULL) {
  //currP->setJustDrawn(i);
  //}
}


void DrawCurves::save(FILE* fp, int version) const { 
  _paths.SetIterationHead();
  const DrawPath* curr;
  while ((curr = _paths.IterateNext()) != NULL) 
    curr->save(fp, version);  

  /*
  _patches.SetIterationHead();
  const DrawPatch* currP;
  while ((currP = _patches.IterateNext()) != NULL) 
    currP->save(fp, version);  
  */
}

void DrawCurves::saveqt(QDataStream* fp, int version) const { 
  _paths.SetIterationHead();
  const DrawPath* curr;
  while ((curr = _paths.IterateNext()) != NULL) 
    curr->saveqt(fp, version);  
}

void DrawCurves::load(FILE* fp, int version) { 
    // get data in
  DrawPath* curr;
  _paths.SetIterationHead();
  while ((curr = _paths.IterateNext()) != NULL) {
    curr->load(fp, version);
  }

  /*
  DrawPatch* currP;
  _patches.SetIterationHead();
  while ((currP = _patches.IterateNext()) != NULL) {
    currP->load(fp, version);
  }
  */
}

void DrawCurves::loadqt(QDataStream* fp, int version) { 
    // get data in
  DrawPath* curr;
  _paths.SetIterationHead();
  while ((curr = _paths.IterateNext()) != NULL) {
    curr->loadqt(fp, version);
  }

}

void DrawCurves::documentSelf() { 
  int i=0;
  _paths.SetIterationHead();
  DrawPath* curr;
  while ((curr=_paths.IterateNext()) != NULL) {
    curr->setPtrDoc(_frame,i);
    i++;
  }

  /*
  i=0;
  _patches.SetIterationHead();
  DrawPatch *currC;
  while ((currC=_patches.IterateNext()) != NULL) {
    currC->setPtrDoc(_frame, i);
    i++;
  }
  */
}

void DrawCurves::addNPaths(const int n) {
  for (int i=0; i<n; i++)
    _paths.AddToTail(new DrawPath());
}

//void DrawCurves::addNPatches(const int n) {
//  for (int i=0; i<n; i++)
//    _patches.AddToTail(new DrawPatch());
//}


DrawPath* DrawCurves::getCurve(const int i) {
  _paths.SetIterationHead();
  DrawPath* res = _paths.IterateForward(i);
  assert(res);
  return res;
}

//DrawPatch* DrawCurves::getPatch(const int i) {
//  _patches.SetIterationHead();
//  DrawPatch* res = _patches.IterateForward(i);
//  assert(res);
// return res;
//}


Vec2i DrawCurves::calcStat() {
  Vec2i res;
  DrawPath* curr;
  _paths.SetIterationHead();
  while ((curr=_paths.IterateNext()) != NULL) {
    res.Inc(0,1);
    if (curr->fixed()) 
      res.Inc(1,0);
  }
  return res;
}



DrawPath* DrawCurves::cycle(DrawPath* p) {
  if (p->next) 
    return p->next;
  else
    return _paths.HeadPtr();
}


void DrawCurves::MoveToTop(DrawPath* dp) { _paths.MoveToTail(dp); }
void DrawCurves::MoveToBottom(DrawPath* dp) { _paths.MoveToHead(dp); }

void DrawCurves::MoveOneUp(DrawPath* dp) { _paths.MoveOneUp(dp); }
void DrawCurves::MoveOneDown(DrawPath* dp) { _paths.MoveOneDown(dp); }

DrawPath* DrawCurves::IterateNext() { 
  return _paths.IterateNext(); }

void DrawCurves::addPath(DrawPath* dp) {
  _paths.AddToTail(dp);
}

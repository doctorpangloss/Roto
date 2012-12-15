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



#include <float.h>
#include <qcolordialog.h>
#include <qcursor.h>
#include <qinputdialog.h>
#include "draw.h"

//#define _DEF_DIV_ 10.


bool Draw::QtpickColor(Vec3f& col) {
  QColor qcol = QColorDialog::getColor(getCurrColor());
  if (!qcol.isValid()) return false;
  col.Set(float(qcol.red())/255., 
	  float(qcol.green())/255., 
	  float(qcol.blue())/255.);
  return true;
}

QRgb Draw::getCurrColor() {
  QRgb initC = qRgb(int(_currColor.r()*255.f), int(_currColor.g()*255.f),int( _currColor.b()*255.f));
  return initC;
}

Draw::Draw(ImgSequence *is, QGLWidget *parent) : 
  _is(is), _parent(parent) {
  _w = is->width();
  _h = is->height();
  _frame = is->getStart();
  _mousePressed = 0;
  _currDPath = NULL;
  _selected = NULL;
  //_selectedPatch = NULL;
  _popupPicked = NULL;
  //_popupPatchPicked = NULL;
  //_propPatch = NULL;
  preRegister = false;
  _tool = D_DRAW; 

  _dc = _is->getRelativeDrawCurves(0);
  _corrShow = 0;
  _currColor.Set(0,0,0);
  //_capturePatch = false;
  //_captureDPatch = NULL;
  //DrawPatch::_globalh = _h;
  _sens = 30.f;
  _pressure = 20.f;
  _canNum = 0;
  _spliceFlag=-1;
  _spliceIndex=-1;
  _dAlpha=1.f;
  _zoom = 1.f;
  _copyPtr=NULL;
  _copyFrame=-1;
  _dragDirty = false;
  _justPasted = false;
  _renderMode=0;
  _freeDraw = true;
  _captureColor = false;

  mutualInit();



}


void Draw::paintGL() {
  //printf("in paintGL %d\n",_time.restart());
  // SPEED: enable renderRecent for hertzStrokes
  //if (!_minimalRenderMode && _currDPath) {
  //_currDPath->renderRecent();
  //return;
  //}
  //glLineWidth(3);

  if (_dAlpha==0) {
    if (_selected)    _selected->render(_renderMode);
    if (_currDPath) 
      _currDPath->render(_renderMode,false);
    return;
  }

  //glColor3f(0,0,0);
  _dc->renderAllPaths(_renderMode);
  //_dc->renderAllPatches();

  // debug SKQ -- ??? making corresponded strokes aqua green 
  glColor3f(1,0,0);
  if (_corrShow && _currDPath) {
    _currDPath->renderCorr();
  }
  
  if (_showAllCorr)
    _dc->renderAllCorr();
  
  if (_selected) {
    glColor4f(1,0,0, _dAlpha);
    _selected->renderSelected();
  }

  /*
  if (_selectedPatch) {
    glColor4f(1,0,0, _dAlpha);
    _selectedPatch->renderSelected();
    }*/

  if (_currDPath) 
    _currDPath->render(_renderMode,false);

  if (_tool==D_OVERDRAW && _selected) { // overdraw, direction vis
    glColor4f(1,0,0, _dAlpha);
    _selected->renderDirection();
  }

  if (_captureColor)
    captureColor();
  //printf("ended paintGL %d\n",_time.restart());
}

void Draw::captureColor() {
  glReadPixels(_captureColorLoc.x(), _captureColorLoc.y(),
	       1,1, GL_RGB, GL_FLOAT, (void*)&_currColor);
  printf("captured %f %f %f\n",_currColor.r(), _currColor.g(), _currColor.b());
  QRgb rgb = qRgb(int(_currColor.r()*255.), int(_currColor.g()*255.), int(_currColor.b()*255.));
  QColorDialog::setCustomColor(15, rgb);
  emit currColorChanged();
}

void Draw::frameChange(int i) {
  if (i!=_frame) {
    if (_corrShow) {
      _corrShow = 0;
      _currDPath = NULL;
    }

    // Interpolate forwards any just drawn curves
    if (_propMode == 0) {
      if (i-_frame == 1) {
	DrawPath* curr;
	_dc->startDrawPathIterator();
	while ((curr=_dc->IterateNext()) != NULL)
	  if (curr->justDrawn() && !curr->nextC())
	    curr->interpolateForwards(i);
	//DrawPatch *currP;
	//_dc->startDrawPatchIterator();
	//while ((currP=_dc->IterateNextPatch()) != NULL) {
	//if (currP->justDrawn() && !currP->nextP())
	//  interpolatePatchForwards(currP);
	//}
	
	if (_selected)
	  _selected->interpolateForwards(i);
	//if (_selectedPatch)
	//interpolatePatchForwards(_selectedPatch);
      }
    }
    else if (_propMode==1) { // keyframe
      //if (_selectedPatch) {
      //_propPatch = _selectedPatch;
      //_propFrame = _frame;
      //printf("have prop patch\n");
      //}
    }
    
    if (_tool==D_OVERDRAW && _spliceFlag==-2 && _selected)
      finishSplicedCurve();
    _spliceFlag=-1;

    if (_tool==D_SELECT && _justPasted)
      finishJustPasted();
    
    _dc->setJustDrawn(0);
    if (_selected!=NULL) {
      _selected = NULL;
      emit selectChanged();
    }
    //_selectedPatch = NULL;

    _frame = i;
    _dc = _is->getDrawCurves(_frame);
  }
}



void Draw::keyPressEvent ( QKeyEvent * e ) {
  if (_selected) {
    if (e->key() == Qt::Key_Backspace) {
      _dc->deletePath(_selected);
      if (_selected) {
	_selected = NULL;
	emit selectChanged();
      }
     
    }
    else if (e->key() == Qt::Key_N) {
      _selected = _dc->cycle(_selected);
      emit selectChanged();
      _parent->updateGL();
    }
    else if (e->key() == Qt::Key_Comma && _selected->getStroke()) {
      _selected->decrementThickness(1);
      _selected->freshenAppearance();
    }
    else if (e->key() == Qt::Key_Period && _selected->getStroke()) {
      _selected->incrementThickness(1);
      _selected->freshenAppearance();
    }
    else
      e->ignore();
    _parent->updateGL();
    /*
    if (e->key() == Qt::Key_Comma && _selected->getStroke()) 
      _selected->decrementThickness(1);
    else if (e->key() == Qt::Key_Period && _selected->getStroke()) 
      _selected->incrementThickness(1);

    else if (e->key() == Qt::Key_A && _selected->prevC()) {
      _selected->prevC()->fair();
      printf("A_prime\n");
      interpolateForwards(_selected->prevC(),_frame);
    }
    else if (e->key() == Qt::Key_P && _selected->prevC()) {
      RotoPath* A_prime = _selected->prevC()->getCorrespondRoto();
      assert(A_prime);
      A_prime->fair();
      printf("A_prime\n");
      interpolateForwards(_selected->prevC(),_frame);
    }
    else if (e->key() == Qt::Key_B && _selected->prevC()) {
      assert(_selected->getCorrespondRoto());
      _selected->getCorrespondRoto()->fair();
      printf("B_prime\n");
      interpolateForwards(_selected->prevC(),_frame);
    }
    else if (e->key() == Qt::Key_Backspace) {
      _dc->deletePath(_selected);
      _selected = NULL;
    }
    else if (e->key() == Qt::Key_S) 
      _selected->fair();

    else
      e->ignore();
    _parent->updateGL();
    */
  }

  else
    e->ignore();
}


void Draw::propagateLook(DrawPath* path) {
  assert(path);
  printf("propagating look\n");
  DrawPath* curr = path;
  while ((curr = curr->nextC()) != NULL) {
    curr->copyLook(path);
    curr->freshenAppearance();
  }
  curr = path;
  while ((curr = curr->prevC()) != NULL) {
    curr->copyLook(path);
    curr->freshenAppearance();
  }
}

/*
void Draw::propagatePathEverywhere(DrawPath* path) {
  RotoPath *theRoto = path->getCorrespondRoto(), *currR;
  if (!theRoto) return;
  DrawPath* currD;
  assert(theRoto);
  int frame = _frame;
  
  // forwards
  currR = theRoto->nextC();
  currD = path;
  if (path->nextC()==NULL) {
    while (currR) {
      frame++;
      DrawPath* newPath = new DrawPath();
      newPath->copyLook(path);
      _is->getDrawCurves(frame)->addPath(newPath);
      newPath->setPrevC(currD);
      currD->setNextC(newPath);
      //DrawPath::fillForwardInterpolatedCurve(newPath,path,currR,theRoto);
      currD = newPath;
      currR = currR->nextC();
    } // end while
  }
  
  
  // backwards
  
  currR = theRoto->prevC();
  currD = path;
  frame = _frame;
  if (path->prevC()==NULL) {
    while (currR) {
      frame--;
      DrawPath* newPath = new DrawPath();
      newPath->copyLook(path);
      _is->getDrawCurves(frame)->addPath(newPath);
      newPath->setNextC(currD);
      currD->setPrevC(newPath);
      DrawPath::fillBackwardInterpolatedCurve(newPath,path,currR,theRoto);
      currD = newPath;
      currR = currR->prevC();
    } // end while
  }
  B
  
  
}
*/

/*
notes:
1. aPath and bPath are drawPaths from the first and last frame over which you are optimizing
the shape.  bPath may have a different number of points than aPath, but it is assumed that 
all the inbetween DrawPaths have the same number of points as aPath.
2. aPath should be a keyframe, and won't be modified.  The argument pinLast controls whether bPath
should be fixed (pinned), or whether it's shape should also be modified.  If you wish to modify bPath's
shape, it should have the same number of points as aPath.  It not, it will not be changed.
3. aFrame and bFrame are the frame numbers of aPath and bPath.
4. aPath and bPath better form a chain of DrawPath's linked up, so that bPath is eventually reached
by following  nextC() pointers in (bFrame-aFrame) steps.
*/
void Draw::optimizeDrawShape(DrawPath* aPath, DrawPath* bPath, int aFrame, int bFrame, bool pinLast) {
  assert(bFrame > aFrame);
  int numP = aPath->getNumElements(), numFrames = bFrame-aFrame, numVec, t, nc;

  if (pinLast) {
    if (numFrames<2) return;
    numVec = numP*(numFrames-1);
    if (bPath->getNumElements()!=aPath->getNumElements())
      nc = numFrames-1;
    else
      nc=numFrames;
  }
  else {
    numVec=numP*numFrames;
    nc = numFrames;
  }

  Vec2f *Plocs = new Vec2f[numP];
  Vec2f *Z = new Vec2f[numP * numFrames], *tZ;
  Vec2f *Sb = new Vec2f[numP * numVec], *tSb;
  memcpy(Plocs,aPath->getData(), numP*sizeof(double));
  DrawPath* curr = aPath;

  for (t=0, tZ=Z, tSb=Sb; t<nc; ++t, tZ+=numP, tSb+=numP) {
    curr = curr->nextC();
    assert(curr);
    assert(curr->getNumElements() == numP);
    assert(curr->getShouldBe() && curr->getShouldBe()->getNumElements()==numP);
    memcpy(tZ, curr->getData(), numP*sizeof(Vec2f)); 
    memcpy(tSb, curr->getShouldBe()->getData(), numP*sizeof(Vec2f));
  }
  
  if (nc==numFrames-1) { // last frame might have different number of points
    curr = curr->nextC();
    assert(curr && curr==bPath);
    int theirnp;
    int* t1d1 = curr->getT1D1(&theirnp);
    assert(t1d1 && theirnp==numP);
    for (int j=0; j<numP; ++j, ++tZ) 
      (*tZ) = bPath->getElement(t1d1[j]);
      
  }
  
  _is->globalTC.longCurveReshape(numFrames, numP, Plocs, Z, Sb, pinLast);
  
  

  curr = aPath;
  for (t=0, tZ=Z; t<nc; ++t, tZ+=numP) {
    curr = curr->nextC();
    memcpy(curr->getData(), tZ, numP*sizeof(Vec2f));
    curr->redoStroke();
    curr->freshenAppearance();
  }

  delete[] Plocs; delete[] Z; delete[] Sb;
}

void Draw::rebuildStrokeSpan(DrawPath* path) {
  if (path->fixed()) {
    printf("Cannot rebuild stroke span from a fixed curve\n");
    return;
  }

  // find stopping points
  DrawPath *aPath, *bPath, *curr = path;
  int aFrame = _frame-1, bFrame = _frame+1;
  
  curr = path->prevC();
  while (curr && !curr->fixed() && curr->prevC() != NULL) {
    curr = curr->prevC();
    aFrame--;
  }
  aPath = curr;
  
  curr = path->nextC(); 
  while (curr && !curr->fixed() && curr->nextC() != NULL) {
    curr = curr->nextC();
    bFrame++;
  }
  bPath = curr;
  if (bPath==NULL) {
    bPath = path; --bFrame;
  }

  //optimizeDrawShape(aPath, bPath, aFrame, bFrame, bPath->fixed()); G!
}


void Draw::propagatePathEverywhere(DrawPath* path) {
  _parent->setCursor(Qt::WaitCursor);
  //RotoPath **currRotos = new RotoPath*[path->getNumElements()];
  bool kosher = true;
  
  DrawPath* currD;
  int frame = _frame, shit, i;
   
  // forwards
  kosher = path->getCorrRoto()->nextRotosExist();
  //for (i=0; i<path->getNumElements(); i++) {
  //currRotos[i] = path->getRotoCorr(i,shit)->nextC();
  //if (!currRotos[i]) kosher = false;
  //}
  currD = path;
  if (path->nextC()==NULL) {
    while (kosher) {
      frame++;      
      DrawPath* newPath = new DrawPath();
      newPath->setFixed(false);
      newPath->copyLook(path);
      _is->getDrawCurves(frame)->addPath(newPath);
      newPath->setPrevC(currD);
      currD->setNextC(newPath);
      DrawPath::fillForwardInterpolatedCurve(newPath,currD);
      currD = newPath;      
      kosher = currD->getCorrRoto()->nextRotosExist();

      currD->redoStroke(); // HM?
      currD->freshenAppearance();

      //for (i=0; i<path->getNumElements(); i++) {
      //currRotos[i] = currRotos[i]->nextC();
      //if (!currRotos[i]) kosher = false;
      //}
    } // end while
  }
  if (frame > _frame && 0) //G!
    optimizeDrawShape(path, currD, _frame, frame, false);

  /*
  // backwards
  kosher = true;
  for (i=0; i<path->getNumElements(); i++) {
    currRotos[i] = path->getRotoCorr(i,shit)->prevC();
    if (!currRotos[i]) kosher = false;
  }  
  currD = path;
  frame = _frame;
  if (path->prevC()==NULL) {
    while (kosher) {
      frame--;
      DrawPath* newPath = new DrawPath();
      newPath->setFixed(false);
      newPath->copyLook(path);
      _is->getDrawCurves(frame)->addPath(newPath);
      newPath->setNextC(currD);
      currD->setPrevC(newPath);
      DrawPath::fillBackwardInterpolatedCurve(newPath,currD,currRotos);
      currD = newPath;
      for (i=0; i<path->getNumElements(); i++) {
	currRotos[i] = currRotos[i]->prevC();
	if (!currRotos[i]) kosher = false;
      }
    } // end while
  }
  */
  emit restoreCursor();
}

void Draw::makeDrawPath(const RotoPath* rp) {
  DrawPath* newp = new DrawPath();
  newp->setColor(_currColor);
  newp->initHertzStroke();
  for (int i=0; i<rp->getNumElements(); ++i)
    newp->addVertex(rp->getElement(i).x(), rp->getElement(i).y(), mapPressure());
  _dc->addPath(newp);
  newp->freshenAppearance();
}


void Draw::mousePressEvent( QMouseEvent *e ) {
  //  printf("mouse press event %d\n", e->button());
  // draw
  _mousePressed = 1;
  if (e->button() == Qt::LeftButton && (_tool == D_DRAW || _tool == D_PAINT)) {
    if (_corrShow) {
      _corrShow = 0;
      _currDPath = NULL;
    }
    if (_currDPath != NULL)
      return;
    _pressure = 20.f;
    //printf("mouse press pressure %f %f\n",_pressure, _sens);
    
    _currDPath = new DrawPath();
    _currDPath->setColor(_currColor);


    _currDPath->initHertzStroke();        // comment out to avoid Hertzman strokes
    Vec2f loc = unproject(e->x(), e->y(), _h);
    //_currDPath->addVertex(loc.x(), loc.y(), mapPressure());
    //_minimalRenderMode = 1;
  }

  else if (_tool==D_SELECT) { // select
    int which; // 0 for patch, 1 for path, -1 for neither
    void* selected = _dc->pickPrimitive(e->x(), e->y(), &which);
    _dragLoc = unproject(e->x(), e->y(), _h);
    if (selected) {
      if (which==0) {
	//_selectedPatch = (DrawPatch*) selected;
	_selected = NULL;
      }
      else {
	_selected = (DrawPath*) selected;
	//_selectedPatch = NULL;
      }
    }
    else {
      if (_justPasted)
	finishJustPasted();
      //_selectedPatch = NULL;
      _selected = NULL;
    }
    emit selectChanged();
  }

  else if (e->button() == Qt::LeftButton && _tool==D_OVERDRAW) { // overdraw
    Vec2f loc = unproject(e->x(), e->y(), _h);
    if (_spliceFlag > -1) { // already drawing
    }
    else {
      if (!_selected) {
	printf("Please select a curve to edit, first.\n");
	return;
      }
      
      double d1 = _selected->distanceTo2(loc, _spliceIndex);
      if (d1<16) { // starting on curve
	_spliceFlag=1;
      }
      else {      // new beginning to curve
	_spliceFlag=0;
	_spliceIndex=0;
      }

      assert(!_currDPath);
      _currDPath = new DrawPath();
      _currDPath->setColor(_selected->getColor());
      _currDPath->initHertzStroke();        
      //_currDPath->addVertex(loc.x(), loc.y(), mapPressure());
    }
  }

  else if (e->button() == Qt::LeftButton && _tool==D_NUDGE) { // nudge
    if (!_selected) {
      printf("Please select a curve to edit, first.\n");
      _spliceIndex=-1;
      return;
    }
    Vec2f loc = unproject(e->x(), e->y(), _h);
    double d1 = _selected->distanceTo2(loc, _spliceIndex);
    if (d1>9) {
      _spliceIndex=-1;
      return;
    }
    
  }

  _parent->updateGL(); 
}


void Draw::mouseReleaseEvent( QMouseEvent *e ) {
  //  printf("mouse release event %d\n", e->button());
  // draw
  _mousePressed = 0;
  if (_tool == D_DRAW && _currDPath) {  // draw
    //_pressure = 2.*_sens; 
    if (_currDPath->getNumElements() > 1) {
      _parent->setCursor(Qt::WaitCursor);
      Vec2f loc = unproject(e->x(), e->y(), _h);
      _currDPath->addVertex(loc.x(), loc.y(), mapPressure());
      _dc->addPath(_currDPath);

      //for (int i=0; i<2; i++)
      //_currDPath->fair();
      _currDPath->resample(NULL);
      _currDPath->redoStroke();

      bool res=false;
      if (!_freeDraw)
	res = _currDPath->correspondRoto2(_is->getRotoCurves(_frame), _w, _h);
      _currDPath->calculateStrokeDisplayList();
      if (res) {
	_corrShow = 1;
	propagatePathEverywhere(_currDPath);
      }
      else _currDPath = NULL;
      emit restoreCursor();
    }
    else {
      delete _currDPath;
      _currDPath = NULL;
    }
  }

  /*
  else if (_tool == D_PAINT && _currDPath) { // paint
    Bboxf2D bbox = _currDPath->calcBbox();
    DrawPatch *dch = new DrawPatch((int)floor(bbox.lower.x()),(int)floor(bbox.lower.y()),
				   (int)ceil(bbox.width()), (int)ceil(bbox.height()),
				   _currColor);
    dch->takeStroke(_currDPath->giveStroke());
    _currDPath->forgetStroke();
    _dc->addPatch(dch);
    //printf("patch corner: %f %f, window: %f %f\n",dch->x(), dch->y(),
    //   bbox.width(), bbox.height());
    delete _currDPath;
    _currDPath = NULL;

    if (_propPatch && _propMode==1)
      keyframeInterpolatePatch(dch);
  }
  */

  else if (_tool == D_SELECT) { // select
    if (_selected && _dragDirty) {
      regeneratePath(_selected);
      _selected->fixLoc();
      _dragDirty = false;
    }
  }
  else if (_tool == D_OVERDRAW && _currDPath) {// overdraw 
    if (_currDPath->getNumElements() > 1) {
      _currDPath->resample(NULL);
      spliceDrawCurve();
    }
    else {
      delete _currDPath; _currDPath=NULL;
      _spliceFlag=-1; 
    }
  }

  else if (_tool==D_DROPPER) {
    _captureColor = true;
    _captureColorLoc.Set(e->x(), _h-e->y());
    _parent->updateGL();
    _captureColor = false;
  }

  else if (_tool==D_NUDGE) { // nudge
    _spliceIndex=-1;
  }

  //_minimalRenderMode = 0;
  _parent->updateGL();
}


void Draw::mouseMoveEvent( QMouseEvent *e ) {
  // printf("called mouse Move %d %d %d\n",_time.restart(), e->x(), e->y());
  // draw
  if ((_tool==D_DRAW || _tool==D_OVERDRAW || _tool==D_PAINT) && _mousePressed && _currDPath) {
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    //printf("%f %f\n",newLoc.x(), newLoc.y());
    if (_currDPath->getNumElements()==0 || _currDPath->distToLast2(newLoc) > 0)
      _currDPath->addVertex(newLoc.x(),newLoc.y(), mapPressure());
    QPoint curs = _parent->mapFromGlobal(QCursor::pos());
    if ((e->x() - curs.x())*(e->x() - curs.x()) +
	(e->y() - curs.y())*(e->y() - curs.y()) < 25)
      _parent->updateGL();
  }
  else if (_tool == D_SELECT && _mousePressed && _selected) {
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    Vec2f delta(newLoc, _dragLoc);
    _selected->translate(delta);
    _dragLoc = newLoc;
    _dragDirty = true;
    _parent->updateGL();
  }
  else if (_spliceIndex>0 && _selected && _tool==D_NUDGE) { // nudge
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    _selected->nudge(_spliceIndex, newLoc);
  }
}


void Draw::spliceDrawCurve() {
  assert(_currDPath && _spliceFlag > -1 && _selected);
  Vec2f loc = _currDPath->getElement(_currDPath->getNumElements()-1);  
  if (_selected->corrToRoto()) 
    _spliceFlag=-2;
  else
    _spliceFlag=-1;

  int ind2;
  double d2 = _selected->distanceTo2(loc, ind2);
  if (d2<16) {  // ending on curve
    _selected->spliceIn(_currDPath,_spliceIndex,ind2);
  }
  else {      // new ending of curve
    _selected->spliceIn(_currDPath,_spliceIndex,
			_selected->getNumElements()-1);
  }

  _selected->setFixed(true);
  delete _currDPath; _currDPath=NULL;
}

void Draw::regeneratePath(DrawPath* path) {
  _parent->setCursor(Qt::WaitCursor);
  printf("regenerate path\n");
  DrawPath *aPath, *bPath, *curr = path;
  int aFrame = _frame-1, bFrame = _frame+1;
  
  curr = path->prevC();
  while (curr && !curr->fixed() && curr->prevC() != NULL) {
    curr = curr->prevC();
    aFrame--;
  }
  aPath = curr;
  
  curr = path->nextC(); 
  while (curr && !curr->fixed() && curr->nextC() != NULL) {
    curr = curr->nextC();
    bFrame++;
  }
  bPath = curr;

  // 1. regenerate forwards and backwards curves
  if (aPath && aPath!= _selected && aPath->nextC()!=_selected)
    biGenerateSpan(aPath, aFrame, _selected, _frame);
  if (bPath && bPath!=_selected && bPath->prevC()!=_selected) {
    if (bPath->fixed())
      biGenerateSpan(_selected, _frame, bPath, bFrame);
    else

      forwGenerateSpan(_selected, aFrame);
  }
  
  
  emit restoreCursor();
  _selected->setFixed(true);
}


void Draw::finishSplicedCurve() {
  printf("finishSplivedCUrve\n");
  // find stopping points
  DrawPath *aPath, *bPath, *curr = _selected;
  int aFrame = _frame-1, bFrame = _frame+1;
  
  curr = _selected->prevC();
  while (curr && !curr->fixed() && curr->prevC() != NULL) {
    curr = curr->prevC();
    aFrame--;
  }
  aPath = curr;
  
  curr = _selected->nextC(); 
  while (curr && !curr->fixed() && curr->nextC() != NULL) {
    curr = curr->nextC();
    bFrame++;
  }
  bPath = curr;
  
  // 2. Correspond _selected to rotocurves
  _selected->vacateRotoCorr();
  _selected->correspondRoto2(_is->getRotoCurves(_frame), _w, _h);

  // 1. regenerate forwards and backwards curves
  if (aPath && aPath!= _selected && aPath->nextC()!=_selected && aPath->corrToRoto())
    biGenerateSpan(aPath, aFrame, _selected, _frame);
  if (bPath && bPath!=_selected && bPath->prevC()!=_selected && bPath->corrToRoto()) {
    if (bPath->fixed())
      biGenerateSpan(_selected, _frame, bPath, bFrame);
    else
      forwGenerateSpan(_selected, aFrame);
  }

  

  _selected->setFixed(true);
}


//: TO DO:
void Draw::forwGenerateSpan(DrawPath* aPath, int aFrame) {
  /*  DrawPath* curr = aPath->nextC(), *bPath, *prev = aPath;
  RotoPath **currRotos = new RotoPath*[aPath->getNumElements()];
  bool kosher = true;
  int i, shit, frame = aFrame;
  for (i=0; i<aPath->getNumElements(); i++) {
    currRotos[i] = aPath->getRotoCorr(i,shit)->nextC();
    if (!currRotos[i]) kosher = false;
  }

  printf("Forward rebuilding starting at frame %d\n",aFrame);
  while (kosher && curr && !curr->fixed()) {  // questionable
    curr->empty();
    DrawPath::fillForwardInterpolatedCurve(curr,prev, currRotos);
    curr->calculateStrokeDisplayList();
    curr->calculateFillDisplayList();
    for (i=0; i<aPath->getNumElements(); i++) {
      currRotos[i] = currRotos[i]->nextC();
      if (!currRotos[i]) kosher = false;
    }
    bPath = curr;
    prev = curr;
    curr = curr->nextC();
    ++frame;
  }
  printf("Final frame modified is %d\n", frame);
  
  delete[] currRotos;
  optimizeDrawShape(aPath, bPath, aFrame, frame, false);
  */
}


void Draw::biGenerateSpan(DrawPath* aPath, const int aFrame,
			  DrawPath* bPath, const int bFrame) {
  
  /*
  int i, j;

  int numP = aPath->getNumElements();
  int *D0_T0_T1_i = new int[numP];
  RotoPath** D0_T0_T1_p = new RotoPath*[numP];
  int *D0_T0_T_new_i = new int[numP];
  RotoPath** D0_T0_T_new_p = new RotoPath*[numP];

  for (j=0; j<numP; j++) {
    D0_T0_T_new_p[j] = D0_T0_T1_p[j] = aPath->getRotoCorr(j, D0_T0_T1_i[j]);
    D0_T0_T_new_i[j] = D0_T0_T1_i[j];
    for (i=aFrame; i<bFrame; i++) {
      D0_T0_T1_i[j] = D0_T0_T1_p[j]->getNextCorrIndex(D0_T0_T1_i[j]);
      D0_T0_T1_p[j] = D0_T0_T1_p[j]->nextC();
      if (!D0_T0_T1_p[j]) {
	delete[] D0_T0_T1_i; delete[] D0_T0_T1_p;
	delete[] D0_T0_T_new_i; delete[] D0_T0_T_new_p;
	printf("Missing a tracking curve.  Did you delete it?\n");
	return;
      }
    }
  }

  int *T1_D1 = new int[numP];
  AbstractPath T1;
  for (j=0; j<numP; j++)
    T1.add(D0_T0_T1_p[j]->getElement(D0_T0_T1_i[j]));
  RotoCorresponder* myCorr = new RotoCorresponder(&T1,bPath);
  myCorr->calculate();
  myCorr->getSolution(T1_D1);
  delete myCorr;
  bPath->takeT1D1(T1_D1, numP);
  
  DrawPath* currR = aPath;

  for (i=aFrame+1; i<bFrame; i++) {
    for (j=0; j<numP; j++) {
      D0_T0_T_new_i[j] = D0_T0_T_new_p[j]->getNextCorrIndex(D0_T0_T_new_i[j]);
      D0_T0_T_new_p[j] = D0_T0_T_new_p[j]->nextC();
      assert(D0_T0_T_new_p[j]);
    }
    currR = currR->nextC();  assert(currR);
    DrawPath::fillBiInterpolatedCurve(aPath, bPath, currR, double(i-aFrame)/double(bFrame-aFrame), 
				      D0_T0_T1_i, D0_T0_T1_p,
				      D0_T0_T_new_i,D0_T0_T_new_p, T1_D1 );
  }

  delete[] D0_T0_T1_i; delete[] D0_T0_T1_p;
  delete[] D0_T0_T_new_i; delete[] D0_T0_T_new_p;
  //delete[] T1_D1;

  optimizeDrawShape(aPath, bPath, aFrame, bFrame, true);
  */
}

/*
void Draw::longTrackSpan(DrawPatch *dch) {
  assert(dch);
  DrawPatch *aPatch, *bPatch, *curr = dch;
  int aFrame = _frame, bFrame = _frame, numFrames;

  while (curr->prevP() != NULL) {
    curr = curr->prevP();
    aFrame--;
  }
  aPatch = curr;

  curr = dch; 
  while (curr->nextP() != NULL) {
    curr = curr->nextP();
    bFrame++;
  }
  bPatch = curr;
  numFrames = bFrame - aFrame;
  if (numFrames < 1) {
    printf("Too short a span\n");
    return;
  }
  printf("aframe: %d, bframe: %d\n",aFrame, bFrame);
  
  const KLT_FullCPyramid** pyrms = new (const KLT_FullCPyramid*)[numFrames+1];
  double* param = new double[numFrames * 6];
  pyrms[0] = _is->getFullCPyramid(aFrame); 
  int i, offset;

  curr = aPatch;  
  Vec2f loca = aPatch->center(), locb = bPatch->center();

  for (i=0, offset=0; i < numFrames; i++, offset+=6) { // fill param values
    curr = curr->nextP(); assert(curr);
    //if (i==numFrames-1) { // DEBUG
    // param[offset+0] = param[offset+1] = param[offset+2] = param[offset+3] = 0;
    // param[offset+4] = locb.x() - loca.x();
    // param[offset+5] = locb.y() - loca.y();
    //}      
    //else {
    curr->getTrans(param+offset);
    param[offset+0] -=1.; param[offset+3] -=1;  // make deltas
    //}
    pyrms[i+1] = _is->getFullCPyramid(aFrame+i+1); // fill pyrms list 
  }

  printDoubleArray(stdout, param, numFrames, 6);
  int res = _is->globalTC.longPatchTrack(pyrms, numFrames, loca.x(), loca.y(), aPatch->width(), 
					 aPatch->height(), param); 

  if (res==KLT_TRACKED  || res == KLT_MAX_ITERATIONS) {
    curr = aPatch;
    int ul = numFrames-1;
    if (!_is->globalTC.pinLast && bPatch->x() == aPatch->x() && bPatch->y() == aPatch->y()) {
      ul = numFrames;
      printf("moving last frame\n");
    }
    for (i=0, offset=0; i < ul; i++, offset+=6) {
      curr = curr->nextP();
      assert(curr);
      printf("(%.3f,%.3f)",param[offset+4],param[offset+5]);
      //curr->setXY(aPatch->x()+param[offset+4], aPatch->y()+param[offset+5]);
      curr->setTrans(param+offset);
    }
    printf("\n");
  }
  else printf("Failed to track patch\n");
  
  delete[] param;
  delete[] pyrms;
}

void Draw::keyframeInterpolatePatch(DrawPatch *dch) {
  assert(_propMode==1 && _propPatch);
  printf("interpolating keyframe patches\n");
  int numFrames =_frame - _propFrame, res;
  if (numFrames <=0) {
    _propPatch = NULL;
    return;
  }
  

  Vec2f loca, locb;  
  loca = _propPatch->center();
  locb = dch->center();
  Vec2f locDelta(locb, loca); // locb - loca

  // Start by registering first and last frame
  double Z[12];
  if (preRegister) {
    memset(Z,0,12*sizeof(double));    
    Z[10] = locDelta.x(); Z[11] = locDelta.y();
    printf("original %.3f %.3f\n",Z[10], Z[11]);
    res = _is->globalTC.colorTrack(_is->getFullCPyramid(_propFrame), _is->getFullCPyramid(_frame),
			   _propPatch->width(), _propPatch->height(),
			     loca.x(), loca.y(), Z);
  }
  if (preRegister) {
    printf("final %.3f %.3f\n",Z[10], Z[11]);
    if (res == KLT_TRACKED && 
	(fabs(Z[10]-locDelta.x()) < _propPatch->width() &&
	 fabs(Z[11]-locDelta.y()) < _propPatch->height() ) ) {
      dch->setTrans(Z+6);
      dch->gutMe(_propPatch);
    }
    else {
      printf("pre-registration failed\n");
      preRegister = false;
    }
  }
  _parent->updateGL();

  // long Tracking
  int i, offset;
  const KLT_FullCPyramid** pyrms = new (const KLT_FullCPyramid*)[numFrames+1];
  float fnumFrames = float(numFrames);
  double* param = new double[numFrames * 6];
  pyrms[0] = _is->getFullCPyramid(_propFrame); 

  DrawPatch *prev = _propPatch;
  for (i=0, offset=0; i < numFrames; i++, offset+=6) {
    double t = float(i+1)/fnumFrames; assert(t>= 0 && t<=1.);
    param[offset+4] = t*(locb.x() - loca.x()); 
    param[offset+5] = t*(locb.y() - loca.y());
    param[offset+0] = param[offset+1] = param[offset+2] = param[offset+3] = 0;
    pyrms[i+1] = _is->getFullCPyramid(_propFrame+i+1);  
    if (i>0) 
      assert(pyrms[i+1]->nPyramidLevels == pyrms[i]->nPyramidLevels);

    if (i!=numFrames-1) {
      DrawPatch *newdch = new DrawPatch(_propPatch);
      newdch->setTransTranslate(param[offset+4], param[offset+5]);
      _is->getDrawCurves(_propFrame+i+1)->addPatch(newdch);
      newdch->setJustDrawn(0);
      newdch->setPrevP(prev);
      prev->setNextP(newdch);
      prev = newdch;
    }
  }
  if (preRegister) {
    memcpy(param + (numFrames-1)*6, Z+6, 6*sizeof(double));
  }
  dch->setPrevP(prev);
  prev->setNextP(dch);

  printf("set up initial conditions\n");
  res = _is->globalTC.longPatchTrack(pyrms, numFrames, loca.x(), loca.y(), _propPatch->width(), 
					 _propPatch->height(), param);
  
  if (res==KLT_TRACKED || res == KLT_MAX_ITERATIONS) {
    DrawPatch* curr = _propPatch;
    for (i=0, offset=0; i < numFrames-1; i++, offset+=6) {
      curr = curr->nextP();
      assert(curr);
      printf("("); 
      for (int tt=0; tt<6; tt++)
	printf("%.3f,",param[offset+tt]);
      printf(")\n");
      curr->setTrans(param+offset);
    }
    printf("\n");
    dch->setTrans(param+(numFrames-1)*6);  // DEBUG
    dch->gutMe(_propPatch);
  }
  else printf("Failed to track patch\n");
  
  delete[] param;
  delete[] pyrms;
  _propPatch = NULL;
}

void Draw::interpolatePatchForwards(DrawPatch* A) {
  if (A->nextP()) return;

  KLT_FullCPyramid *aPyrm = _is->getFullCPyramid(_frame), 
    *bPyrm = _is->getFullCPyramid(_frame+1);

  // track
  Vec2f center = A->center();
  Vec2f lastT = A->getLastTranslation();
  double Z[12];
  A->getTrans(Z);     Z[0] -= 1.f; Z[3] -= 1.f;
  A->getTrans(Z+6);   Z[6] -= 1.f; Z[9] -= 1.f;
  Z[10] += lastT.x(); Z[11] += lastT.y();
  _is->globalTC.colorTrack(aPyrm, bPyrm, A->width(), A->height(), center.x(), center.y(), Z); 
  printf("translation %f %f, loc %f %f\n\n\n",Z[10], Z[11], center.x()+Z[10], center.y()+Z[11]);

  // apply transform to patch (remember it's a delta)
  DrawPatch *B = new DrawPatch(A);
  B->setLastTranslate(Z[10] - Z[4], Z[11] - Z[5]);
  B->setTrans(Z+6);
  A->setNextP(B);
  B->setPrevP(A);
  _is->getDrawCurves(_frame+1)->addPatch(B);
}

*/
/*
void Draw::visualize(const DrawPath* A, const RotoPath* A_prime,
		     const DrawPath* B, const RotoPath* B_prime) {
  
  FILE* fp = fopen("vis.eps","w");
  double scale = 10;
  fprintf(fp,"%!PS-Adobe-2.0 EPSF 2.0\n");
  fprintf(fp,"%%%%BoundingBox: 0 0 %f %f\n",scale*_w,scale*_h);
  fprintf(fp,"%f setlinewidth\n",0.25);
  fprintf(fp,"%f %f scale\n",scale,scale);



  fprintf(fp, "1 0 0 setrgbcolor\n");
  A->renderCorrPS(fp);



  fprintf(fp, ".75 .75 .75 setrgbcolor\n");
  A_prime->renderNextCorrPS(fp);



  fprintf(fp, "1 1 0 setrgbcolor\n");
  B->renderCorrPS(fp);

  fprintf(fp, "0 0 1 setrgbcolor\n");
  B->renderPS(fp);

  fprintf(fp, "0 0 0 setrgbcolor\n");
  A->renderPS(fp);

  fprintf(fp, "0 1 1 setrgbcolor\n");
  B_prime->renderPS(fp);

  fprintf(fp, "1 0 1 setrgbcolor\n");
  A_prime->renderPS(fp);

  fprintf(fp,"showpage\n");
  fclose(fp);
}
*/



void Draw::toolChange(DrawTool id) {
  if (_tool==D_OVERDRAW && id!=D_OVERDRAW && _spliceFlag==-2 && _selected)
    finishSplicedCurve();
  _spliceFlag=-1;  
  if (_tool==D_SELECT && id!=D_SELECT && _justPasted)
    finishJustPasted();
  
  _tool = id;
  if (_corrShow) {
    _corrShow = 0;
    _currDPath = NULL;
  }

  if (id==D_DRAW) {
    if (_selected) {
    _selected = NULL;
    emit selectChanged();
    }
    //_selectedPatch = NULL;
  }
  _parent->updateGL();
}


CmdMap* Draw::getPopupFunctions(const Vec2i& loc, DynArray<int,5>* enab) {
  CmdMap *cm = new CmdMap();
  DrawPath* dp=NULL; //DrawPatch *dcp=NULL;

  int which;
  void* selected = _dc->pickPrimitive(loc.x(), loc.y(), &which);
  if (selected) {
    //if (which==0) 
    //dcp = (DrawPatch*) selected;
    //else 
      dp = (DrawPath*) selected;
  }
  
  //  if (dcp) {
  //cm->insert(make_pair(LONG_TRACK_SPAN, string("Long track span")));
  //cm->insert(make_pair(DELETE_PATCH, string("Delete patch")));
  //}
  if (dp) {
    cm->insert(make_pair(PROPAGATE_LOOK, string("Propagate look")));
    cm->insert(make_pair(CHANGE_COLOR, string("Change contour color...")));
    cm->insert(make_pair(CHANGE_FILL_COLOR, string("Change fill color...")));
    cm->insert(make_pair(UNFILL, string("Unfill")));
    cm->insert(make_pair(SET_THICKNESS, string("Set to uniform thickness...")));
    if (!dp->filled())
      enab->add(UNFILL);
    cm->insert(make_pair(CORRESPOND_ROTO, string("Correspond to track curve")));
    cm->insert(make_pair(UNCORRESPOND, string("Uncorrespond to track curve")));
    cm->insert(make_pair(PROPAGATE_MORE, string("Propagate curve")));
    if (dp->corrToRoto())
      enab->add(CORRESPOND_ROTO);    
    else {
      enab->add(PROPAGATE_MORE);
      enab->add(UNCORRESPOND);
    }
    cm->insert(make_pair(MOVE_TO_TOP, string("Move to top")));
    cm->insert(make_pair(MOVE_SPAN_TO_TOP, string("Move Span to top")));
    cm->insert(make_pair(MOVE_TO_BOTTOM, string("Move to bottom")));
    cm->insert(make_pair(MOVE_SPAN_TO_BOTTOM, string("Move Span to bottom")));
    cm->insert(make_pair(MOVE_ONE_UP, string("Move one up")));
    cm->insert(make_pair(MOVE_ONE_DOWN, string("Move one down")));
    cm->insert(make_pair(DELETE, string("Delete stroke")));
    cm->insert(make_pair(DELETE_SPAN, string("Delete whole span")));
    cm->insert(make_pair(SELECTD, string("Select")));
    cm->insert(make_pair(FIXD, string("Fix curve location")));
    if (dp->fixed())
      enab->add(FIXD);
    cm->insert(make_pair(REBUILD_STROKE_SPAN, string("Rebuild stroke span")));
    if (!dp->getShouldBe())
      enab->add(REBUILD_STROKE_SPAN);
    cm->insert(make_pair(REDO_FORWARD_SPAN, string("Redo forward span")));
    if (!dp->fixed() || !dp->nextC())
      enab->add(REDO_FORWARD_SPAN);
  }


  cm->insert(make_pair(DELETE_ALL, string("Delete All")));
  //cm->insert(make_pair(ACTIVATE_ALL, string("Activate All")));
  
  _popupPicked = dp;
  //_popupPatchPicked = dcp;
  return cm;
}

void Draw::callPopup(const int which) {
  printf("%d\n",which);
  
  Vec3f col;
  switch(which) {
  case LONG_TRACK_SPAN:
    //longTrackSpan(_popupPatchPicked);
    break;
  case PROPAGATE_MORE:
    propagatePathEverywhere(_popupPicked);
    break;
  case PROPAGATE_LOOK:
    propagateLook(_popupPicked);
    break;
  case CHANGE_COLOR:
    if (QtpickColor(col)) {
      _popupPicked->setColor(col);
      _popupPicked->calculateStrokeDisplayList();
    }
    break;
  case CHANGE_FILL_COLOR:
    if (QtpickColor(col)) 
      _popupPicked->setFillColor(col);
    break;
  case CORRESPOND_ROTO:
    _popupPicked->correspondRoto2(_is->getRotoCurves(_frame), _w, _h);
    break;
  case UNCORRESPOND:
    _popupPicked->vacateRotoCorr();
    _parent->updateGL();
    break;
  case MOVE_TO_TOP:
    _dc->MoveToTop(_popupPicked);
    break;
  case MOVE_ONE_UP:
    _dc->MoveOneDown(_popupPicked);
    break;
  case MOVE_ONE_DOWN:
    _dc->MoveOneUp(_popupPicked);
    break;
  case MOVE_SPAN_TO_TOP:
    moveSpanToTop(_popupPicked);
    break;
  case MOVE_TO_BOTTOM:
    _dc->MoveToBottom(_popupPicked);
    break;
  case MOVE_SPAN_TO_BOTTOM:
    moveSpanToBottom(_popupPicked);
    break;
  case DELETE:
    _dc->deletePath(_popupPicked);
    if (_selected) {
      _selected = NULL;
      emit selectChanged();
    }
    if (_currDPath == _popupPicked) _currDPath = NULL;
    _popupPicked = NULL;
    _parent->updateGL();
    //_corrShow = 0;
    break;
  case DELETE_SPAN:
    deleteSpan(_popupPicked);
    break;
  case DELETE_PATCH:
    //_dc->deletePatch(_popupPatchPicked);
    //_popupPatchPicked = NULL;
    break;
  case DELETE_ALL:
    _dc->deleteAll();
    if (_currDPath) { delete _currDPath; _currDPath=NULL; }
    if (_selected) {
      _selected = NULL;
      emit selectChanged();
    }
    _parent->updateGL();
    break;
  case ACTIVATE_ALL:
    _dc->setJustDrawn(1);
    break;
  case SELECTD:
    _selected = _popupPicked;
    emit selectChanged();
    _parent->updateGL();
    break;
  case FIXD:
    _popupPicked->setFixed(true);
    break;
  case REBUILD_STROKE_SPAN:
    rebuildStrokeSpan(_popupPicked);
    break;
  case REDO_FORWARD_SPAN:
    redoForwardSpan(_popupPicked);
    break;
  case UNFILL:
    _popupPicked->setFilled(false);
    _popupPicked->freshenAppearance();
    break;
  case SET_THICKNESS:
    bool ok;
    int res = QInputDialog::getInteger("yo","Width:",5,0,100,1,&ok, _parent);
    if (ok)
      _popupPicked->setUniformThickness(res);
    break;

  }
}

void Draw::redoForwardSpan(DrawPath* path) {
  path->vacateRotoCorr();
  path->correspondRoto2(_is->getRotoCurves(_frame), _w, _h);
  forwGenerateSpan(path, _frame);
}

void Draw::deleteSpan(DrawPath* p) {
  
  DrawPath *curr=p, *prev;
  int frame = _frame;

  int a=_is->getStart(),b=_is->getEnd();
  if (!getRange(a,b))
    return;

  while (curr->prevC() != NULL && frame>a) {
    curr = curr->prevC();
    frame--;
  }
  
  prev = curr;
  curr = curr->nextC();
  while (curr!=NULL && frame<b) {
    _is->getDrawCurves(frame)->deletePath(prev);
    prev = curr;
    curr = curr->nextC();
    frame++;
  }
  _is->getDrawCurves(frame)->deletePath(prev);

  if (_selected) {
    _selected = NULL;
    emit selectChanged();
  }
  if (_currDPath == _popupPicked) _currDPath = NULL;
  _popupPicked = NULL;
  //_corrShow = 0;
}

void Draw::moveSpanToTop(DrawPath* p) {
  DrawPath *curr=p;
  int frame = _frame;

  int a=_is->getStart(),b=_is->getEnd();
  if (!getRange(a,b))
    return;
  
  while (curr->prevC() != NULL && frame>a) {
    curr = curr->prevC();
    frame--;
  }

  do {
    _is->getDrawCurves(frame)->MoveToTop(curr);
    curr = curr->nextC();
    frame++;
  } while (curr!=NULL && frame<=b);
}

void Draw::moveSpanToBottom(DrawPath* p) {
  DrawPath *curr=p;
  int frame = _frame;

  int a=_is->getStart(),b=_is->getEnd();
  if (!getRange(a,b))
    return;

  while (curr->prevC() != NULL && frame>a) {
    curr = curr->prevC();
    frame--;
  }

  do {
    _is->getDrawCurves(frame)->MoveToBottom(curr);
    curr = curr->nextC();
    frame++;
  } while (curr!=NULL && frame<=b);
}


float Draw::mapPressure() const {
  float res = _pressure * (_sens+5.f)*.0025;
  //printf("map %f %f %f\n",_pressure, _sens, res);
  //return _pressure/(_sens/* * _zoom*/);
  return res;
}


void Draw::runCan() {
  printf("cannum %d\n",_canNum);

  _canNum++;
  _parent->updateGL();
}



void Draw::copy() {
  if (_tool!=D_SELECT) {
    printf("Must be in select mode to copy\n"); return;
  }
  if (_selected==NULL) {
    printf("Select curve to copy\n"); return;
  }
  _copyPtr = _selected;
  _copyFrame = _frame;
  printf("Copy\n");
}

void Draw::finishJustPasted() {
  if (!_freeDraw) {
    assert(_selected);
    bool res = _selected->correspondRoto2(_is->getRotoCurves(_frame), _w, _h);
    if (res)
    propagatePathEverywhere(_selected);
    _justPasted = false;
  }
}

void Draw::paste() {
  if (_tool!=D_SELECT) {
    printf("Must be in select mode to paste\n");
    return;
  }
  if (!_copyPtr) return;

  printf("Paste\n");
  DrawPath* B = new DrawPath(*_copyPtr);
  _dc->addPath(B);
  _selected = B;
  emit selectChanged();
  _justPasted = true;
  _parent->updateGL();

  /*
  if (_copyFrame == _frame) return;

  DrawPath *A=_copyPtr;
  int numP = A->getNumElements();
  int delta = (_frame > _copyFrame) ? 1: -1, i, j;
  RotoPath** BR = new RotoPath*[numP];
  int *BRi = new int[numP];
  for (j=0; j<numP; j++) {
    BR[j] = A->getRotoCorr(j, BRi[j]);
    if (!BR[j]) {
      printf("No go on the paste 1\n");
      delete[] BR; delete[] BRi; return;
    }
  }

  for (i=_copyFrame; i!=_frame; i+=delta) {
    for (j=0; j<numP; j++) {
      if (delta==1) {
	BR[j] = BR[j]->nextC();
	BRi[j] = BR[j]->getNextCorrIndex(BRi[j]);
      }
      else {
	BR[j] = BR[j]->nextC();
	BRi[j] = BR[j]->getNextCorrIndex(BRi[j]);
      }
      if (!BR[j]) {
	printf("No go on the paste 2\n");
	delete[] BR; delete[] BRi; return;
      }
    }
  }

   DrawPath* B= new DrawPath();
   _is->getDrawCurves(_frame)->addPath(B);
   B->copyLook(A);
   DrawPath::fillInterpolatedCurve(B,A, BR,BRi);
  */
}


void Draw::corrPropAll() {
  _parent->setCursor(Qt::WaitCursor);
  _dc->startDrawPathIterator();
  DrawPath* curr;
  while ((curr=_dc->IterateNext()) != NULL) {
    if (!curr->corrToRoto()) {
      curr->correspondRoto2(_is->getRotoCurves(_frame), _w, _h);
    }
    if (curr->nextC() == NULL || curr->prevC() == NULL)
      propagatePathEverywhere(curr);
  }
  emit restoreCursor();
}

void Draw::setCurrTexture ( int i ) { 

  printf("Draw: setting texture to %d\n", i) ; 
  DrawPath::setCurrTexture(i) ; 
  
} 


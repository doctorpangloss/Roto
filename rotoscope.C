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
#include <qcursor.h>
#include <qinputdialog.h>
#include "rotoscope.h"

/*
Rotoscope::Rotoscope(const QImage im, QWidget* parent) :  _parent(parent) {
  _seedPoint.Set(-1,-1);
  _graph = new ImageGraph(im);
  _wire = NULL;
  _w = im.width();
  _h = im.height();
  _is = NULL;
  _currPath = NULL;  
}
*/

Rotoscope::Rotoscope(ImgSequence* is, QGLWidget* parent) :  
  _parent(parent), _is(is) {
  //_graph = _is->getRelativeImageGraph(0);
  _seedPoint.Set(-1,-1);
  _wire = NULL;
  _w = _is->width();
  _h = _is->height();
  _currPath = NULL;
  //_is->loadRotoCurves();  // should this go here?
  _rc = is->getRelativeRotoCurves(0);
  _frame = is->getStart();
  _corrShow = 0;
  //_selected = NULL;
  _iscissors = false;
  _visible = true;
  //_propFrame = -1;
  _canNum = 0;
  //_spliceFlag=-1;
  //_spliceIndex=-1;
  _showTrackPoints=false;
  _nudgeFreezesCurve = false;
  _lastFrameTouched = -1;
  _changeLastFrameTouched = -1;
  _shadowSelected = NULL;
  _drawing = false;
  clearCorrForw();
  _tool = T_MANUAL;
  _ctrlDrag = NULL; _dragCtrlNum = -1;
  //_ccomp=NULL;
  //_trackData = NULL;  
  _tracking = false;
  _manualStage = 0;
  _regionPicked = NULL;
  _viewRegions = true;
  _visEffort = true;
  mutualInit();
}

void Rotoscope::setImageGraph(ImageGraph* ig) {
  delete _graph;  assert(!ig || (ig->width() == _w && ig->height() == _h));
  _graph = ig;
}


void Rotoscope::timerEvent (QTimerEvent *) {
  assert(_tracking); 
  _parent->updateGL();
}

void Rotoscope::paintGL() {
  //glLineWidth(4); 
  if (_tracking) {
    /*assert(_ccomp && _trackData);
    _trackData->_z_mutex.lock();
    _ccomp->copyLocs(_trackData);
    _trackData->_z_mutex.unlock();
    if (_is->globalTC.finished()) {
      delete _ccomp; _ccomp = NULL;
      delete _trackData; _trackData = NULL;
      _tracking = false;
      killTimer(_trackingTimer);
      }*/

    list<TrackGraph*>::iterator tgc;
    //list<MultiTrackData*>::iterator mtc;
    list<MultiSplineData*>::iterator mtc;
    list<KLT_TrackingContext*>::iterator kc;
    for (tgc=_ccompV.begin(), mtc=_trackDataV.begin(), kc=_TCV.begin();
	 tgc!=_ccompV.end() && mtc!=_trackDataV.end() && kc!=_TCV.end(); ) {

      (*mtc)->_z_mutex->lock();
      if ((*kc)->_stateOk)
	(*tgc)->copyLocs(*mtc);
      //else
      //assert(0);  
      (*mtc)->_z_mutex->unlock();

      if ((*kc)->finished()) {
	delete *tgc; tgc = _ccompV.erase(tgc);
	delete *mtc; mtc = _trackDataV.erase(mtc);
	delete *kc;  kc  = _TCV.erase(kc);
      }
      else {
	++tgc; ++mtc; ++kc;
      }

      if (_ccompV.size()==0) {
	assert(_trackDataV.size()==0 && _TCV.size()==0);
	killTimer(_trackingTimer);
	_tracking = false;
	emit selectChanged();
	break;
      }

    }
  }

  // paint older paths
  _rc->renderAllPaths(_showTrackPoints, _visEffort);

  if (_currPath)
    _currPath->renderHandles();

  // paint previous frame active, if it exists
  /*
  RotoPath* prevp = getPrevActive(NULL);
  if (prevp != NULL) {
    glEnable(GL_LINE_STIPPLE);
    glPointSize(5.0);
    glColor3f(.5,.5,.5);
    prevp->render(false);
    glDisable(GL_LINE_STIPPLE);
    }*/
  
  // print old tracking paths grey
  if (_tool==T_TRACK && _propMode==1) {
    glPointSize(4.0);
    glColor3f(.5,.5,.5);
    if (_is->validFrameNum(_lastFrameTouched) && _lastFrameTouched < _frame) {
      if (!_shadowSelected) {
	glEnable(GL_LINE_STIPPLE);
	_is->getRotoCurves(_lastFrameTouched)->renderFreePaths(false);
	glDisable(GL_LINE_STIPPLE);
      }
      else {
	_shadowSelected->render(false);
	glColor3f(1,0,0);
	_shadowSelected->renderDirection();
	//_shadowSelected->renderSelected();
      }
    }
    glPointSize(2.0);
  }

  // draw correspondence lines
  
  if (_corrShow && _currPath && _currPath->prevC()) {

    glColor3f(0,1,1); 

    assert(_currPath->prevC() && _currPath->prevC()->nextC() == _currPath);
    _currPath->prevC()->renderNextCorrLines();
  }

  // paint current path red
  if (_wire != NULL) {
    glColor3f(1,0,0);
    _wire->startTracePath(_freePoint);
    const Vec2i* loc;
    glBegin(GL_LINE_STRIP);
    while ((loc = _wire->prevStep()) != NULL) 
      glVertex2f(loc->x()+.5, loc->y()+.5);
    glEnd();
  }

  // paint all correspondences
  if (_showAllCorr) {
    glColor3f(0,1,0);
    _rc->renderNextCorrLines();
    /*
    _rc->startCurvesIterator();
    const RotoPath* curr;
    while ((curr=_rc->getNextCurve()) != NULL) {
      if (!curr->nextC()) continue;
      curr->renderNextCorrLines();
    }
    */
  }

  if (_viewRegions)
    _rc->renderRegions(); // make switchable

  // render selection
  if (!_selected.empty()) {
    glColor3f(1,0,0);
    for (PathV::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
      (*c)->renderSelected();
  }

  // render transform selection
  if (_tool == T_TRANSFORM && _transformer.transformable()) {
    glColor3f(1,0,0);
    _transformer.render();
    //for (PathV::const_iterator c = _selected.begin(); c!=_selected.end(); ++c)
    //(*c)->renderSelected(); // G!
  }

  // if endpoint dragging, render endpoint in red
  //if (_ctrlDrag)
  //_ctrlDrag->renderEmphasizedControl(_endDragWhich);

  /*if (_tool==T_OVERDRAW && _selected.size()==1) { // overdraw
    glColor3f(1,0,0);
    (*_selected.begin())->renderDirection();
  }

  if (_tool==T_OVERDRAW && _currPath) { // overdraw
    _currPath->render(false);
    }*/


}


void Rotoscope::restartPausedTracking() {
  
  list<TrackGraph*>::iterator tgc;
  //list<MultiTrackData*>::Iterator mtc;
  list<MultiSplineData*>::iterator mtc;
  list<KLT_TrackingContext*>::iterator kc;
  for (tgc=_ccompV.begin(), mtc=_trackDataV.begin(), kc=_TCV.begin();
       tgc!=_ccompV.end() && mtc!=_trackDataV.end() && kc!=_TCV.end(); 
       ++tgc, ++mtc, ++kc) {
   
    /*MultiTrackData* curr = *mtc;
    QMutex* mut = curr->_z_mutex;
    mut->lock();
    *mtc = new MultiTrackData();
    (*mtc)->transferThreadStuff(curr);
    (*mtc)->_numFrames = curr->_numFrames;
    (*tgc)->startMulti(*mtc);
    (*tgc)->finishStartingMulti(*mtc);
    (*mtc)->transferEdgeMins(curr);
    (*kc)->_mt = *mtc;
    delete curr;
    mut->unlock();
    (*mtc)->_z_wait->wakeAll();*/


    
    MultiSplineData* curr = *mtc;
    QMutex* mut = curr->_z_mutex;
    mut->lock();
    *mtc = new MultiSplineData();
    (*mtc)->transferThreadStuff(curr);
    (*mtc)->_numFrames = curr->_numFrames;
    (*tgc)->buildMulti(*mtc);
    (*mtc)->takeMasks(curr);
    //(*mtc)->transferEdgeMins(curr);  We are building these from scratch, I believe
    (*kc)->_mts = *mtc;
    delete curr;
    mut->unlock();
    (*mtc)->_z_wait->wakeAll();

  }
  
  _trackingTimer = startTimer(200);
}

void Rotoscope::signalTracking() {
  killTimer(_trackingTimer);

  // pause ALL trackers (big stick approach)
  list<TrackGraph*>::iterator tgc;
  list<MultiSplineData*>::iterator mtc;
  //list<MultiTrackData*>::iterator mtc;
  list<KLT_TrackingContext*>::iterator kc;
  for (tgc=_ccompV.begin(), mtc=_trackDataV.begin(), kc=_TCV.begin();
       tgc!=_ccompV.end() && mtc!=_trackDataV.end() && kc!=_TCV.end(); 
       ++tgc, ++mtc, ++kc) {
    (*mtc)->_z_mutex->lock();
    (*kc)->_stateOk=false;
    (*mtc)->_z_mutex->unlock();
  }
}


void Rotoscope::mousePressEvent( QMouseEvent *e ) {

  if (_corrShow) {
    assert(_currPath);
    _currPath = NULL;
    _corrShow = 0;
  }

  if (_tool == T_MANUAL) {
    bool curveStart = (_currPath==NULL);
    if (curveStart) {
      _currPath = new RotoPath();
      _rc->addPath(_currPath);
      _currPath->startManualCreation();
      _manualStage = 0;  
      printf("Manual curve start\n");
    }
    _parent->setMouseTracking(1);
    Vec2f newLoc = unproject(e->x(), e->y(), _h);    

    if (_manualStage == 0) {
      printf("Manual add 2\n");
      _currPath->addCtrl(newLoc);
      _currPath->addCtrl(newLoc);
    }
    else {
      printf("Manual add 3\n");
      _currPath->addCtrl(newLoc);
      _currPath->addCtrl(newLoc);
      _currPath->addCtrl(newLoc);
    }

  }

  else if (_tool==T_DRAW || (_tool==T_TRACK && _shadowSelected)) {
    _drawing = true;
    bool curveStart = (_currPath==NULL);
    if (curveStart) {
      _currPath = new RotoPath();
      _rc->addPath(_currPath);
    
      // deal with rotocurve correspondence
      if (_propMode == 1 && _tool==T_TRACK) 
	  _currPath->startPrevCorrespondence(_shadowSelected);
      _parent->setMouseTracking(1);
    }

    /*if (!_graph)
      setImageGraph(_is->getImageGraph(_frame));

    if (_iscissors && _graph) {
      if (_wire != NULL)
	addPathSegment(); // add current path segment to latest path
      
      Vec2f loc = unproject(e->x(), e->y(), _h);
      _seedPoint.Set((int)loc.x(), (int)loc.y());
      _freePoint = _seedPoint;
      
      if (_wire != NULL) delete _wire;
      _wire = new LiveWire(_seedPoint,_graph);
      }*/
    //else if (!curveStart) // not scissoring, thus click ends curve
    //finishRotoCurve();
    
  }



  /*else if (_tool==T_OVERDRAW) {  // overdraw
    Vec2f loc = unproject(e->x(), e->y(), _h);
    if (_spliceFlag > -1) { // already drawing

    }

    else {  // new overdraw segment
      if (_selected.empty()) {
	printf("Please select a curve to edit, first.\n");
	return;
       }
      
      assert(_currPath==NULL);
      _currPath = new RotoPath();
      _drawing = true;

      double d1 = (*_selected.begin())->distanceTo2(loc, _spliceIndex);
      if (d1<9) { // starting on curve
	_spliceFlag=1;
      }
      else {      // new beginning to curve
	_spliceFlag=0;
	_spliceIndex=0;
      }

 
    }

    } // end overdraw tool*/



  else if (e->button() == Qt::LeftButton && _tool==T_NUDGE) { // nudge

    // for now, implement endpoint mouse dragging
    Vec2f loc = unproject(e->x(), e->y(), _h);
    
    _ctrlDrag = _rc->distanceToCtrls2(loc, &_dragCtrlNum);
    if (_ctrlDrag) {
      int whichEnd;
      if (_ctrlDrag->isCtrlEnd(_dragCtrlNum, &whichEnd))
	_ctrlDrag->fixEndpoint(whichEnd); 
      else /*if (_ctrlDrag->_bez->isCtrlAnyEnd(_dragCtrlNum))*/
	_ctrlDrag->fixInternal(_dragCtrlNum);
      //else
      //_ctrlDrag->setFixed(true);
      _ctrlDrag->setCtrlTouched(_dragCtrlNum, true);
      if (_tracking)
	signalTracking();
      _parent->updateGL();
    }

  } // end nudge

  else if (_tool==T_SELECT) { // select
    
    RotoPath* pick = _rc->pickPrimitive(e->x(), e->y());
    if (pick) {
      if (_shift) {  // multiple selection
	if (find(_selected.begin(), _selected.end(), pick) == _selected.end())
	  _selected.push_back(pick);
      }
      else {         // not multiple selection
	_selected.clear();
	_selected.push_back(pick);
      }
    }
    else {
      _selected.clear();
      emit selectChanged();
    }

    if (_selected.size() > 0) {
      _parent->setMouseTracking(1);
      _dragLoc = unproject(e->x(), e->y(), _h);
    }

    if (_corrForw && _frame-_corrForwFrame == 1 && _selected.size()==1)
      _corrForw->handleCorrForw(*_selected.begin());
    _parent->updateGL();
    _changeLastFrameTouched=_frame;
    emit selectChanged();
  }
  
  else if (_tool == T_TRANSFORM) {
    Vec2f loc = unproject(e->x(), e->y(), _h);
    if (_shift && !_control) { // add to selected
      float which;
      RotoPath* rp = _rc->pickSplineSegment(loc, which);
      if (rp)
	_transformer.addSplineSegment(rp, which);
      if (_tracking)
	signalTracking();
      _parent->updateGL();
    }
    else if (_shift && _control) {
      RotoPath* pick = _rc->pickPrimitive(e->x(), e->y());
      if (pick) {
	_transformer.addSpline(pick);
	if (_tracking)
	  signalTracking();
	_parent->updateGL();
      }
    }
    else {
      if (!_transformer.transformable()) return;
      _transformer.startTransform(loc);
    }
  }
}



void Rotoscope::mouseMoveEvent( QMouseEvent *e ) {
  if (_tool == T_MANUAL && _currPath) {
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    int n = _currPath->getNumControls();
    if (_manualStage == 0)
      _currPath->setControl(newLoc, n-1);
    else {
      _currPath->setControl(newLoc, n-1);
      // set -3 to reflection of -1 across -2
      Vec2f rloc (_currPath->getCtrlRef(n-2), _currPath->getCtrlRef(n-1));
      rloc += _currPath->getCtrlRef(n-2);
      _currPath->setControl(rloc, n-3);
    }
    _parent->updateGL();
  }


  else if (_drawing && (_tool == T_DRAW || /*_tool==T_OVERDRAW ||*/ _tool==T_TRACK) && _currPath) {
    bool renderMe=false;

    if (e->state() & Qt::LeftButton) { // dragging
      if (_iscissors && _graph) {
	addPathSegment();
	delete _wire; _wire = NULL;
	_iscissors = false;
      }
      Vec2f newLoc = unproject(e->x(), e->y(), _h);
      _currPath->addVertex(newLoc.x(), newLoc.y());
      renderMe=true;
    }
    else {
      if (_iscissors && _graph) {
	assert(_wire);
	assert(_seedPoint.x() != -1 && _wire != NULL);  // should be wiring
	Vec2f loc = unproject(e->x(), e->y(), _h);
	_freePoint.Set((int)loc.x(), (int)loc.y());
	_freePoint.Clamp(0,_w-1,_h-1);
	_wire->registerInterest(_freePoint, _parent);
	renderMe=true;
      }
    }
  

    /*
    if (_iscissors && _graph) {
      assert(_wire);
      assert(_seedPoint.x() != -1 && _wire != NULL);  // should be wiring
      Vec2f loc = unproject(e->x(), e->y(), _h);
      _freePoint.Set((int)loc.x(), (int)loc.y());
       _freePoint.Clamp(0,_w-1,_h-1);
      _wire->registerInterest(_freePoint, _parent);
      renderMe=true;
    }
    else if (_currPath) {
      Vec2f newLoc = unproject(e->x(), e->y(), _h);
      _currPath->addVertex(newLoc.x(), newLoc.y());
      renderMe=true;
    }
    else
      printf("Why mouse move\n");
    */

    if (renderMe) {
      QPoint curs = _parent->mapFromGlobal(QCursor::pos());
      if ((e->x() - curs.x())*(e->x() - curs.x()) +
	  (e->y() - curs.y())*(e->y() - curs.y()) < 25)
	_parent->updateGL();
    }
  }

  else if (_tool==T_NUDGE && _ctrlDrag!=NULL) { // endpoint nudge
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    _ctrlDrag->setControl(newLoc, _dragCtrlNum);
    //modifyEnd(_endDragWhich, newLoc);
    int whichEnd;
    if (_ctrlDrag->isCtrlEnd(_dragCtrlNum, &whichEnd))
      _ctrlDrag->reconcileOneJointToMe(whichEnd);
    //_endDrag->reconcileOneJointToMe(_endDragWhich); 
    _ctrlDrag->handleNewBezCtrls();
    _parent->updateGL();
  }

  else if (_tool == T_SELECT && _selected.size() > 0) {
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    Vec2f delta(newLoc, _dragLoc);
    PathV::iterator i;
    for (i=_selected.begin(); i!=_selected.end(); ++i) {
      (*i)->translate(delta);
      //(*i)->handleNewBezCtrls();
    }
    _dragLoc = newLoc;
    _parent->updateGL();
  }

  else if (_tool == T_TRANSFORM && _transformer.transforming()) {
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    _transformer.handleNewLoc(newLoc);
    _parent->updateGL();
  }
}




void Rotoscope::mouseReleaseEvent( QMouseEvent *e ) {

  if (_tool == T_MANUAL && _currPath) {
    _parent->setMouseTracking(0);
    Vec2f newLoc = unproject(e->x(), e->y(), _h);

    if (_manualStage == 0) {
      //_currPath->addCtrl(newLoc);
      _manualStage=1;
    }

    else {
      //_currPath->setControl(newLoc, _currPath->getNumControls()-2);
      //_manualStage=0;
    }
  }

  else if (_drawing && (_tool == T_DRAW || _tool==T_TRACK) && _currPath) {
    if (!_iscissors) {
      _iscissors = true;
      assert(_wire == NULL);
      Vec2f loc = unproject(e->x(), e->y(), _h);
      _seedPoint.Set((int)loc.x(), (int)loc.y());
      //QPoint curs = _parent->mapFromGlobal(QCursor::pos());
      //_seedPoint.Set(curs.x(),curs.y());
      _freePoint = _seedPoint;
      if (!_graph)
	setImageGraph(_is->getImageGraph(_frame));
      _wire = new LiveWire(_seedPoint,_graph);
    }
    else {
      if (!_graph)
	setImageGraph(_is->getImageGraph(_frame));
      
      if (_iscissors && _graph) {
	if (_wire != NULL)
	  addPathSegment(); // add current path segment to latest path
	
	Vec2f loc = unproject(e->x(), e->y(), _h);
	_seedPoint.Set((int)loc.x(), (int)loc.y());
	_freePoint = _seedPoint;
	
	if (_wire != NULL) delete _wire;
	_wire = new LiveWire(_seedPoint,_graph);
      }
    }
  }

  /*else if (_drawing && (_tool == T_OVERDRAW)) {
    
    Vec2f loc = unproject(e->x(), e->y(), _h);
    if (_spliceFlag > -1) { // already drawing

	if (_wire != NULL)
	  addPathSegment();
	_seedPoint.Set((int)loc.x(), (int)loc.y());
	_freePoint = _seedPoint;
	if (_wire != NULL) delete _wire;
	if (!_graph)
	  setImageGraph(_is->getImageGraph(_frame));
	_wire = new LiveWire(_seedPoint,_graph);
	_iscissors = true;
	_parent->setMouseTracking(1);            
    }

    else {  // new overdraw segment
      if (_selected.empty()) {
	printf("Please select a curve to edit, first.\n");
	return;
      }
      

      if (!_graph)
	setImageGraph(_is->getImageGraph(_frame));


      
      if (_graph) {
	assert(_currPath);
	_seedPoint.Set((int)loc.x(), (int)loc.y());
	_freePoint = _seedPoint;
	if (_wire != NULL) delete _wire;
	_wire = new LiveWire(_seedPoint,_graph);
	_iscissors = true;
      }
      _parent->setMouseTracking(1);            
    }


    }*/
  

  else if (_tool == T_TRACK  && _propMode==1) { // Track
    if (_is->validFrameNum(_lastFrameTouched) && _lastFrameTouched < _frame) {
      if (_shadowSelected==NULL) {
	_shadowSelected = _is->getRotoCurves(_lastFrameTouched)->pickFreePrimitive(e->x(), e->y());
	if (_shadowSelected)
	  _parent->updateGL();
      }
    }
  }
  else if (_ctrlDrag != NULL &&_tool==T_NUDGE) { // endpoint nudge
    //_spliceIndex=-1;
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    _ctrlDrag->setControl(newLoc, _dragCtrlNum);
    int whichEnd;
    if (_ctrlDrag->isCtrlEnd(_dragCtrlNum, &whichEnd))
      _ctrlDrag->reconcileOneJointToMe(whichEnd);
    //_crtlDrag->modifyEnd(_endDragWhich, newLoc);
    //_endDrag->reconcileOneJointToMe(_endDragWhich); 
    _ctrlDrag->handleNewBezCtrls();
    if (_tracking)
      restartPausedTracking();
    _ctrlDrag = NULL;
    _parent->updateGL();
  }

  else if (_tool == T_SELECT && _selected.size() > 0) {
    PathV::iterator i;
    for (i=_selected.begin(); i!=_selected.end(); ++i) {
      (*i)->handleNewBezCtrls();
    }
    //for (int i=0; i<_selected.size(); ++i)
    // _selected[i]->handleNewBezCtrls();
    _parent->setMouseTracking(0);
    _parent->updateGL();
  }

  else if (_tool == T_TRANSFORM && _transformer.transforming()) {
    Vec2f newLoc = unproject(e->x(), e->y(), _h);
    _transformer.endTransform(newLoc);    
    _parent->updateGL();
  }
  
}


void Rotoscope::keyPressEvent ( QKeyEvent *e ) {
  if (_tool==T_DRAW || _tool==T_TRACK) {  // draw or track
    if (e->key() == Qt::Key_Return) {  // close path
      _drawing = false;
      if (!_currPath || _currPath->getNumElements()==0) return;
      Vec2f beg = _currPath->getElement(0);
      _freePoint.Set((int)beg.x(), (int)beg.y());
      if (_wire == NULL) {
	QPoint curs = _parent->mapFromGlobal(QCursor::pos());
	_seedPoint.Set(curs.x(),curs.y());
	_wire = new LiveWire(_seedPoint,_graph);
      }
      addPathSegment();
      _currPath->addVertex(float(_freePoint.x())+.5, float(_freePoint.y())+.5);
      finishRotoCurve();
    }
    else if (e->key() == Qt::Key_Backspace && _currPath) { // end path
      _drawing = false;
      _iscissors =false;
      finishRotoCurve();
    }
    /*  else if (e->key() == Qt::Key_S) {
	_rc->getTail()->fair();
	printf("smoothed\n");
    _parent->updateGL();
    }*/
    else if (e->key() == Qt::Key_C) {
      killAnyRoto();
      _iscissors =false;
      _parent->updateGL();
    }
    /*else if (e->key() == Qt::Key_I) {
      _iscissors = !_iscissors;
      if (!_iscissors && _wire != NULL) {
	addPathSegment();
	delete _wire; _wire = NULL;
      }
      if (_iscissors && _wire == NULL && _graph) {
	QPoint curs = _parent->mapFromGlobal(QCursor::pos());
	_seedPoint.Set(curs.x(),curs.y());
	_freePoint = _seedPoint;
	_wire = new LiveWire(_seedPoint,_graph);
      }
    } */
    else
      e->ignore();
  }
  else if (_tool == T_SELECT && _selected.size()>1) {  // select
    if (e->key() == Qt::Key_Backspace) {
      for (PathV::const_iterator c = _selected.begin(); c!=_selected.end(); ++c) 
	_rc->deletePath(*c);
      _selected.clear();
      emit selectChanged();
      _parent->updateGL();
    }
    else if (e->key() == Qt::Key_N && _selected.size()==1) {
      RotoPath* s = _rc->cycle(*(_selected.begin()));
      _selected.clear(); _selected.push_back(s);
      emit selectChanged();
      _parent->updateGL();
    }
    else
      e->ignore();
  }


  else
      e->ignore();

}

#define ABSMIN(a,b) ((fabs(a))<(fabs(b)) ? (a) : (b))
void Rotoscope::intrinsicBlend(Vec2f* out, Vec2f* lengths, Vec2f* angles, int n, double t) {
  assert(n>0);
  double length = lengths[n-1].x() + t*(lengths[n-1].y() - lengths[n-1].x());
  double b=angles[n-1].y(), a=angles[n-1].x();
  double angle = a + t* ABSMIN(b-a, b - ( (a>0) ? (a-2.*M_PI) : (a+2.*M_PI) ) );
  //printf("%f t from %f to %f yields %f\n", t, a*180./M_PI, b*180./M_PI, angle*180./M_PI);
  out[n] = out[n-1];
  Vec2f unit(cos(angle), sin(angle));
  unit *= (float) length;
  out[n] += unit;
}

void Rotoscope::intrinsicBlend2(Vec2f* out, const Vec2f* earlier, const Vec2f& length, 
				const Vec2f& angle, double t) {
  double l = length.x() + t*(length.y() - length.x());
  double b=angle.y(), a=angle.x();
  double an = a + t* ABSMIN(b-a, b - ( (a>0) ? (a-2.*M_PI) : (a+2.*M_PI) ) );
  //printf("%f t from %f to %f yields %f\n", t, a*180./M_PI, b*180./M_PI, angle*180./M_PI);
  *out = *earlier;
  Vec2f unit(cos(an), sin(an));
  unit *= (float) l;
  *out += unit;
}

void Rotoscope::trackForward() {
  PathV::const_iterator c;
  int i = _frame+1;
  c = _selected.begin();
  //for (c=_selected.begin(); c!=_selected.end(); ++c) {
    (*c)->takeBez((*c)->prevC()->getBez());
    (*c)->handleNewBezCtrls();
    RotoPath* newpath = new RotoPath(**c);
    newpath->buildTouched(false);
    newpath->setFixed(false);
    _is->getRotoCurves(i)->addPath(newpath);
    (*c)->setNextC(newpath);
    (*c)->buildINextCorrs();
    newpath->setPrevC(*c);
    newpath->buildIPrevCorrs();
    _toTrack.push_back(newpath);
    //}

  RotoPath::buildForwardReconcileJoints(_selected, _frame, i);
  performTracks(_frame-1,i,false,true);

  newpath->takeBez((*c)->getBez());
  newpath->handleNewBezCtrls();
  //for (c=_selected.begin(); c!=_selected.end(); ++c) {
    //(*c)->nextC()->takeBez((*c)->getBez());
    //(*c)->nextC()->handleNewBezCtrls();
    (*c)->setFixed(true);
    //}
}

void Rotoscope::rejigWrapper(bool startOver) {
  if (_selected.empty()) return;
  PathV paths, oPaths;
  //RotoPath::copySetToVector(&paths, &_selected);
  paths = _selected;
  unsigned int i;
  bool fixed = false;

  for (i=0; i<paths.size(); ++i) {
    if (paths[i]->fixed())
      fixed = true;
  }

  if (!fixed) {
    rejig(&paths, _frame, startOver);
  }
  else {
    oPaths = paths;
    bool good=true;
    for (i=0; i<oPaths.size(); ++i) {
      oPaths[i] = oPaths[i]->prevC();
      if (oPaths[i]==NULL)
	good=false;
    }
    if (good)
      rejig(&oPaths, _frame-1, startOver);
    
    good=true;
    oPaths = paths;
    for (i=0; i<oPaths.size(); ++i) {
      oPaths[i] = oPaths[i]->nextC();
      if (oPaths[i]==NULL)
	good=false;
    }
    if (good)
      rejig(&oPaths, _frame+1, startOver);
  }

  _selected.clear();
  emit selectChanged();
}


void Rotoscope::rejig(PathV* paths, int frame, bool startOver) {
  if (paths->empty()) return;
  PathV aPaths, bPaths;
  aPaths = *paths;
  bPaths = aPaths;
  int aFrame = frame, bFrame = frame, numFrames;
  bool keepGoing = true;
  PathV::const_iterator c, c2;
  RotoPath *tmp; int i;

  while (keepGoing) {
    for (c = aPaths.begin(); c!=aPaths.end(); ++c) {
      if ((*c)->fixed() || (*c)->prevC()==NULL)
	keepGoing=false;
    }
    if (keepGoing) {  // we're going back in time, step pointers, delete the left behind paths
      for (c = aPaths.begin(), i=0; c!=aPaths.end(); ++c, ++i) {
	tmp = *c;
	aPaths[i] = tmp->prevC();
	//if (aFrame!=frame && startOver)  // don't delete frame paths, since next loop forwards will
	//_is->getRotoCurves(aFrame)->deletePath(tmp);
      }
      --aFrame;
    }
  }
  printf("aFrame is %d\n",aFrame);
  assert(aPaths.size()==paths->size());

  keepGoing = true;
  
  while (keepGoing) {
    for (c = bPaths.begin(); c!=bPaths.end(); ++c) {
      if ((*c)->fixed() || (*c)->nextC()==NULL)
	keepGoing=false;
    }
    if (keepGoing) {  // we're going forward in time, step pointers, delete the left behind paths
      for (c = bPaths.begin(), i=0; c!=bPaths.end(); ++c, ++i) {
	tmp = *c;
	bPaths[i] = tmp->nextC();
	//if (startOver)
	//_is->getRotoCurves(bFrame)->deletePath(tmp);
      }
      ++bFrame;
    }
  }

  printf("bFrame is %d\n",bFrame);
  assert(bPaths.size()==paths->size());
  numFrames = bFrame - aFrame;
  if (numFrames < 2) {
    printf("Too short a span\n");
    _toTrack.clear();
    return;
  }

  assert(_toTrack.empty());
  /*
  if (startOver) {
    for (c = aPaths.begin(); c!=aPaths.end(); ++c)
      (*c)->vacateNextRotoCorr();
    for (c = bPaths.begin(); c!=bPaths.end(); ++c)
      (*c)->vacatePrevRotoCorr();
  }
  */
  for (c = aPaths.begin(), c2 = bPaths.begin(); c!=aPaths.end(); ++c, ++c2) {
    /*
    if (startOver) {
      if (_tool != T_MANUAL) { //G! Need better interface to control this
	(*c2)->startPrevCorrespondence(*c);
	(*c2)->finishPrevCorrespondence();
      }
      else {
	(*c2)->startPrevCorrespondence(*c);
	(*c)->buildINextCorrs();
	(*c2)->buildIPrevCorrs();
      }
    }
    */
    _toTrack.push_back(*c2);
  }
  //performTracks(aFrame, bFrame, startOver);

  if (startOver) // redoTrack
    performTracks(aFrame, bFrame, true, true);
  else
    performTracks(aFrame, bFrame, false, true);

}


void Rotoscope::performTracks(const int aFrame, const int bFrame, bool doInterp, bool useExistingInbetweens) {
  printf("toTrack size %d, a %d b %d\n",_toTrack.size(), aFrame, bFrame);
  bool* done = new bool[_toTrack.size()];
  memset(done,0,_toTrack.size()*sizeof(bool));
  PathV::const_iterator c2;
  PathV::iterator c;

  
  // interpolate somehow
  if (doInterp && !useExistingInbetweens) {
    for (c2 = _toTrack.begin(); c2 != _toTrack.end(); ++c2) {
      keyframeSedInterp((*c2)->prevC(), aFrame, *c2, bFrame);
    }
    printf("Interpolated\n");
  }
  
  
  // build & reconcile joints
  if (doInterp && !useExistingInbetweens)
    RotoPath::buildBackReconcileJoints(_toTrack,bFrame,aFrame);
    //ccomp->buildBackReconcileJoints(bFrame, aFrame);
  printf("Joints done\n");
  

  int i=0, j;
  for (c = _toTrack.begin(); c!= _toTrack.end(); ++c, ++i) {

    if (!done[i]) {
      // get a full list beg to end of linked up paths (stop at already interpolated ones,
      // ones not in toTrack list)
      TrackGraph* ccomp = new TrackGraph();
      MultiSplineData* mts = new MultiSplineData();

      (*c)->buildccomp(ccomp, &_toTrack);

      // start multiTrack
      mts->_numFrames = bFrame - aFrame;
      printf("started multi\n");

      ccomp->buildMulti(mts);
      
      addMasksToMulti(mts, ccomp->getKey0Paths(), aFrame);

      // track
      // pyrms and pyrmsE get deleted by longCurveTrack2
      const KLT_FullCPyramid** pyrms = new (const KLT_FullCPyramid*)[mts->_numFrames+1];
      const KLT_FullPyramid** pyrmsE = new (const KLT_FullPyramid*) [mts->_numFrames+1];

      for (j=0; j<=mts->_numFrames; j++) {
	if (_is->globalTC.useImage) {
	  pyrms[j] = _is->getFullCPyramid(aFrame+j);   
	  assert(pyrms[j]);
	}
	else
	  pyrms[j] = NULL;
	if (_is->globalTC.useEdges) {
	  pyrmsE[j] = _is->getFullEdgePyramid(aFrame+j); 
	  assert(pyrmsE[j]);
	}
	else
	  pyrmsE[j] = NULL;
      }

      KLT_TrackingContext* tc = new KLT_TrackingContext(); // transfer global track settings
      tc->copySettings(&(_is->globalTC));
      //tc->setupLongCurveTrack2(pyrms, pyrmsE, trackData, doInterp); 
      tc->setupSplineTrack(pyrms, pyrmsE, mts, doInterp);

      _ccompV.push_back(ccomp);
      _trackDataV.push_back(mts);
      _TCV.push_back(tc);

      
      tc->start();  // ONE
      //tc->runNoThread();  // For profiling // ONE
      //paintGL();
      //ccomp->copyLocs(mts); // ONE
      //std::exit(0);  
     
      
      if (!_tracking) {
	_tracking = true; 
	emit selectChanged(); 
	_trackingTimer = startTimer(200); 
      }
      
      // copy over results
      //_ccomp->copyLocs(_trackData); // Not for tracking
      
      // Let done know these are processed
      for (c2 = ccomp->paths()->begin(); c2 != ccomp->paths()->end(); ++c2) {
	PathV::const_iterator which = find(_toTrack.begin(), _toTrack.end(), *c2);
	assert(which != _toTrack.end());
	done[which - _toTrack.begin()] = true;
      }
      
      
  // ONE
      /* 
      delete ccomp; _ccompV.clear();
      delete mts; _trackDataV.clear();
      delete tc; _TCV.clear();
      */
    }
  }

  
  _toTrack.clear();
  delete[] done;
  //delete _trackData; delete _ccomp;
}

void Rotoscope::addMasksToMulti(MultiSplineData* mts, const PathV& key0, const int frame0) {
  PathV currPaths = key0;
  PathV::iterator c;
  mts->setMaskDims(_w, _h);
  for (int i=0; i<mts->_numFrames+1; ++i) { // iterate over frames
    unsigned char* mask = _is->getRotoCurves(frame0+i)->makeCummMask(currPaths);
    mts->addMask(mask);
    for (c=currPaths.begin(); c!=currPaths.end(); ++c) // advance paths 1 frame
      (*c) = (*c)->nextC();
  }
}


void Rotoscope::copySplinesAcrossTime() {

  int a=_frame+1,b=_is->getEnd(), i, j;
  if (!getRange(a,b))
    return;
  if (a!=_frame+1) return;

  PathV::const_iterator c;
  PathV::const_iterator v;
  PathV prev;
  for (c=_selected.begin(); c!=_selected.end(); ++c)
    prev.push_back(*c);
  PathV currPaths = prev, next = prev;

  // copy forward
  //a = MAX(a,_frame+1);
  for (i=a; i<=b; ++i) {
    for (v=prev.begin(), j=0; v!=prev.end(); ++v, ++j) {
      RotoPath* newpath = new RotoPath(**v);
      newpath->buildTouched(false);
      newpath->setFixed(false);
      _is->getRotoCurves(i)->addPath(newpath);
      (*v)->setNextC(newpath);
      (*v)->buildINextCorrs();
      newpath->setPrevC(*v);
      newpath->buildIPrevCorrs();
      prev[j] = newpath;
    }
  }
  
  RotoPath::buildForwardReconcileJoints(currPaths, _frame, b);

}

void Rotoscope::keyframeSedInterp(RotoPath* aPath, int aFrame, RotoPath *bPath, int bFrame) { 
  assert(aPath && bPath);
  assert(_propMode);
  int numFrames = bFrame - aFrame;
  assert(numFrames > 0);
  //int numP = aPath->getNumElements();
  int t;

   bPath->setLowHeight(aPath->lowHeight());
   bPath->setHighHeight(aPath->highHeight());
   bPath->setTrackEdges(aPath->trackEdges());
   printf("KeyframeSedInterp\n");


  RotoPath* prev = aPath;
  for (t=1; t < numFrames; ++t) {
    RotoPath *newpath = new RotoPath(*aPath);
    //newpath->setCan(t); 
    newpath->setFixed(false);
    newpath->buildTouched(false);
    _is->getRotoCurves(aFrame+t)->addPath(newpath);
    newpath->setPrevC(prev);
    prev->setNextC(newpath);
    if (t==0)
      newpath->buildINextCorrs();
    else if (t==numFrames-1 && _is->globalTC.pinLast) {
      //newpath->_nextCoors = aPath->_nextCoors;
      newpath->takeNextCont(aPath);
      newpath->buildIPrevCorrs();
    }
    else {
      newpath->buildINextCorrs();
      newpath->buildIPrevCorrs();
    }
    prev=newpath;
  }

  bPath->setPrevC(prev);
  prev->setNextC(bPath);
  //aPath->forgetNextCont();
  aPath->buildINextCorrs();
  assert(aPath->nextC());

}

/*


void Rotoscope::keyframeSedInterp(RotoPath* aPath, int aFrame, RotoPath *bPath, int bFrame) { 
  assert(aPath && bPath);
  assert(_propMode);
  int numFrames = bFrame - aFrame;
  assert(numFrames > 0);
  int numP = aPath->getNumElements();
  int n, t;

   bPath->setLowHeight(aPath->lowHeight());
   bPath->setHighHeight(aPath->highHeight());
   bPath->setTrackEdges(aPath->trackEdges());

  float fnumFrames = float(numFrames);
  Vec2f *lengths = new Vec2f[numP-1], *angles = new Vec2f[numP-1];
  const Vec2f *aP = aPath->getData(), *bP, *bP1; 
  for (n=0; n<numP-1; ++n, ++aP) {
    bP =  bPath->getPointerToElement(aPath->getNextCorrIndex(n));
    bP1 = bPath->getPointerToElement(aPath->getNextCorrIndex(n+1));
    Vec2f vec1(*(aP+1), *aP), vec2(*bP1, *bP);
    lengths[n].Set(vec1.Length(), vec2.Length());

    int bPi = aPath->getNextCorrIndex(n+1);
    while (vec2.x()==0 && vec2.y()==0 && bPi>0)
      Vec2f_Sub(vec2, *bP1, bPath->getElement(--bPi));
    angles[n].Set(atan2(vec1.y(), vec1.x()), atan2(vec2.y(), vec2.x()));
  }

  RotoPath* prev = aPath;
  for (t=1; t < numFrames; ++t) {
    RotoPath *newpath = new RotoPath(*aPath);
    newpath->setFixed(false);
    _is->getRotoCurves(aFrame+t)->addPath(newpath);
    newpath->setPrevC(prev);
    prev->setNextC(newpath);
    if (t==0)
      newpath->buildINextCoors();
    else if (t==numFrames-1 && _is->globalTC.pinLast) {
      newpath->_nextCoors = aPath->_nextCoors;
      newpath->buildIPrevCoors();
    }
    else {
      newpath->buildINextCoors();
      newpath->buildIPrevCoors();
    }
    prev=newpath;

    double s = float(t)/fnumFrames; 
    Vec2f_LinInterp(*(newpath->getPointerToElement(0)), *(aPath->getPointerToElement(0)), 
		    *(bPath->getPointerToElement(0)), s);
    for (n=1; n<numP; ++n) {
      intrinsicBlend2(newpath->getPointerToElement(n), newpath->getPointerToElement(n-1),
		      lengths[n-1], angles[n-1], s);
    }
  }

  bPath->setPrevC(prev);
  prev->setNextC(bPath);
  aPath->_nextCoors = NULL;
  aPath->buildINextCoors();
  assert(aPath->nextC());
  delete[] lengths; delete[] angles;
}

*/

/*
void Rotoscope::keyframeInterpolatePath(RotoPath* aPath, int aFrame, RotoPath *bPath, int bFrame) { 
  _parent->setCursor(Qt::WaitCursor);
  assert(aPath && bPath);
  assert(_propMode);
  int numFrames = bFrame - aFrame;
  assert(numFrames > 0);
  int numP = aPath->getNumElements();

  int n, t, offset,offset2;
  const KLT_FullCPyramid** pyrms = new (const KLT_FullCPyramid*)[numFrames+1];
  float fnumFrames = float(numFrames);
  Vec2f *Plocs = new Vec2f[numP];
  Vec2f *Z = new Vec2f[numP*numFrames];

  if (_is->globalTC.useImage) {
    pyrms[0] = _is->getFullCPyramid(aFrame);
    pyrms[numFrames] = _is->getFullCPyramid(bFrame);
  }
  offset = numP*(numFrames-1);
  for (n=0; n<numP; n++,offset++) {
    Plocs[n] = aPath->getElement(n);
    Z[offset] = bPath->getElement(aPath->getNextCorrIndex(n));
  }

  // calculate lengths, angles (from x-axis, clockwise?)
  offset2 = numP*(numFrames-1);
  Vec2f *lengths = new Vec2f[numP-1], *angles = new Vec2f[numP-1];
  for (n=0; n<numP-1; ++n, ++offset2) {

    Vec2f vec1(Plocs[n+1],Plocs[n]), vec2(Z[offset2+1], Z[offset2]);
    lengths[n].Set(vec1.Length(), vec2.Length());

    t=offset2;
    while (vec2.x()==0 && vec2.y()==0 && t>0) {  // deal with 2 A points mapping to same B point
      Vec2f_Sub(vec2, Z[offset2+1], Z[--t]);
    }
    angles[n].Set(atan2(vec1.y(), vec1.x()), atan2(vec2.y(), vec2.x()));
  }
  
  offset=0;     
  RotoPath* prev = aPath;
  for (t=1; t < numFrames; ++t) {
    if (_is->globalTC.useImage)
      pyrms[t] = _is->getFullCPyramid(aFrame+t);
      
    RotoPath *newpath = new RotoPath(*aPath);
    newpath->setFixed(false);
    _is->getRotoCurves(aFrame+t)->addPath(newpath);
    newpath->setPrevC(prev);
    prev->setNextC(newpath);
    if (t==0)
      newpath->buildINextCoors();
    else if (t==numFrames-1 && _is->globalTC.pinLast) {
      newpath->_nextCoors = aPath->_nextCoors;
      newpath->buildIPrevCoors();
    }
    else {
      newpath->buildINextCoors();
      newpath->buildIPrevCoors();
    }
    prev=newpath;

    double s = float(t)/fnumFrames; 
    Vec2f* startArray = Z+offset;
    offset2 = numP*(numFrames-1);
    for (n=0; n<numP; ++n,++offset,++offset2) {
      //Vec2f_LinInterp(Z[offset], Plocs[n],Z[offset2], s);
      if (n==0)
	Vec2f_LinInterp(Z[offset], Plocs[n],Z[offset2], s);
      else
	intrinsicBlend(startArray, lengths, angles, n, s);
      newpath->actualSetLoc(n,Z[offset]);
      //assert(n==0 || Z[offset-1] != Z[offset]);      
    }

  }
  bPath->setPrevC(prev);
  prev->setNextC(bPath);
  aPath->_nextCoors = NULL;
  aPath->buildINextCoors();

  
  //FILE* fp = fopen("Z","w");   
  //fwrite(Z,sizeof(Vec2f),numP*numFrames,fp);
  //fclose(fp);
  

  int res = _is->globalTC.longCurveTrack( pyrms,numFrames,numP, aPath->trackWidth(),
					  aPath->lowHeight(), aPath->highHeight(), 
					  Plocs, Z );
  
  if (res==KLT_TRACKED || res == KLT_MAX_ITERATIONS) {
    //fprintf(stderr,"final translations\n");
    offset=0;
    RotoPath *curr=aPath;
    for (t=1; t < numFrames; t++) {
      curr = curr->nextC();
      assert(curr);
      for (n=0; n<numP; n++,offset++) {
	curr->actualSetLoc(n, Z[offset]);
      }
    }
    if (!_is->globalTC.pinLast) {
      bPath->resetElements();
      for (n=0; n<numP; n++,offset++) {
	bPath->add(Z[offset]);
      }
    }
    //fprintf(stderr,"\n");
  }
  else
    fprintf(stderr,"Failed to track\n");
  
  delete[] pyrms;
  delete[] Plocs;
  delete[] Z;
  delete[] lengths; delete[] angles;
  _parent->setCursor(Qt::ArrowCursor);
}
*/

/*
void Rotoscope::longTrackSpan(RotoPath* rp) { 
  assert(rp);
  if (rp->fixed()) {
    printf("This curve has been fixed\n");
    return;
  }
  _parent->setCursor(Qt::WaitCursor);
  RotoPath *aPath, *bPath, *curr = rp, *prev=NULL;
  int aFrame = _frame, bFrame = _frame, numFrames;

  printf("back\n");
  while (!curr->fixed() && curr->prevC() != NULL) {
    printf("frame %d, %d points\n",aFrame, curr->getNumElements());
    curr = curr->prevC();
    aFrame--;
  }
  aPath = curr;
  printf("frame %d, %d points\n",aFrame, curr->getNumElements());
  assert(aPath->getNumElements()==rp->getNumElements());

  curr = rp; 
  printf("forwards\n");
  while (!curr->fixed() && curr->nextC() != NULL) {
    printf("frame %d, %d points\n",bFrame, curr->getNumElements());
    curr = curr->nextC();
    bFrame++;
  }
  bPath = curr;
  printf("frame %d, %d points\n",bFrame, curr->getNumElements());
  numFrames = bFrame - aFrame;
  if (numFrames < 1) {
    printf("Too short a span\n");
    return;
  }
  printf("aframe: %d, bframe: %d\n",aFrame, bFrame);


  int numP = aPath->getNumElements();

  int n, t, offset;
  const KLT_FullCPyramid** pyrms = new (const KLT_FullCPyramid*)[numFrames+1];
  Vec2f *Plocs = new Vec2f[numP];
  Vec2f *Z = new Vec2f[numP*numFrames];
  pyrms[0] = _is->getFullCPyramid(aFrame);
  //pyrms[numFrames] = _is->getFullCPyramid(bFrame);

  offset = numP*(numFrames-1);
  for (n=0; n<numP; n++,offset++) 
    Plocs[n] = aPath->getElement(n);

  offset=0; 
  curr = aPath;
  for (t=1; t <= numFrames; t++) {
    pyrms[t] = _is->getFullCPyramid(aFrame+t);  
    prev=curr;
    curr = curr->nextC();
    assert(curr && prev && prev->_nextCoors);
    printf("%d %d %d\n",numP, curr->getNumElements(), prev->getNumElements());
    assert(prev->getNumElements()==numP);

    for (n=0; n<numP; n++,offset++)
      Z[offset] = curr->getElement(prev->getNextCorrIndex(n));
  }

  int res = _is->globalTC.longCurveTrack( pyrms,numFrames,numP, aPath->trackWidth(),
					  aPath->lowHeight(), aPath->highHeight(), 
					  Plocs, Z );
 
  if (res==KLT_TRACKED || res == KLT_MAX_ITERATIONS) {
    //fprintf(stderr,"final translations\n");
    offset=0;
    RotoPath *curr=aPath;
    for (t=1; t < numFrames; t++) {
      curr = curr->nextC();
      assert(curr);
      for (n=0; n<numP; n++,offset++) {
	//fprintf(stderr,"%.3f %.3f, ",Z[offset].x(), Z[offset].y());
	curr->actualSetLoc(n, Z[offset]);
      }
      //fprintf(stderr,"\n");
    }
    if (!_is->globalTC.pinLast) {
      bPath->resetElements();
      for (n=0; n<numP; n++,offset++) {
	//Vec2f v = Plocs[n];
	//v += Z[offset];
	//fprintf(stderr,"%.3f %.3f, ",Z[offset].x(), Z[offset].y());
	bPath->add(Z[offset]);
      }
    }
    //fprintf(stderr,"\n");
  }
  else
    fprintf(stderr,"Failed to track patch\n");
  
  delete[] pyrms;
  delete[] Plocs;
  delete[] Z;
  _parent->setCursor(Qt::ArrowCursor);
}
*/
/*
RotoPath* Rotoscope::getPrevActive(int* frame) { 
  int curr = _frame;
  RotoPath* active=NULL;
  do {
    curr--;
  } while (curr >= _is->getStart() &&
	   (active = _is->getRotoCurves(curr)->getActive()) == NULL);
  if (active && frame)
    *frame = curr;
  return active;
}
*/

void Rotoscope::finishManualRotoCurve() {
  _currPath->startFinishManualCreation();
  _rc->setupJoints(_currPath);
  _currPath->finishManualCreation();
  _currPath = NULL;
}


void Rotoscope::finishRotoCurve() {
  _parent->setMouseTracking(0);
  _seedPoint.Set(-1,-1);

  if (_currPath->getNumElements()>2) {


    //_rc->printJoints();


    _currPath->fair();
    assert(_currPath->okEnds());
    _currPath->goBezier(12., NULL);  

    // set up joints
    _rc->setupJoints(_currPath);
    _currPath->fillFromBez();
    

    
    if (_currPath->prevC()) {
      _corrShow = 1;
      assert(_is->validFrameNum(_lastFrameTouched));
      _currPath->finishPrevCorrespondence();
      _parent->updateGL();
      if (_propMode == 1)
	_toTrack.push_back(_currPath);
	//keyframeInterpolatePath(_currPath->prevC(), _lastFrameTouched, _currPath, _frame); 
    }
    else
      _currPath = NULL;
  
    //if (_is->getRotoCurves(_propFrame) != NULL) 
    //_is->getRotoCurves(_propFrame)->setActive(NULL);
    //_propFrame = -1;

    _changeLastFrameTouched = _frame;
  }
  else {
    _currPath->clearPrevCorrespondence();
    _rc->deletePath(_currPath);
    _currPath=NULL;
  }
  
  if (_wire != NULL) delete _wire;
  _wire = NULL;
  _shadowSelected=NULL;
  
  _parent->updateGL();
}

void Rotoscope::killAnyRoto() {
  _parent->setMouseTracking(0);
  _seedPoint.Set(-1,-1);
  if (_currPath) {
    _currPath->clearPrevCorrespondence();
    if (!_toTrack.empty() && _toTrack.back() == _currPath)
      _toTrack.pop_back();
  }
  //  if (_drawing) {
  if (_currPath) {
    if (_tool!=T_OVERDRAW)  // if overdrawing, path has not been added to rotocurves
      _rc->deletePath(_currPath);
    else {
      delete _currPath;
      //_spliceFlag=-1;
    }
    _currPath = NULL;
  }
  if (_wire) {
    delete _wire;
    _wire = NULL;
  }
  //}
  _corrShow = false;
  _drawing = false;
  _shadowSelected=false;
}

Rotoscope::~Rotoscope() {
  delete _graph;
  if (_wire!=NULL) delete _wire;
}



void Rotoscope::addPathSegment() { 
  assert(_currPath && _wire);
  DynArray<Vec2i,100> backwards;
  _wire->startTracePath(_freePoint);
  const Vec2i* loc;
  while ((loc = _wire->prevStep()) != NULL)
    backwards.add(*loc);

  // reverse point order and add to _currPath
  _currPath->ensureCapacity(_currPath->getNumElements() +
			    backwards.getNumElements());

  // make sure not to repeat seed point twice
  int lastElement;
  if (_currPath->getNumElements() != 0)
    lastElement = backwards.getNumElements()-1;
  else
    lastElement = backwards.getNumElements()-2;
  for (int i = lastElement; i>-1; i--) {
    Vec2i vi = backwards.getElement(i);
    _currPath->addVertex(float(vi.x())+.5, float(vi.y())+.5);
  }
  assert(_currPath->okEnds());
} 

// CORRECTNESS: need to fill in correspondences
void Rotoscope::trackForwards(RotoPath* A) {
   if (A->nextC()) return;

   /*
   KLT_FullCPyramid *aPyrm = _is->getFullCPyramid(_frame), 
     *bPyrm = _is->getFullCPyramid(_frame+1);

   
   Vec2f center, dim;
   A->centroidBbox(&center,&dim);
   Vec2f lastT = A->getLastTranslation();
   double Z[12];
   A->getTrans(Z);     Z[0] -= 1.f; Z[3] -= 1.f;
   A->getTrans(Z+6);   Z[6] -= 1.f; Z[9] -= 1.f;
   Z[10] += lastT.x(); Z[11] += lastT.y();
   _is->globalTC.colorTrack(aPyrm, bPyrm, (int)dim.x(), (int)dim.y(), center.x(), center.y(), Z); 
   printf("translation %f %f, loc %f %f\n\n\n",Z[10], Z[11], center.x()+Z[10], center.y()+Z[11]);
   
   // apply transform to patch (remember it's a delta)
   RotoPath *B = new RotoPath(*A);
   B->setLastTranslate(Z[10] - Z[4], Z[11] - Z[5]);
   B->setTrans(Z+6);
   A->setNextC(B);
   B->setPrevC(A);
   _is->getRotoCurves(_frame+1)->addPath(B);
   */

   Vec2f *pts = new Vec2f[A->getNumElements()], 
     *delta = new Vec2f[A->getNumElements()]; // should initialize to 0,0 automatically
   int i;
   for (i=0; i<A->getNumElements(); i++)
     pts[i] = A->getElement(i);
   //_is->globalTC.curveTrack(aPyrm, bPyrm, A->getNumElements(), pts, 1,  G!
   //		    A->lowHeight(), A->highHeight(), delta);

   for (i=0; i<A->getNumElements(); i++)
     pts[i] += delta[i];


   RotoPath *B = new RotoPath();
   B->setFixed(false);
   B->setLowHeight(A->lowHeight());
   B->setHighHeight(A->highHeight());
   for (i=0; i<A->getNumElements(); i++)
     B->addVertex(pts[i]);
   A->setNextC(B);
   B->setPrevC(A);
   B->buildIPrevCorrs();
   A->buildINextCorrs();
   //B->setTrackSubsampling(A->trackSubsampling());
   _is->getRotoCurves(_frame+1)->addPath(B);  
}

void Rotoscope::frameChange(int i) {

  if (_visible)
    setImageGraph(NULL);
  
  _shadowSelected = NULL;
  if (_selected.size()==1) {
    if (_propMode == 0  && (i-_frame == 1))
      trackForwards(*(_selected).begin()); 
  }

  /*if (_tool==T_OVERDRAW && _spliceFlag==-2 && _selected.size()==1) {
    finishSplicedCurve();
    rejigWrapper(false);
  }
  _spliceFlag=-1;
  */

  _transformer.clear();
  if (!_selected.empty()) {
    _selected.clear();
    emit selectChanged();
  }

  if (_tool == T_MANUAL && _currPath) {
    finishManualRotoCurve();
  }
  else if (_tool==T_TRACK && !_toTrack.empty()) {
    performTracks(_lastFrameTouched, _frame);
  }

  if (_changeLastFrameTouched>-1) {
    _lastFrameTouched = _changeLastFrameTouched;
    _changeLastFrameTouched=-1;
  }

  // execute change
  _rc = _is->getRotoCurves(i);
  _frame = i;


  if (_corrShow) {
    _corrShow = 0;
    _currPath = NULL;
  }

}

void Rotoscope::toolChange(TrackTool id) {
  //assert (id>=0 && id<2);

  if (_tool == T_MANUAL && _currPath) {
    finishManualRotoCurve();
  }

  else if (_tool==T_TRACK && id != T_OVERDRAW && !_toTrack.empty()) {
    performTracks(_lastFrameTouched, _frame);
  }
  else
    killAnyRoto();

  if (_tool == T_TRANSFORM && id != T_TRANSFORM)
    _transformer.clear();

  if (id== T_TRANSFORM && _tool!= T_TRANSFORM) {
      _selected.clear();
      emit selectChanged();
  }

  /*if (_tool==T_OVERDRAW && id!=T_OVERDRAW && _spliceFlag==-2 && _selected.size()==1) {
    finishSplicedCurve();
    rejigWrapper(false);
  }
  _spliceFlag=-1;
  
  if (id==T_OVERDRAW && _tool!=T_OVERDRAW && _selected.size()>1) {
    _selected.clear();
    emit selectChanged();
  }
  */

  _tool = id;
  if (_corrShow) {
    _corrShow = 0;
    _currPath = NULL;
  }
  
  if (id==T_DRAW)
    if (!_selected.empty()) {
      _selected.clear();
      emit selectChanged();
    }
  clearCorrForw();
  _parent->updateGL();
}


void Rotoscope::makeVisible() {
  if (_visible) return;
  _visible = true;
  //if (!_graph)
  //setImageGraph(_is->getImageGraph(_frame));
}

void Rotoscope::makeInvisible() {
  if (!_visible) return;
  _visible = false;
  delete _graph; _graph = NULL;
}


// enab lists item menus that should be DISABLED
CmdMap* Rotoscope::getPopupFunctions(const Vec2i& loc, DynArray<int,5>* enab) {
  CmdMap *cm =  new CmdMap();
  _popupPicked = _rc->pickPrimitive(loc.x(), loc.y());
  _popupLoc = unproject(loc.x(), loc.y(), _h);
  if (_popupPicked) {
    _popupCtrl = _popupPicked->pickCtrl(_popupLoc.x(), _popupLoc.y());
    cm->insert(make_pair(TRACK_RIGHT, string("Track right only")));
    cm->insert(make_pair(TRACK_LEFT, string("Track left only")));
    cm->insert(make_pair(TRACK_BOTH, string("Track both sides")));
    //cm->insert(make_pair(LONG_TRACK_RSPAN, string("Long track span")));
    cm->insert(make_pair(DELETE_PATH, string("Delete path")));
    cm->insert(make_pair(DELETE_RSPAN, string("Delete span")));
    cm->insert(make_pair(SELECTR, string("Select")));
    cm->insert(make_pair(TRI_WIDTH, string("Trimap Width")));
    cm->insert(make_pair(MAKE_DRAW_PATH, string("Make drawpath")));
    if (!(_popupPicked->nextC() && _popupPicked->nextCorrs()))
      cm->insert(make_pair(CORRFORW, string("Correspond forwards")));
    //if (_popupPicked->nextC() && _popupPicked->nextCorrs())
    //enab->add(CORRFORW);
    if (_popupPicked->trackEdges())
      cm->insert(make_pair(TRACK_EDGES, string("Don't track edges")));
    else
      cm->insert(make_pair(TRACK_EDGES, string("Track edges")));
    if (_popupPicked->fixed())
      cm->insert(make_pair(R_FIX, string("Unfix curve")));
    else
      cm->insert(make_pair(R_FIX, string("Fix curve")));
    cm->insert(make_pair(RPRINT, string("Print info")));
    if (_popupPicked->distanceToEnds2(_popupLoc) < 10) {
      cm->insert(make_pair(MAKE_JOINTS, string("Make joints")));
      cm->insert(make_pair(BREAK_JOINTS, string("Break joints")));
    }
    if (_popupPicked->distanceToInternalEnds2(_popupLoc, NULL) < 10)
      cm->insert(make_pair(SPLIT_CURVE, string("Break curve")));
    else
      cm->insert(make_pair(SPLIT_SEGMENT, string("Split segment")));
    if (_popupCtrl>=0 && _popupPicked->somehowFixed(_popupCtrl))
      cm->insert(make_pair(UNFIX_POINT, string("Unfix point")));
    //printf("dist %f\n",_popupPicked->distanceToInternalEnds2(_popupLoc, NULL));
    //cm->insert(make_pair(REJIG, string("Rejig")));
    if (_selected.size() > 1)
      cm->insert(make_pair(FORM_REGION, string("Form Region")));
  }
  else {
  
    _regionPicked = _rc->pickRegion(_popupLoc.x(), _popupLoc.y());
    if (_regionPicked) {
      printf("region picked\n");
      cm->insert(make_pair(PROPAGATE_REGION, string("Propagate Region")));
    }

  }
  return cm;
}


void Rotoscope::callPopup(const int which) {
  switch(which) {
    /*case CLEAR_ACTIVES:
    clearActives();
    break;
    */
  case TRACK_RIGHT:
    _popupPicked->setLowHeight(0);
    _popupPicked->setHighHeight(4);
    break;
    
  case TRACK_LEFT: 
    _popupPicked->setLowHeight(-4);
    _popupPicked->setHighHeight(0);
    break;
  
  case TRACK_BOTH:
    _popupPicked->setLowHeight(-4);
    _popupPicked->setHighHeight(4);
    break;

  case TRACK_EDGES:
    _popupPicked->flipTrackEdges();      
    break;

  case R_FIX:
    _popupPicked->flipFixed();
    break;

  case DELETE_PATH:
    if (!_selected.empty() && _popupPicked==*(_selected.begin()))
      if (_selected.size()>0) {
	_selected.clear();
	emit selectChanged();
      }
    if (_popupPicked==_currPath) {
      killAnyRoto();
    }
    else {
      _rc->deletePath(_popupPicked);
    }
    _parent->updateGL();
    break;
  case DELETE_RSPAN:
    deleteSpan(_popupPicked);
    _parent->updateGL();
    break;
  case LONG_TRACK_RSPAN:
    //longTrackSpan(_popupPicked);
    break;
  case SELECTR:
    if (_shift) {   // multiple selection
      if (find(_selected.begin(), _selected.end(), _popupPicked) == _selected.end())
      _selected.push_back(_popupPicked);
    }
    else {        // not multiple selection
      _selected.clear();
      _selected.push_back(_popupPicked);
    }
    emit selectChanged();
    _changeLastFrameTouched=_frame;
    _parent->updateGL();
    break;
  case CORRFORW:
    _corrForw = _popupPicked;
    _corrForwFrame = _frame;
    break;
  case RPRINT:
    _popupPicked->printInfo();
    break;
  case REJIG:
    //rejig(_popupPicked);
    break;    
  case MAKE_JOINTS:
    tryToMakeJoints();
    break;
  case BREAK_JOINTS:
    tryToBreakJoints();
    break;
  case SPLIT_CURVE:
    splitCurve();
    break;
  case SPLIT_SEGMENT:
    splitSegment();
    break;
  case FORM_REGION:
    formRegion();
    break;
  case PROPAGATE_REGION:
    propagateRegion();
    break;
  case UNFIX_POINT:
    _popupPicked->unfixCtrl(_popupCtrl);
    _parent->updateGL();
    break;
  case TRI_WIDTH:
    setTrimapWidth();
    break;
  case MAKE_DRAW_PATH:
    makeDrawPath();
    break;
  }

  _popupPicked = NULL;
  _regionPicked = NULL;
}

void Rotoscope::makeDrawPath() {
  assert(_popupPicked);
  emit rotoToDrawPath();
}

void Rotoscope::setTrimapWidth() {
  assert(_popupPicked);
  bool ok;
  int res = QInputDialog::getInteger("yo", "Width:", 5, 0, 100, 1, &ok, _parent);
  if (ok)
    _popupPicked->setTrimapWidth(res);
}

void Rotoscope::propagateRegion() {
  assert(_regionPicked);
  
  //_rc->copyForwardForever(_regionPicked);
  RotoRegion *r = _regionPicked, *c;
  int i=1;
  while (r) {
    c = r->tryToCopyForward();
    if (c) {
      _is->getRotoCurves(_frame+i)->addRotoRegion(c);
      //assert(_is->getRotoCurves(_frame+i)->numRegions() == 1); // G!
    }
    ++i;
    r = c;
  }
}

void Rotoscope::formRegion() {
  assert(_selected.size() > 1);
  RotoRegion* rr = new RotoRegion(_selected);
  _rc->addRotoRegion(rr);
}

void Rotoscope::deleteSpan(RotoPath* path) {
  PathV::iterator which = find(_selected.begin(), _selected.end(), _popupPicked);
  if (which != _selected.end()) {
    _selected.erase(which);
    emit selectChanged();
  }
  //if (_popupPicked==_selected) {
  //_selected=NULL;
  //emit selectChanged();
  //}
  if (_popupPicked==_currPath) { // G! probably a problem for manual drawing
    _currPath=NULL;
    _corrShow=false;
  }

  int a=_is->getStart(),b=_is->getEnd();
  if (!getRange(a,b))
    return;

  RotoPath *curr=path, *prev;
  int frame = _frame;
  while (curr->prevC() != NULL && frame>a) {
    curr = curr->prevC();
    frame--;
  }
  while (curr->nextC() != NULL && frame<a) {
    curr = curr->nextC();
    frame++;
  }
  
  prev = curr;
  curr = curr->nextC();
  while (curr!=NULL && frame < b) {
    _is->getRotoCurves(frame)->deletePath(prev);
    prev = curr;
    curr = curr->nextC();
    frame++;
  }
  _is->getRotoCurves(frame)->deletePath(prev);
}

/*
void Rotoscope::clearActives() {
  for (int i=_is->getStart(); i<_is->getEnd(); i++)
    _is->getRotoCurves(i)->setActive(NULL);
  
}
*/

void Rotoscope::runCan() {

 
  RotoPath* a = new RotoPath();
  a->buildCan(true);
  _is->getRelativeRotoCurves(0)->addPath(a);
  RotoPath* b = new RotoPath();
  b->buildCan(false);
  _is->getRelativeRotoCurves(6)->addPath(b);
  b->startPrevCorrespondence(a);
  b->finishPrevCorrespondence();
  //a->buildINextCorrs();
  //b->buildIPrevCorrs();
  _lastFrameTouched = 1;
  _toTrack.push_back(b);
  _parent->updateGL();
 
  /*

  _canNum++;
  RotoPath* tmpp = _is->getRotoCurves(1)->getCurve(0);
  RotoPath* rp = new RotoPath(*tmpp );
  _rc->addPath(rp);
  _currPath = rp;
  _currPath->startPrevCorrespondence(tmpp);
  finishRotoCurve();
  _parent->updateGL();
  */
  /*
  if (_canNum==0) {
    for (int i=70; i<=146; i++) {
      rp->addVertex(i,89);
    }
    _currPath = rp;
    finishRotoCurve();
  }
  else {
    _currPath = rp;
    RotoPath* tmpp = getPrevActive(&_propFrame);//_is->getRotoCurves(_frame-1); 
    if ( tmpp != NULL) 
      _currPath->startPrevCorrespondence(tmpp);
    for (int i=75; i<=151; i++) {
      rp->addVertex(i,94);
    }
    finishRotoCurve();
  }
  
  _canNum++;
  _parent->updateGL();
  */
}


bool Rotoscope::modTrackPossible() const {
  if (_selected.empty() || _tracking) return false;
  //for (set<RotoPath*>::const_iterator c = _selected.begin(); c!=_selected.end(); ++c) 
  //if (!(*c)->fixed()) return true;
  //return false;
  
  return true;
}


void Rotoscope::tryToMakeJoints() {
  assert(_popupPicked);
  //Vec2f _popupLoc = unproject(_popupLoc.x(), _popupLoc.y(), _h);
  //int whichSide;
  //RotoPath* whichPath = _rc->distanceToEnds2(_popupLoc, whichSide);
  //if (!whichPath) return;

  float d0 =  _popupPicked->distanceToSide2(_popupLoc, 0);
  float d1 = _popupPicked->distanceToSide2(_popupLoc, 1);
  if (d0 < d1 && d0 < 10) {
    _rc->setupJointsOneSide(_popupPicked, 0);
    makeJointsConsistent(_popupPicked->joints(0));
  }
  else if (d1 < 10) {
    _rc->setupJointsOneSide(_popupPicked, 1);
    makeJointsConsistent(_popupPicked->joints(1));
  }

  _rc->checkAllJoints();
}

void Rotoscope::tryToBreakJoints() {
  assert(_popupPicked);
  float d0 =  _popupPicked->distanceToSide2(_popupLoc, 0);
  float d1 = _popupPicked->distanceToSide2(_popupLoc, 1);
  if (d0 < d1 && d0 < 10) {
    _popupPicked->destroyAllCoincidentJoints(0);
  }
  else if (d1 < 10) {
    _popupPicked->destroyAllCoincidentJoints(1);
  }

  _rc->checkAllJoints();
}


void Rotoscope::splitCurve() {
  assert(_popupPicked);
  int whichCtrl=-1;
  double d = _popupPicked->distanceToInternalEnds2(_popupLoc, &whichCtrl);
  if (d>10) return;
  assert(whichCtrl != -1);

  PathV::iterator which = find(_selected.begin(), _selected.end(), _popupPicked);
  if (which != _selected.end()) {
    _selected.erase(which);
    emit selectChanged();
  }

  _rc->splitCurve(_popupPicked, whichCtrl);
}

void Rotoscope::splitSegment() {
  assert(_popupPicked && !_popupPicked->samplingDirty());
  _popupPicked->splitSegment(_popupPicked->findClosestT(_popupLoc));
  _rc->checkAllJoints();
}


void Rotoscope::setViewRegion(const bool s) { 
  _viewRegions = s; _parent->updateGL();
}

void Rotoscope::setVisEffort(const bool s) { 
  _visEffort = s; _parent->updateGL();
}


void Rotoscope::calcLerp() {
  for (int i=_is->getStart(); i<_is->getEnd(); i++)
    _is->getRotoCurves(i)->calcLerp();
}

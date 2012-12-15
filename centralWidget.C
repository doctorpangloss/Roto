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



//#include <Q3PopupMenu.h>
#include <qcursor.h>
#include <qmenu.h>

#include "centralWidget.h"

CentralWidget::CentralWidget(const QGLFormat &format, ImgSequence* is, QWidget *parent) : QGLWidget(format, parent,NULL), 
  _is(is) {
  _w = is->width();
  _h = is->height();
  setFixedSize(_w,_h);
  setMinimumSize(_w,_h);
  setMaximumSize(_w,_h);
  setFocusPolicy(QWidget::focusPolicy());
  _frame = _is->getStart();
  _imageShow = true;
  _propMode = 1;
  AbstractPath::_is = is; 
  RotoRegion::_w = _w;
  RotoRegion::_h = _h;
  //DrawPatch::_is = is;
  //Joint::_gid=0;
  //_toolNum = 0;
  _zoomMode = false;

  QImage buf;
  buf = _is->getImage(_frame);
  _glback = QGLWidget::convertToGLFormat(buf);

  _roto = new Rotoscope(_is,this);
  _draw =  new Draw(is,this);
  _module = _roto;
  _mode = 0;

  _center.Set(_w/2.f, _h/2.f);
  _zoom = 1.f;
  _iAlpha=100;
  _dAlpha=100;
  _ctrlDown=false; _altDown=false;
  _off = false;

  connect(_draw, SIGNAL(restoreCursor()), this, SIGNAL(restoreCursor()));
  connect(_roto, SIGNAL(rotoToDrawPath()), this, SLOT(rotoToDrawPath()));
}

void CentralWidget::createDlist() {
  glNewList(_backDlist, GL_COMPILE);
  glClear(GL_COLOR_BUFFER_BIT);
  glRasterPos2i(0,_h);
  glDrawPixels(_w,_h,GL_RGBA,GL_UNSIGNED_BYTE,_glback.bits());
  glEndList();
}


void CentralWidget::initializeGL() { 
  //glClearColor( 0.0, 0.0, 0.0, 0.0 ); //G!
  glClearColor( 1.0, 1.0, 1.0, 0.0 );
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glClearDepth(1.0);
  glClear(GL_DEPTH_BUFFER_BIT);
  glDisable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  glClearStencil(0x0);
  glDisable(GL_STENCIL_TEST);
  //glEnable(GL_LINE_SMOOTH); // G!
  glLineStipple(2,0xAAAA);
  //glDisable(GL_DEPTH_TEST);
  _backDlist = glGenLists(1);
  //DrawPath::_dAlpha=1.;
  //DrawPatch::_dAlpha=1.;
  //glEnable(GL_BLEND); // G!
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // SKQ 5-19-03
  //glEnable( GL_BLEND ) ; // G!

}

void CentralWidget::resizeGL(int width, int height) {
  assert(width==_w && height == _h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0f, (GLfloat) width, 0.0f, 
	  (GLfloat) height);
  glMatrixMode(GL_MODELVIEW);
  myGlLoadIdentity();
  glViewport(0, 0, width, height);
  createDlist();
}

void CentralWidget::paintGL() {
  glPushMatrix();
  zoomUpdateGL();
  if (!_module->getMinimalRenderMode()) {
    if (_imageShow && _iAlpha>0) {
      if (!_off) {
	glCallList(_backDlist);   
	//printf("error = %d\n", glGetError()); 
	assert(!glGetError());
      }
      else {	
	glClear(GL_COLOR_BUFFER_BIT);
	glPushMatrix(); glLoadIdentity(); // we'll do this in GL space
	Vec2f myCenter(_center.x(), float(_h)-_center.y());
	Vec2f fd(_w, _h); // float dimensions
	
	Vec2f lowerLeft(myCenter.x() - fd.x()/(2.f*_zoom), 
			myCenter.y() - fd.y()/(2.f*_zoom));
	glRasterPos2f(MAX( fd.x()/2.f - myCenter.x()*_zoom ,0),
		      MAX( fd.y()/2.f - myCenter.y()*_zoom ,0));
	//int newWidth = 2*_w - MAX(int(myCenter.x() + fd.x()/_zoom), _w);
	//printf("new width %d\n",newWidth);
	glPopMatrix();
	lowerLeft.Clamp(0,_w,_h);
	int offset = 4 * (int(lowerLeft.y())*_w + (int)lowerLeft.x());
	//assert(offset>=0 && offset < _w*_h); // DEBUG
	assert((int)lowerLeft.y() < _h);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, _w);
	glDrawPixels(_w-(int)ceil(lowerLeft.x()), _h-(int)ceil(lowerLeft.y()), 
		     GL_RGBA,GL_UNSIGNED_BYTE,_glback.bits() + offset);
	assert(!glGetError());
      }
    }
    else
      glClear(GL_COLOR_BUFFER_BIT);
  }

  _module->paintGL();
  glPopMatrix();
  //glFlush();
}


void CentralWidget::resetZoom() {
  _zoom = 1.f;
  _draw->reportZoom(_zoom);
  glPixelZoom(_zoom, _zoom);
  _center.Set(_w/2.f, _h/2.f); //(363,118.5);//(_w/2.f, _h/2.f);
  _off = false; 
  myGlLoadIdentity();
  updateGL();
}


void CentralWidget::mousePressEvent( QMouseEvent *e ) {
  //printf("mouse press event %d\n", e->button());
  if (e->button() == Qt::RightButton || e->button() == Qt::MidButton) {
    Vec2i loc(e->x(), e->y());
    QMenu pm;

    DynArray<int,5> enab;
    CmdMap* commands = _module->getPopupFunctions(loc, &enab);
    for (CmdMap::const_iterator c = commands->begin(); 
	 c != commands->end(); ++c) {
      pm.insertItem(QString(c->second.c_str()), c->first);
      if (enab.findElement(c->first)!=-1)   // speed, should keep enab sorted
	pm.setItemEnabled(c->first,false);
      else
	pm.setItemEnabled(c->first,true);
    }
    delete commands;

    int res = pm.popup(QCursor::pos());
    if (_module == _roto && (res==FORM_REGION || res==PROPAGATE_REGION))
      resetZoom();
    _module->callPopup(res);
  }


  else if (e->button() & Qt::LeftButton) {
    if (_ctrlDown && !_altDown) {
      _zoomStart.Set(e->x(),e->y());
      _off = true;
      
    }
    else if (_zoomMode) {
      
      //_zoomStart.Set(e->x(),e->y());
      float zoomLast = _zoom;
      if (_altDown)//((e->state() & AltButton) || (e->state() & ShiftButton))
	_zoom -= 1.f;
      else
	_zoom +=1.f;
      if (_zoom <=1.f)
	resetZoom();
      else {
	Vec2f lowerLeft(_center.x() - float(_w)/(2.f*zoomLast), _center.y()-float(_h)/(2.f*zoomLast));
	_center.Set(lowerLeft.x() + float(e->x())/zoomLast, lowerLeft.y() + float(e->y())/zoomLast);
	_off = true;
	_draw->reportZoom(_zoom);
	glPixelZoom(_zoom, _zoom);
	zoomUpdateGL();
	updateGL();
      }
    }
    else
      _module->mousePressEvent(e);
  }
}

void CentralWidget::mouseReleaseEvent( QMouseEvent *e ) {
  if (_zoomMode) { 

  }
  else
    _module->mouseReleaseEvent(e);
}

void CentralWidget::zoomUpdateGL() {
  if (_off) {
    myGlLoadIdentity();
    glTranslatef(_w/2., _h/2.,0);
    glScalef(_zoom, _zoom, 1.);
    glTranslatef(-_center.x(), -_center.y(),0);
    //printf("center: %f %f\n", _center.x(), _center.y());
  }
}

void CentralWidget::mouseMoveEvent( QMouseEvent *e ) {
  if (((_ctrlDown && !_altDown) || _zoomMode) && 
      (e->state() & Qt::LeftButton)) { // Zoom tool
    Vec2i newLoc(e->x(), e->y());
    if (_ctrlDown && !_altDown) {
      Vec2f delta(_zoomStart.x()-newLoc.x(), _zoomStart.y()-newLoc.y());
      delta /= _zoom;
      _center+=delta;
      _zoomStart=newLoc;
      
    }
    else {
      /*
      _zoom = _zoomLast + sqrt(newLoc.distanceTo2(_zoomStart)) / 100.f;
      _draw->reportZoom(_zoom);
      glPixelZoom(_zoom, _zoom);
      */
      //printf("Zoom is %f\n",_zoom);
    }
    
    zoomUpdateGL();
    QPoint curs = mapFromGlobal(QCursor::pos());
    if ((e->x() - curs.x())*(e->x() - curs.x()) +
	(e->y() - curs.y())*(e->y() - curs.y()) < 25)
      updateGL();
  }
  else
    _module->mouseMoveEvent(e);
}
void CentralWidget::keyPressEvent ( QKeyEvent * e ) {
  if (e->key() == Qt::Key_Control) {
    _ctrlDown = true; 
    _module->controlChanged(true);
    //setCursor(Qt::PointingHandCursor);
    emit restoreCursor();
  }
  if (e->key() == Qt::Key_Alt || e->key() == Qt::Key_Shift) {
    _altDown = true;
    _module->shiftChanged(true);
    emit restoreCursor();
  }
  else if (e->key() == Qt::Key_Up) {
    emit internalFrameChange(_frame+1);
    //frameChange(_frame+1);
  }
  else if (e->key() == Qt::Key_Down) {
    emit internalFrameChange(_frame-1);
    //frameChange(_frame-1);
  }
  else  
    _module->keyPressEvent(e);
}

void CentralWidget::keyReleaseEvent ( QKeyEvent * e ) {
  if (e->key() == Qt::Key_Control) {
    _ctrlDown = false; 
    _module->controlChanged(false);
    //setCursor(Qt::ArrowCursor);
    emit restoreCursor();
  }
  else if (e->key() == Qt::Key_Alt || e->key() == Qt::Key_Shift) {
    _altDown = false; printf("alt up\n");
    _module->shiftChanged(false);
    emit restoreCursor();
  }
  else  
    e->ignore();
}


void CentralWidget::tabletEvent( QTabletEvent *e) {
  //if (e->device() == QTabletEvent::Stylus)
  //printf("stylus: ");
  //else
  //printf("unknown: ");
  int k = e->pressure();
  //if ((k >=10) {
  _module->reportPressure(k);
    //printf("pressure: %d\n",e->pressure());
    //}
  e->accept();
  
}


void CentralWidget::frameChange(int i) { 
  _frame = i;
  _roto->frameChange(i);
  _draw->frameChange(i);
  QImage buf;
  buf = _is->getImage(_frame);
  _glback = QGLWidget::convertToGLFormat(buf);
  createDlist();
  updateGL(); 
}

void CentralWidget::showChanged(int s) { 
  _roto->toggleShow(s); 
  _draw->toggleShow(s);
  updateGL();
}

void CentralWidget::imageShowChanged(int s) {
  if (s==2)
    _imageShow = true;
  else
    _imageShow = false;
  updateGL();
}

void CentralWidget::propModeChanged(int s) {
  assert(s==0 || s==1);
  _propMode = s;
  _draw->setPropMode(s);
  _roto->setPropMode(s);
}

void CentralWidget::renderModeChanged(int s) {
  if (s==1)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  _draw->setRenderMode(s);
  updateGL();
}


void CentralWidget::setModule(const int i) {
  if (i!=_mode) {
    _mode = i;
    if (_mode==0) {
      _roto->makeVisible();
      _module = _roto;
      if (_roto->getTool()!=T_ZOOM)
	_zoomMode=false;
      emit rotoNow();
    }
    else {
      _roto->killAnyRoto();
      _roto->makeInvisible();
      _module = _draw;
      if (_draw->getTool()!=D_ZOOM)
	_zoomMode=false;
      emit drawNow();
    }
   
    //_roto->toolChange(0);
    //_draw->toolChange(0);
  }
  updateGL();
}


void CentralWidget::drawToolChange(DrawTool id) {
  if (id==D_ZOOM) 
    _zoomMode = true;
  else
    _zoomMode = false;
  
  _draw->toolChange(id);
}

void CentralWidget::trackToolChange(TrackTool id) {
  if (id==T_ZOOM) 
    _zoomMode = true;
  else
    _zoomMode = false;
  
  _roto->toolChange(id);
}


void CentralWidget::saveImages(const char* fileRoot) {
  //QPixmap pmp = renderPixmap(0,0,true);
  //pmp.save(fileRoot,"PNG");
  char name[200];
  for (int i = _frame; i<=_is->getEnd(); i++) {
    sprintf(name,"%s%.3d.png", fileRoot,i);
    frameChange(i);
    QPixmap pmp = renderPixmap(0,0,true);
    pmp.save(name,"PNG");
  }
  
}


void CentralWidget::myGlLoadIdentity() {
  glLoadIdentity();
  glTranslatef(0,_h,0);
  glScalef(1.,-1.,1.);
}


void CentralWidget::imageAlphaChanged(int val) { 
  _iAlpha=val; 
  //printf("image alpha %d\n",_iAlpha); 
  glPixelTransferf(GL_ALPHA_SCALE, float(_iAlpha)/100.f);
  if ((_iAlpha==0 || _iAlpha==100) && (_dAlpha==0 || _dAlpha==100))
    glDisable(GL_BLEND);
  else
    glEnable(GL_BLEND);
  updateGL(); 
}


void CentralWidget::drawAlphaChanged(int val) {
  _dAlpha=val;
  float floatVal = float(val)/100.f;
  _draw->setdAlpha(floatVal);
  DrawPath::_dAlpha = floatVal;
  //DrawPatch::_dAlpha = floatVal;
  if ((_iAlpha==0 || _iAlpha==100) && (_dAlpha==0 || _dAlpha==100))
    glDisable(GL_BLEND);
  else
    glEnable(GL_BLEND);
  updateGL();
}

void CentralWidget::rotoToDrawPath() {
  _draw->makeDrawPath(_roto->getPopupPicked());
}


void CentralWidget::runCan() {
  //_roto->runCan();  
  _is->calculateStatistics();
}

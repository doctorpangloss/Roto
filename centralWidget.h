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



#ifndef CENTRALWIDGET_H
#define CENTRALWIDGET_H

#include <qpixmap.h>
#include <qwidget.h>
#include <qevent.h>
#include <qpainter.h>
#include <qimage.h>
#include <qpoint.h>
#include <qrect.h>
#include <stdio.h>
#include <qgl.h>
#include "Texture.h"
#include "imgSequence.h"
#include "interModule.h"
#include "rotoscope.h"
#include "draw.h"

class CentralWidget : public QGLWidget {
  
  Q_OBJECT

public:

  CentralWidget(const QGLFormat &format, ImgSequence* is, QWidget *parent=0);

  void setModule(const int i); // 0 for roto, 1 for draw

  ImageGraph* getImageGraph(int i) {
    return _is->getImageGraph(i);
  }

  void drawToolChange(DrawTool id);
  void trackToolChange(TrackTool id);

  void setColor(const Vec3f& col) { _draw->setColor(col); }
  QRgb getCurrColor() { return _draw->getCurrColor(); }

  void setTexture( int i ) { _draw->setCurrTexture(i); }  

  void saveImages(const char* fileRoot);

  int getWidth() const { return _w;}
  int getHeight() const { return _h; }

  void runCan();
  
  bool drawSelected() { return (_draw->selected() != NULL); }
  int rotoSelectedNum() { return _roto->selectedSize(); }
  Draw* getDraw() { return _draw; }
  Rotoscope* getRoto() { return _roto;}
  
  bool modTrackPossible() const { return _roto->modTrackPossible(); }
  void redoTrack(bool startOver) { assert(_module==_roto); _roto->rejigWrapper(startOver); }
  void copySplinesAcrossTime() { assert(_module==_roto); _roto->copySplinesAcrossTime(); }

  //int getTool() const { return _toolNum; }

  public slots:
    
  void frameChange(int i);
  void showChanged(int s);
  void imageShowChanged(int s);
  void propModeChanged(int s); // 0 for forward, 1 for keyframe
  void renderModeChanged(int s); // 0 for normal, 1 for wireframe, 2 for points

  void calcLerp() { _roto->calcLerp(); updateGL(); }
  void cwToggle(bool s) { _roto->setVisEffort(s); }        //_is->globalTC.D2mode = s; printf("D2mode set to %d\n",s);}
  void plToggle(bool s) { _roto->setViewRegion(s); }
  void usToggle(bool s) { _is->globalTC.useImage = s; printf("useS set to %d\n",s); }
  void unToggle(bool s) { _is->globalTC.useNormalEq = s; printf("useN set to %d\n",s); }
  void psToggle(bool s) { _is->globalTC.printSingulars = s; printf("printSing set to %d\n",s); }
  void upsToggle(bool s) { _is->globalTC.usePseudo = s; printf("usePseudo set to %d\n",s); }
  void dwToggle(bool s) { _is->globalTC.dumpWindows = s; printf("dumpWindows set to %d\n",s); }
  void prToggle(bool s) { _draw->preRegister = s; printf("preRegister set to %d\n",s); }
  void stpToggle(bool s) { _roto->_showTrackPoints = s;  updateGL(); }
  void fixToggle(bool s) { _roto->_nudgeFreezesCurve = s; }
  void fhdToggle(bool s) { _draw->setFreeDraw(s); }
  void corrPropAll() { _draw->corrPropAll(); }
  
  void rotoToDrawPath();


  // texture slots defined. 
  // later this function will take a parameter but now we will have a constant 
  // texture file name. SKQ 
  //void setTexFName( void ) { _draw->setCurrTexFName("myTexture.bmp");  } 
  //void useTex ( bool u ) { _draw->setUseTex( u ) ; printf("called central widget useTex\n");  } 
  
  void imageAlphaChanged(int val);
  void drawAlphaChanged(int val);
  void resetZoom();
  void sensChanged(int val) { _draw->sensChanged(float(val)); }

  void copy() { if (_module==_draw) _draw->copy(); }
  void paste() { if (_module==_draw) _draw->paste(); }  

  void wireValChanged() {
    updateGL();
  }

  int zoomOutIn() const { // 0 translate, 1 zoom out, 2 zoom in
    if (_ctrlDown && !_altDown) return 0;
    else if (_altDown) return 1; 
    else return 2; 
  }

 signals:
  
  void internalFrameChange(int i);
  void rotoNow();
  void drawNow();
  void restoreCursor();

 protected:
  void initializeGL();

  void resizeGL(int width, int height);

  void paintGL();


  void mousePressEvent( QMouseEvent *e );
  void mouseReleaseEvent( QMouseEvent *e );
  void mouseMoveEvent( QMouseEvent *e );
  void keyPressEvent ( QKeyEvent * e );
  void keyReleaseEvent ( QKeyEvent * e );
  void tabletEvent( QTabletEvent *e);


 private:

  void createDlist();
  void zoomUpdateGL();
  void myGlLoadIdentity();

  int _frame;
  InterModule* _module;
  Rotoscope* _roto;
  Draw* _draw;  
  int _w, _h;
  ImgSequence *_is;
  int _mode; // 0 for roto, 1 for draw
  QImage _glback;
  bool _imageShow;
  int _backDlist;
  bool _minimalRenderMode;
  int _propMode;
  Vec2f _center;
  float _zoom;//, _zoomLast;
  //int _toolNum;
  bool _zoomMode;
  Vec2i _zoomStart;
  int _iAlpha, _dAlpha;
  bool _ctrlDown, _altDown, _off;  // control is for translating, alt/shift for mag glass inversion
};

#endif

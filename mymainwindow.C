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



#include "mymainwindow.h"

#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qcolordialog.h>
#include <qapplication.h>
#include <qspinbox.h>
#include <qinputdialog.h>
#include <qlayout.h>
#include <qhbox.h>
#include <qvbox.h>
#include <qlabel.h>
#include <qfiledialog.h>
#include <qstring.h>
#include <qtooltip.h>
#include <qbitmap.h>
#include "pixmaps.h"



MyMainWindow::MyMainWindow(ImgSequence* is) {
  _is = is;
  setFocusPolicy(QWidget::StrongFocus);

  QGLFormat f;
  //f.setAlpha(true);
  f.setStencil(true);
  _cw = new CentralWidget(f,is,this);
  printf("%d %d\n",_cw->format().alpha(),  _cw->format().stencil());
  setMinimumSize(_cw->width()+200, _cw->height()+60);
  //assert(_cw->format().alpha());
  assert(_cw->format().stencil());
  setCentralWidget(_cw);
  _mode = 0;
  AbstractPath::_nudgeRadius = 15;

  

  file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file,0 );  
  file->insertItem( "Quit", qApp, SLOT( quit() ), CTRL + Key_Q );
  file->insertItem("Save Curves",this,SLOT( saveCurves() ), CTRL + Key_S);
  file->insertItem("Save Curves As...",this,SLOT( saveCurvesAs() ));
  file->insertItem("Save Images",this,SLOT( saveImages() ));
  file->insertItem("Save Roto XML",this,SLOT(saveXML()));
  file->insertItem("Get statistics", this, SLOT(runCan()));
  file->insertItem("Save mattes", this, SLOT(saveMattes()));
  file->insertItem("Show lerp", _cw, SLOT(calcLerp()));
  /*
  tool = new QPopupMenu( this );
  tool->setCheckable(true);
  menuBar()->insertItem("&Tool", tool, 1);
  buildToolMenu();
  connect(tool,SIGNAL(activated(int)), this, SLOT(toolChanged(int)));
  */

  edit = new QPopupMenu(this);
  menuBar()->insertItem("Edit", edit, 2);
  edit->insertItem("Copy", _cw, SLOT(copy()), CTRL + Key_C);
  edit->insertItem("Paste", _cw, SLOT(paste()), CTRL + Key_V);


  tracking = new QPopupMenu(this);
  menuBar()->insertItem("Tracking", tracking, 3);
  tracking->insertItem("Load/Calc Tracking Data...",this,SLOT(loadCalc()), CTRL + Key_L,0);
  tracking->insertItem("Write Tracking Data...",this,SLOT(writeTrack()),CTRL + Key_D,1);
  tracking->insertItem("Copy spline through time...",this,SLOT(copyTrack()), CTRL + Key_H, 2);
  tracking->insertItem("Improve track", this, SLOT(improveTrack()),CTRL + Key_I,3);
  tracking->insertItem("Redo track", this, SLOT(redoTrack()),CTRL + Key_W, 4);
  //tracking->insertItem("Track forward",this,SLOT(trackForward()), CTRL + Key_T, -1, 5);

  drawing = new QPopupMenu(this);
  menuBar()->insertItem("Drawing", drawing, 4);
  drawing->insertItem("Correspond & Propagate All Free", _cw, SLOT(corrPropAll()));
  drawing->insertItem("Delete all drawings across time", this, SLOT(deleteAllTime()));

  /* new style pop up menu - lets user select the drawing style */
  /*style = new QPopupMenu(this); 
  menuBar()->insertItem("Style", style, 5);  

  QRgb black =  qRgb(0, 0, 0); 
  QPixmap pen_pixmap( 22, 22 ); 
  pen_pixmap.fill(black);

  style->insertItem(QIconSet(pen_pixmap), "Pen", this, SLOT(pickPenTexture()) ); 
  style->insertItem(QIconSet(QPixmap(charcoal_xpm)), "Charcoal", this, SLOT(pickCharcoalTexture()) ); 
  style->insertItem(QIconSet(QPixmap(pastel_xpm)),"Pastel", this, SLOT(pickPastelTexture()) ) ;  
  style->insertItem(QIconSet(QPixmap(pencil_texture_xpm)), "Pencil", this, SLOT(pickPencilTexture()) );  
  style->insertItem(QIconSet(QPixmap(pencil_texture_xpm)), "Pencil Stroke", this, SLOT(pickPencilStrokeTexture())); 
  style->insertItem(QIconSet(QPixmap(charcoal_xpm)), "Charcoal Stroke", this, SLOT(pickCharcoalStrokeTexture()) ); 
  style->insertItem(QIconSet(QPixmap(pastel_xpm)),"Pastel Stroke", this, SLOT(pickPastelStrokeTexture()) ) ;  
  */

  //  connect( style, SIGNAL(activated(int)), this, SLOT(pickTexture(int)) ); 
    /*  

  connect(_show,SIGNAL(stateChanged(int)), _cw,SLOT(showChanged(int)));

  col.Set(float( 0.0, 0.0, 0.0)) ; 
, 
	  float(qcol.green())/255., 
	  float(qcol.blue())/255.);

  QRgb rgb = _cw->getCurrColor();
  QPixmap newLook(22,22);

  _draw = new QToolButton(tbar);
  _draw->setIconSet(QIconSet(QPixmap(mode_xpm)));
  _draw->setToggleButton(1);
  _draw->setTextLabel("Draw/Track mode", true);
  connect(_draw,SIGNAL(stateChanged(int)), this, SLOT(drawChanged(int)));
  */


  tbar = new QToolBar (this, "top");
  extWin = new QDockWindow(QDockWindow::InDock, this);
  extWin->setOrientation(Qt::Vertical);
  extWin->setFixedExtentWidth(200);
  moveDockWindow(extWin,Qt::DockRight);
  QSpinBox* spbox = new QSpinBox(_is->getStart(), _is->getEnd(), 1, tbar);
  QToolTip::add(spbox, "Frame number");

  connect(spbox, SIGNAL(valueChanged(int)), this, SLOT(frameChange(int)));
  connect(_cw,SIGNAL(internalFrameChange(int)), spbox, SLOT(setValue(int)));

  //_roto = new QToolButton(QIconSet(QPixmap("pencil.xpm")),"Roto","Roto",
  //this,SLOT(rotoPressed()),tbar);
//_roto->setToggleButton(1);
  _draw = new QToolButton(tbar);
  _draw->setIconSet(QIconSet(QPixmap(mode_xpm)));
  _draw->setToggleButton(1);
  _draw->setTextLabel("Draw/Track mode", true);
  connect(_draw,SIGNAL(stateChanged(int)), this, SLOT(drawChanged(int)));

  _image = new QToolButton(tbar);
  _image->setIconSet(QIconSet(QPixmap(monalisa_xpm)));
  _image->setToggleButton(1);
  _image->setOn(1);
  _image->setTextLabel("Show image", true);
  connect(_image,SIGNAL(stateChanged(int)), _cw,SLOT(imageShowChanged(int)));

  _reset = new QPushButton("Reset Zoom", tbar);
  connect(_reset, SIGNAL(clicked()), _cw, SLOT( resetZoom()));

  _freehand = new QCheckBox("Freehand", extWin);
  extWin->boxLayout()->addWidget(_freehand);
  _freehand->setChecked(true);
  connect(_freehand, SIGNAL(toggled(bool)), _cw, SLOT(fhdToggle(bool)));


  _show = new QCheckBox("Show Correspondence", extWin);
  extWin->boxLayout()->addWidget(_show);
  _show->setChecked(false);
  connect(_show,SIGNAL(stateChanged(int)), _cw,SLOT(showChanged(int)));

  /*propagateMode = new QListBox(extWin);
  propagateMode->setFocusPolicy(QWidget::NoFocus);
  extWin->boxLayout()->addWidget(propagateMode);
  propagateMode->insertItem("Forward",0);
  propagateMode->insertItem("KeyFrame",1);
  propagateMode->setCurrentItem(1);
  connect(propagateMode, SIGNAL(highlighted(int)), _cw, SLOT(propModeChanged(int)));

  renderMode = new QListBox(extWin);
  renderMode->setFocusPolicy(QWidget::NoFocus);
  extWin->boxLayout()->addWidget(renderMode);
  renderMode->insertItem("Normal",0);
  renderMode->insertItem("Wireframe",1);
  renderMode->insertItem("Points",2);
  renderMode->setCurrentItem(0);
  connect(renderMode, SIGNAL(highlighted(int)), _cw, SLOT(renderModeChanged(int)));*/

  QRgb rgb = _cw->getCurrColor();
  QPixmap newLook(22,22);
  newLook.fill(rgb);
  _color = new QToolButton(QIconSet(newLook), "Current color",  // dummy xpm, will make one later
			   "Current color", this, SLOT(pickColor()), tbar);
  

  

  constantWindow = new QRadioButton("Visualize effort",extWin);
  constantWindow->setChecked(true);
  extWin->boxLayout()->addWidget(constantWindow);
  connect(constantWindow,SIGNAL(toggled(bool)),_cw, SLOT(cwToggle(bool)));
  
  /*pinLast = new QRadioButton("View regions",extWin);
  pinLast->setChecked(true);
  extWin->boxLayout()->addWidget(pinLast);
  connect(pinLast,SIGNAL(toggled(bool)),_cw, SLOT(plToggle(bool)));*/

  /*
  useNormalEq = new QRadioButton("Use normal equations",extWin);
  useNormalEq->setChecked(true);
  extWin->boxLayout()->addWidget(useNormalEq);
  connect(useNormalEq,SIGNAL(toggled(bool)),_cw, SLOT(unToggle(bool)));
  */

  useSmooth = new QRadioButton("Use Image Term",extWin);
  useSmooth->setChecked(true);
  extWin->boxLayout()->addWidget(useSmooth);
  connect(useSmooth,SIGNAL(toggled(bool)),_cw, SLOT(usToggle(bool)));

  /*
  printSingulars = new QRadioButton("Print Singulars",extWin);
  printSingulars->setChecked(false);
  extWin->boxLayout()->addWidget(printSingulars);
  connect(printSingulars,SIGNAL(toggled(bool)),_cw, SLOT(psToggle(bool)));

  usePseudo = new QRadioButton("Use Pseudoinverse",extWin);
  usePseudo->setChecked(false);
  extWin->boxLayout()->addWidget(usePseudo);
  connect(usePseudo,SIGNAL(toggled(bool)),_cw, SLOT(upsToggle(bool)));

  dumpWindows = new QRadioButton("Dump windows",extWin);
  dumpWindows->setChecked(false);
  extWin->boxLayout()->addWidget(dumpWindows);
  connect(dumpWindows,SIGNAL(toggled(bool)),_cw, SLOT(dwToggle(bool)));

  preRegister = new QRadioButton("Preregister",extWin);
  preRegister->setChecked(false);
  extWin->boxLayout()->addWidget(preRegister);
  connect(preRegister,SIGNAL(toggled(bool)),_cw, SLOT(prToggle(bool)));
  */

  showTrackPoints = new QRadioButton("Show Tracking samples\n", extWin);
  showTrackPoints->setChecked(false);
  extWin->boxLayout()->addWidget(showTrackPoints);
  connect(showTrackPoints,SIGNAL(toggled(bool)),_cw, SLOT(stpToggle(bool)));

  /*
  fixNudge = new QRadioButton("Nudge freezes whole curve\n", extWin);
  fixNudge->setChecked(true);
  extWin->boxLayout()->addWidget(fixNudge);
  connect(fixNudge,SIGNAL(toggled(bool)),_cw, SLOT(fixToggle(bool)));*/

  QHBox *nlbox = new QHBox(extWin);
  extWin->boxLayout()->addWidget(nlbox);
  QLabel *nltext = new QLabel( "Pyramid Levels:",nlbox);
  nlevels = new QLineEdit(nlbox);
  //nlevels->resize(15,nlevels->sizeHint().height());
  nlevels->setMaxLength(1);
  QString num;
  num.setNum(_is->globalTC.nPyramidLevels);
  nlevels->setText(num);
  connect(nlevels,SIGNAL(returnPressed()), this, SLOT(nlevelsChanged()));

  
  QHBox *shbox = new QHBox(extWin);
  extWin->boxLayout()->addWidget(shbox);
  QLabel *shtext = new QLabel( "Draw Curve Smoothness:",shbox);
  smoothShape = new QLineEdit(shbox);
  num.setNum(_is->globalTC.shape2Deriv);
  smoothShape->setText(num);
  connect(smoothShape,SIGNAL(returnPressed()), this, SLOT(shape2DerivChanged()));

  QHBox *sabox = new QHBox(extWin);
  extWin->boxLayout()->addWidget(sabox);
  QLabel *satext = new QLabel( "Curve D0 Smoothness:",sabox);
  smoothAlpha = new QLineEdit(sabox);
  num.setNum(_is->globalTC.smooth0Deriv);
  smoothAlpha->setText(num);
  connect(smoothAlpha,SIGNAL(returnPressed()), this, SLOT(smooth0DerivChanged()));

  QHBox *s1box = new QHBox(extWin);
  extWin->boxLayout()->addWidget(s1box);
  QLabel *s1text = new QLabel( "Curve D1 Smoothness:",s1box);
  smooth1Deriv = new QLineEdit(s1box);
  num.setNum(_is->globalTC.smooth1Deriv);
  smooth1Deriv->setText(num);
  connect(smooth1Deriv,SIGNAL(returnPressed()), this, SLOT(smooth1DerivChanged()));

  QHBox *s2box = new QHBox(extWin);
  extWin->boxLayout()->addWidget(s2box);
  QLabel *s2text = new QLabel( "Curve D2 Smoothness:",s2box);
  smooth2Deriv = new QLineEdit(s2box);
  num.setNum(_is->globalTC.smooth2Deriv);
  smooth2Deriv->setText(num);
  connect(smooth2Deriv,SIGNAL(returnPressed()), this, SLOT(smooth2DerivChanged()));

  QHBox *edgebox = new QHBox(extWin);
  extWin->boxLayout()->addWidget(edgebox);
  QLabel *edgetext = new QLabel( "Edge weight:",edgebox);
  edgeWeight = new QLineEdit(edgebox);
  num.setNum(_is->globalTC.edgeWeight);
  edgeWeight->setText(num);
  connect(edgeWeight,SIGNAL(returnPressed()), this, SLOT(edgeWeightChanged()));
  

  QVBox *slider1box = new QVBox(extWin);
  extWin->boxLayout()->addWidget(slider1box);
  QLabel *slider1text = new QLabel( "Image visibility (%):", slider1box);
  slider1 = new QSlider(0,100,50,100,Qt::Horizontal, slider1box);
  connect(slider1, SIGNAL(valueChanged(int)), _cw, SLOT(imageAlphaChanged(int)));
  
  QVBox *slider2box = new QVBox(extWin);
  extWin->boxLayout()->addWidget(slider2box);
  QLabel *slider2text = new QLabel( "Drawing visibility (%):", slider2box);
  slider2 = new QSlider(0,100,50,100,Qt::Horizontal, slider2box);
  connect(slider2, SIGNAL(valueChanged(int)), _cw, SLOT(drawAlphaChanged(int)));

  /*QVBox *ssliderbox = new QVBox(extWin);
  extWin->boxLayout()->addWidget(ssliderbox);
  QLabel *sslidertext = new QLabel( "Tablet sensitivity (%):", ssliderbox);
  sensSlider = new QSlider(0,100,50,30,Qt::Horizontal, ssliderbox);
  connect(sensSlider, SIGNAL(valueChanged(int)), _cw, SLOT(sensChanged(int)));*/
  
  _draw->setOn(0);
  rotoNow();
  connect(_cw,SIGNAL(rotoNow()), this, SLOT(rotoNow()));
  connect(_cw,SIGNAL(drawNow()), this, SLOT(drawNow()));

  _drawToolGroup = new QActionGroup(this, "draw_tools", TRUE);
  _d_draw = new QAction("Draw", QIconSet(QPixmap(paintbrush_xpm)), "Draw", Key_D, _drawToolGroup, 0, TRUE);
  _d_select = new QAction("Select", QIconSet(QPixmap(arrow_xpm)), "Select", Key_S, _drawToolGroup, 0, TRUE);
  _d_overdraw = new QAction("Overdraw", QIconSet(QPixmap(overdraw_xpm)), "Overdraw", Key_O, _drawToolGroup, 0, TRUE);
  _d_nudge = new QAction("Nudge", QIconSet(QPixmap(nudge_xpm)), "Nudge", Key_N, _drawToolGroup, 0, TRUE);
  _d_zoom = new QAction("Zoom", QIconSet(QPixmap(zoom_xpm)), "Zoom", Key_Z, _drawToolGroup, 0, TRUE);
  _d_dropper = new QAction("Dropper", QIconSet(QPixmap(dropper_xpm)), "Dropper", Key_P, _drawToolGroup, 0, TRUE);
    
  _trackToolGroup = new QActionGroup(this, "track_tools", TRUE);
  //_t_draw = new QAction("Draw", QIconSet(QPixmap(scissor_xpm)), "Draw", Key_D, _trackToolGroup, 0, TRUE);
  _t_manual = new QAction("Draw", QIconSet(QPixmap(scissor_xpm)), "Draw", Key_M, _trackToolGroup, 0, TRUE);
  //_t_track = new QAction("Track", QIconSet(QPixmap(scissor2_xpm)), "Track", Key_T, _trackToolGroup, 0, TRUE);
  _t_select = new QAction("Select", QIconSet(QPixmap(arrow_xpm)), "Select", Key_S, _trackToolGroup, 0, TRUE);
  _t_transform = new QAction("Transform", QIconSet(QPixmap(overdraw_xpm)), "Transform", Key_R, _trackToolGroup, 0, TRUE);
  //_t_overdraw = new QAction("Overdraw", QIconSet(QPixmap(overdraw_xpm)), "Overdraw", Key_O, _trackToolGroup, 0, TRUE);
  _t_nudge = new QAction("Nudge", QIconSet(QPixmap(nudge_xpm)), "Nudge", Key_N, _trackToolGroup, 0, TRUE);
  _t_zoom = new QAction("Zoom", QIconSet(QPixmap(zoom_xpm)), "Zoom", Key_Z, _trackToolGroup, 0, TRUE);

  _zoomoutCurs = new QCursor(QBitmap(32,32,myzoomout_bits, true), QBitmap(32,32,myzoomoutmask_bits, true));
  _zoominCurs = new QCursor(QBitmap(32,32,myzoomin_bits, true), QBitmap(32,32,myzoominmask_bits, true));
  _dropperCurs = new QCursor(QBitmap(32,32,dropper_bits, true), QBitmap(32,32,droppermask_bits, true), 8, 23);
  _paintCurs = new QCursor(QBitmap(32,32,paintbrush_bits, true), QBitmap(32,32,paintbrushmask_bits, true), 13, 18);
  _nudgeCurs = new QCursor(QBitmap(32,32,nudge_bits, true), QBitmap(32,32,nudgemask_bits, true));
  _transformCurs = new QCursor(QBitmap(32,32,overdraw_bits, true), QBitmap(32,32,overdrawmask_bits, true));
  _overdrawCurs = new QCursor(QBitmap(32,32,overdraw_bits, true), QBitmap(32,32,overdrawmask_bits, true));
  _pencilCurs = new QCursor(QBitmap(32,32,pencil_bits, true), QBitmap(32,32,pencilmask_bits, true));
  
  
  _d_draw->setOn(true);
  _t_manual->setOn(true);
  tbar->addSeparator();
  _trackToolGroup->addTo(tbar);
  _drawToolGroup->setEnabled(false);
  enableTrack();
  connect(_drawToolGroup, SIGNAL(selected(QAction*)), this, SLOT(aToolChanged(QAction*)));
  connect(_trackToolGroup, SIGNAL(selected(QAction*)), this, SLOT(aToolChanged(QAction*)));  
  connect(_cw->getDraw(), SIGNAL(selectChanged()), this, SLOT(selectChanged()));
  connect(_cw->getRoto(), SIGNAL(selectChanged()), this, SLOT(selectChanged()));
  connect(_cw->getDraw(), SIGNAL(currColorChanged()), this, SLOT(currColorChanged()));
  connect(_cw, SIGNAL(restoreCursor()), this, SLOT(setCursor()));
  _cw->setFocus();

  currColorChanged();
  setCursor();
  //testSuite();
}
/*
void MyMainWindow::buildToolMenu() {
  tool->clear();
  tool->insertTearOffHandle();
  if (_mode == 1) {
    tool->insertItem("Draw",(int)0,-1);
    tool->setAccel(Key_D, 0);
    //tool->insertItem("Paint",(int)2,-1);
    //tool->setAccel(Key_P, 2);
    tool->insertItem("Select",(int)1,-1);
    tool->setAccel(Key_S, 1);
    tool->insertItem("Overdraw",(int)4,-1);
    tool->setAccel(Key_O, 4);
    tool->insertItem("Nudge",(int)6,-1);
    tool->setAccel(Key_N, 6);
    tool->insertItem("Zoom",(int)3,-1);
    tool->setAccel(Key_Z, 3);
    tool->insertItem("Dropper",(int)5,-1);
    tool->setAccel(Key_P, 5);
    tool->setItemChecked(0,0==_cw->getTool());
    tool->setItemChecked(1,1==_cw->getTool());
    //tool->setItemChecked(2,2==_cw->getTool());
    tool->setItemChecked(3,3==_cw->getTool());
    tool->setItemChecked(4,4==_cw->getTool());
    tool->setItemChecked(5,5==_cw->getTool());
    tool->setItemChecked(6,6==_cw->getTool());
  }
  else {
    tool->insertItem("Draw",(int)0,-1);
    tool->setAccel(Key_D, 0);
    tool->insertItem("Track",(int)2,-1);
    tool->setAccel(Key_T, 2);
    tool->insertItem("Select",(int)1,-1);
    tool->setAccel(Key_S, 1);
    tool->insertItem("Overdraw",(int)4,-1);
    tool->setAccel(Key_O, 4);
    tool->insertItem("Nudge",(int)6,-1);
    tool->setAccel(Key_N, 6);
    tool->insertItem("Zoom",(int)3,-1);
    tool->setAccel(Key_Z, 3);
    tool->setItemChecked(0,0==_cw->getTool());
    tool->setItemChecked(1,1==_cw->getTool());
    tool->setItemChecked(2,2==_cw->getTool());
    tool->setItemChecked(3,3==_cw->getTool());
    tool->setItemChecked(4,4==_cw->getTool());
    tool->setItemChecked(6,6==_cw->getTool());
  }
}
*/

void MyMainWindow::frameChange(int i) {
  _cw->setFocus();
  _cw->frameChange(i);
}

void MyMainWindow::currColorChanged() {
  QRgb rgb = _cw->getCurrColor();
  QPixmap newLook(22,22);
  newLook.fill(rgb);
  _color->setIconSet(QIconSet(newLook));
}

void MyMainWindow::pickColor() {
  QColor col = QColorDialog::getColor(_cw->getCurrColor());
  if (!col.isValid()) return;
  Vec3f myCol(float(col.red())/255.f, 
	      float(col.green())/255.f, 
	      float(col.blue())/255.f);
  _cw->setColor(myCol);
  currColorChanged();
}

/****** figure out how to combine these slots *********/
void MyMainWindow::pickPenTexture() { 
  _cw->setTexture( PEN ) ; 
} 

void MyMainWindow::pickCharcoalTexture() { 
  _cw->setTexture( CHARCOAL ) ; 
} 

void MyMainWindow::pickPastelTexture() { 
  _cw->setTexture( PASTEL ) ; 
} 

void MyMainWindow::pickPencilTexture() { 
  _cw->setTexture( PENCIL ) ; 
}

void MyMainWindow::pickPencilStrokeTexture() { 
  _cw->setTexture( PENCIL_STROKE ) ; 
} 


void MyMainWindow::pickCharcoalStrokeTexture() { 
  _cw->setTexture( CHARCOAL_STROKE ) ; 
} 


void MyMainWindow::pickPastelStrokeTexture() { 
  _cw->setTexture( PASTEL_STROKE ) ; 
} 

/*****************
void MyMainWindow::pickTexture( int i ) { 

  printf("i chosen = %d\n", i) ; 
  assert (i < NUM_TEXTURES && i >=0 ); 
  _cw->setTexture( i ) ; 
} 
**************/

void MyMainWindow::enableTrack() {
  if (_cw->rotoSelectedNum()!=1) {
    //_t_nudge->setEnabled(false);
    //_t_overdraw->setEnabled(false);
  }
  else {
    //_t_nudge->setEnabled(true);
    //_t_overdraw->setEnabled(true);
  }
  if (_cw->modTrackPossible()) {
      tracking->setItemEnabled(2, true);
      tracking->setItemEnabled(3, true);
      tracking->setItemEnabled(4, true);
  }
  else {
    tracking->setItemEnabled(2, false);
    tracking->setItemEnabled(3, false);
    tracking->setItemEnabled(4, false);
  }
}

void MyMainWindow::enableDraw() {
  if (!_cw->drawSelected()) {
    _d_nudge->setEnabled(false);
    _d_overdraw->setEnabled(false);
  }
  else {
    _d_nudge->setEnabled(true);
    _d_overdraw->setEnabled(true);
  }
}

void MyMainWindow::selectChanged() {
  if (_mode==0)
    enableTrack();
  else
    enableDraw();
}


void MyMainWindow::drawChanged(int state) {
  int newMode;
  
  printf("state = %d\n", state ) ; 

  if (state == 2) newMode = 1;
  else newMode = 0;
  if (newMode != _mode) {
    _mode = newMode;
    //buildToolMenu();
    if (newMode==0) {
      _drawToolGroup->removeFrom(tbar);
      _trackToolGroup->addTo(tbar);
      _drawToolGroup->setEnabled(false);
      _trackToolGroup->setEnabled(true);
      enableTrack();
    }
    else {
      _trackToolGroup->removeFrom(tbar);
      _drawToolGroup->addTo(tbar);
      _drawToolGroup->setEnabled(true);
      _trackToolGroup->setEnabled(false);
      enableDraw();
    }
    _cw->setModule(_mode);
  }
  else
    assert(0);
  setCursor();
}


void MyMainWindow::setCursor() {  // cursors don't look right on macs
#ifndef __APPLE__
  int zoomer = _cw->zoomOutIn();
  if (zoomer==0)
    _cw->setCursor(Qt::PointingHandCursor);
  else if (_mode==0) { // roto
    TrackTool tool = _cw->getRoto()->getTool();
    if (tool==T_ZOOM) {
      if (zoomer==1)
	_cw->setCursor(*_zoomoutCurs); 
      else
	_cw->setCursor(*_zoominCurs); 
    }    
    else if (tool==T_NUDGE)
      _cw->setCursor(*_nudgeCurs);
    else if (tool==T_OVERDRAW)
      _cw->setCursor(*_overdrawCurs);
    else if (tool==T_DRAW || tool==T_TRACK || tool==T_MANUAL)
      _cw->setCursor(Qt::CrossCursor);
    else
      _cw->setCursor(Qt::ArrowCursor);
  }
  else if (_mode==1) { // draw
    DrawTool tool = _cw->getDraw()->getTool();
    if (tool==D_ZOOM) {
     if (zoomer==1)
	_cw->setCursor(*_zoomoutCurs); 
     else
       _cw->setCursor(*_zoominCurs);
    }
    else if (tool==D_DRAW)
      _cw->setCursor(*_paintCurs);
    else if (tool==D_DROPPER)
      _cw->setCursor(*_dropperCurs);
    else if (tool==D_NUDGE)
      _cw->setCursor(*_nudgeCurs);
    else if (tool==D_OVERDRAW)
      _cw->setCursor(*_overdrawCurs);
    
    else 
      _cw->setCursor(Qt::ArrowCursor);
  }
#else
  _cw->setCursor(Qt::ArrowCursor);
#endif
  
}


void MyMainWindow::aToolChanged(QAction* act) {
  if (_mode==0) { // roto
    //if (act==_t_draw) 
    //_cw->trackToolChange(T_DRAW);
    /*else*/ if (act==_t_manual)
      _cw->trackToolChange(T_MANUAL);
    //else if (act==_t_track)
    //_cw->trackToolChange(T_TRACK);
    else if (act==_t_select)
      _cw->trackToolChange(T_SELECT);
    else if (act==_t_transform)
      _cw->trackToolChange(T_TRANSFORM);
    ///else if (act==_t_overdraw)
    //_cw->trackToolChange(T_OVERDRAW);
    else if (act==_t_nudge)
      _cw->trackToolChange(T_NUDGE);
    else if (act==_t_zoom) 
      _cw->trackToolChange(T_ZOOM);
    
    else
      assert(0);
  }

  else if (_mode==1) { // draw
    if (act==_d_draw)
      _cw->drawToolChange(D_DRAW); 
    else if (act==_d_select)
      _cw->drawToolChange(D_SELECT);
    else if (act==_d_overdraw)
      _cw->drawToolChange(D_OVERDRAW);
    else if (act==_d_nudge)
      _cw->drawToolChange(D_NUDGE);
    else if (act==_d_zoom) 
      _cw->drawToolChange(D_ZOOM);
    else if (act==_d_dropper)
      _cw->drawToolChange(D_DROPPER);
    else
      assert(0);
  }
  
  setCursor();
}

/*
void MyMainWindow::toolChanged(int id) {
  if (id<0) return;
  printf("toolChanged %d\n",id);
  for (unsigned int i=0; i<tool->count(); i++) {
    if (id==i)
      tool->setItemChecked(i,true);
    else
      tool->setItemChecked(i,false);
  }
  _cw->toolChange(id);
}
*/

/*
void MyMainWindow::testSuite() { // DEBUG
  RotoCurves *rc0 = _is->getRotoCurves(0), *rc1 = _is->getRotoCurves(1);
  DrawCurves *dc0 = _is->getDrawCurves(0);

  double ratio = 4*M_PI / 199.;
  RotoPath *newrp0 = new RotoPath(), *newrp1;
  int i;
  for (i=50; i<250; i++)
    newrp0->add(Vec2f(i,50));
  rc0->addPath(newrp0);
  newrp1 = new RotoPath();
  for (i=50; i<250; i++)
    newrp1->addVertex(Vec2f(i, 100. + 10*sin(double(i)*ratio)));
  rc1->addPath(newrp1);
  newrp1->startCorrespondence(newrp0);
  newrp1->finishCorrespondence();
  
  DrawPath *newdp = new DrawPath();
  
  for (i=50; i<250; i+=1)
    newdp->addVertex(i,(int)rint(50. + 20.*sin(double(i)*ratio)));
  
  //for (i=50; i<250; i+=1)
  //newdp->addVertex(i,25);
  dc0->addPath(newdp);
  newdp->correspondRoto(rc0);
  _cw->_draw->interpolateForwards(newdp,1);
  _cw->updateGL();
}
*/

void MyMainWindow::saveCurvesAs() {
  //bool ok = FALSE;
  //QString root = QInputDialog::getText(QString::null,"Enter file root",
  //			       QLineEdit::Normal, _is->getFileRoot(),
  //			       &ok,this); 
  QString s = QFileDialog::getExistingDirectory(NULL,
						this, 0,  "Choose a directory for the clip",
						TRUE );
  
  //if (ok && !root.isEmpty())
  printf("%s\n",s.ascii()); 
  if (!s.isEmpty()) {
    s.append(QString("/frame"));
    _is->saveCurvesAs(s);
  }
}

void MyMainWindow::saveCurves() {
  if (!_is->fileRootInited())
    saveCurvesAs();
  else
    _is->saveCurvesQt();
}


void MyMainWindow::saveImages() {
  bool ok = FALSE;
  QString root = QInputDialog::getText(QString::null,"Enter file root",
				       QLineEdit::Normal, _is->getFileRoot(),
				       &ok,this); 
  
  if (ok && !root.isEmpty())
    _cw->saveImages(root);
}

void MyMainWindow::saveMattes() {
  _is->saveMattes();
}

void MyMainWindow::saveXML() {
  //bool ok = FALSE;
  QString res = QFileDialog::getSaveFileName(QString::null, QString::null, this);
  
  if (!res.isEmpty()) {
    std::ofstream fp(res);
    _is->saveRotoXML(fp);
  }
}


void MyMainWindow::loadCalc() {
  Vec2i range;
  RangeDialog rd(this, _is->getStart(), _is->getEnd());
  if (rd.exec() == QDialog::Accepted) {
    range = rd.getRange();
    if (range.x() > range.y() || range.x() < _is->getStart() ||
	range.y() > _is->getEnd())
      return;
    printf("Loading %d to %d\n",range.x(), range.y());
    _is->loadCalcFullCPyramids(range.x(),range.y());
    _is->loadCalcFullEdgePyramids(range.x(),range.y());
  }
}

void MyMainWindow::writeTrack() {
  Vec2i range;
  // select range
  RangeDialog rd(this, _is->getStart(), _is->getEnd());
  if (rd.exec() == QDialog::Accepted) {
    range = rd.getRange();
    if (range.x() > range.y() || range.x() < _is->getStart() ||
	range.y() > _is->getEnd())
      return;
    for (int i=range.x(); i<=range.y(); i++) {
      _is->writeFullCPyramid(i);
      _is->writeFullEdgePyramid(i); 
    }
  }
}

void MyMainWindow::nlevelsChanged() {
  QString out = nlevels->text();
  bool ok;
  int res = out.toInt(&ok);
  assert(ok);
  _is->nlevelsChanged(res);
  nlevels->setText(out);
}

void MyMainWindow::smooth0DerivChanged() {
  QString out = smoothAlpha->text();
  bool ok;
  float res = out.toFloat(&ok);
  assert(ok);
  _is->smooth0DerivChanged(res);
  smoothAlpha->setText(out);
}

void MyMainWindow::smooth2DerivChanged() {
  QString out = smooth2Deriv->text();
  bool ok;
  float res = out.toFloat(&ok);
  assert(ok);
  _is->smooth2DerivChanged(res);
  smooth2Deriv->setText(out);
}


void MyMainWindow::smooth1DerivChanged() {
  QString out = smooth1Deriv->text();
  bool ok;
  float res = out.toFloat(&ok);
  assert(ok);
  _is->smooth1DerivChanged(res);
  smooth1Deriv->setText(out);
}

void MyMainWindow::edgeWeightChanged() {
  QString out = edgeWeight->text();
  bool ok;
  float res = out.toFloat(&ok);
  assert(ok);
  _is->edgeWeightChanged(res);
  edgeWeight->setText(out);
}


void MyMainWindow::shape2DerivChanged() {
  QString out = smoothShape->text();
  bool ok;
  float res = out.toFloat(&ok);
  assert(ok);
  _is->shape2DerivChanged(res);
  smoothShape->setText(out);
}


void MyMainWindow::rotoNow() {
  _freehand->setEnabled(false);
  _color->setEnabled(false);
  //sensSlider->setEnabled(false);
  menuBar()->setItemEnabled(2,false);
  menuBar()->setItemEnabled(3,true);
  menuBar()->setItemEnabled(4,false);
  //renderMode->setEnabled(false);
  //propagateMode->setEnabled(true);
  
  nlevels->setEnabled(true);
  /*
  smoothAlpha->setEnabled(true);
  smooth1Deriv->setEnabled(true);
  smooth2Deriv->setEnabled(true);
  smoothShape->setEnabled(false);*/ // G!
}

void MyMainWindow::drawNow() {
  _freehand->setEnabled(true);
  _color->setEnabled(true);
  //sensSlider->setEnabled(true);
  menuBar()->setItemEnabled(2,true);
  menuBar()->setItemEnabled(3,false);
  menuBar()->setItemEnabled(4,true);
  //renderMode->setEnabled(true);
  //propagateMode->setEnabled(true);

  nlevels->setEnabled(false);
  /*
  smoothAlpha->setEnabled(false);
  smooth1Deriv->setEnabled(false);
  smooth2Deriv->setEnabled(false);
  smoothShape->setEnabled(true); */ // G!
}


void MyMainWindow::deleteAllTime() {
  for (int i=_is->getStart(); i<_is->getEnd(); i++) {
    _is->getDrawCurves(i)->deleteAll();
    //_is->getRotoCurves(i)->assertNoDraws(); 
  }
}

void MyMainWindow::improveTrack() {
  _cw->redoTrack(false);
}


void MyMainWindow::redoTrack() {
  _cw->redoTrack(true);  
}

void MyMainWindow::copyTrack() {
  _cw->copySplinesAcrossTime();
}

void MyMainWindow::trackForward() {
  _cw->getRoto()->trackForward();
}

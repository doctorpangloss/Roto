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



#ifndef MYMAINWINDOW_H
#define MYMAINWINDOW_H

#include <qmainwindow.h>
#include <qstring.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qcheckbox.h>
#include <qradiobutton.h>
#include <qlineedit.h>
#include <qslider.h>
#include <Qt3Support/q3listbox.h>
#include <Qt3Support/q3popupmenu.h>
#include <qaction.h>
#include <qcursor.h>
#include <Qt3Support/q3dockwindow.h>
#include "imgSequence.h"
#include "centralWidget.h"
#include "draw.h"
#include "rangeDialog.h"

class Q3PopupMenu;

class MyMainWindow: public QMainWindow {
  Q_OBJECT

    public:
  
  MyMainWindow(ImgSequence* is);    
    
  //  void testSuite();

  Q3DockWindow* extWin;

  public slots:
    
  void saveCurvesAs();
  void saveCurves();
  void saveXML();
  void saveImages();
  void drawChanged(int state);
  //void toolChanged(int id);
  void pickColor();

  // void pickTexture(int i) ; 
  // fix later to combine all these into one slot: 
  void pickPenTexture() ; 
  void pickCharcoalTexture(); 
  void pickPastelTexture() ; 
  void pickPencilTexture() ; 
  void pickPencilStrokeTexture() ; 
  void pickCharcoalStrokeTexture() ; 
  void pickPastelStrokeTexture() ; 

  void loadCalc();
  void writeTrack();
  void improveTrack();
  void redoTrack();
  void trackForward();
  void copyTrack();
  void nlevelsChanged();
  void smooth0DerivChanged();
  void smooth2DerivChanged();
  void smooth1DerivChanged();
  void edgeWeightChanged();
  void shape2DerivChanged();
  void runCan() { _cw->runCan(); }
  void rotoNow();
  void drawNow();
  void deleteAllTime();
  void aToolChanged(QAction* act);
  void selectChanged();
  void currColorChanged();
  void setCursor();

  void saveMattes();

  void frameChange(int i);


 private:
  void enableDraw();
  void enableTrack();
  //void buildToolMenu();

  CentralWidget* _cw;
  ImgSequence* _is;
  int _mode; // 0 for roto, 1 for draw
  unsigned int _tool; 

  QActionGroup *_drawToolGroup, *_trackToolGroup;
  QAction *_d_draw, *_d_select, *_d_overdraw, *_d_nudge, *_d_zoom, *_d_dropper;
  QAction *_t_draw, *_t_manual, *_t_track, *_t_select, *_t_transform, *_t_overdraw, *_t_nudge, *_t_zoom;
  QCursor *_zoomoutCurs, *_zoominCurs, *_dropperCurs, *_paintCurs, *_nudgeCurs, *_transformCurs,
    *_overdrawCurs, *_pencilCurs;

  QToolButton *_draw, *_image, *_color ;
  QCheckBox  *_show, *_freehand;
  Q3PopupMenu *tool, *file, *tracking, *drawing, *edit, *style ;
  QRadioButton *constantWindow, *pinLast, *useSmooth, *useNormalEq,
    *printSingulars, *usePseudo, *preRegister, *dumpWindows, *showTrackPoints, *fixNudge;
  QLineEdit *nlevels, *smoothAlpha, *smooth2Deriv, *smooth1Deriv, *smoothShape, *edgeWeight;
  QPushButton *_reset;
  QSlider *slider1, *slider2, *sensSlider;
  QToolBar *tbar;
  Q3ListBox *propagateMode, *renderMode;

};

#endif

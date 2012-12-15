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



#ifndef DRAW_H
#define DRAW_H


#ifdef __APPLE__
using namespace std;
#endif
#include <qwidget.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <qgl.h>
#include <qimage.h>
#include <qdatetime.h>
#include <vector>
#include <string>
#include "Texture.h" 
#include "imgSequence.h"
#include "interModule.h"
//#include "drawPatch.h"
#include "KLT/klt.h"



enum DrawTool { D_DRAW, D_SELECT, D_OVERDRAW, D_NUDGE, D_ZOOM, D_DROPPER, D_PAINT };

enum DrawPopupFunction { DELETE_PATCH, LONG_TRACK_SPAN,CHANGE_COLOR, CHANGE_FILL_COLOR, UNFILL, 
			 CORRESPOND_ROTO, UNCORRESPOND, PROPAGATE_LOOK, PROPAGATE_MORE, REBUILD_STROKE_SPAN,
			 REDO_FORWARD_SPAN, MOVE_TO_TOP, MOVE_TO_BOTTOM, MOVE_ONE_UP, MOVE_ONE_DOWN, MOVE_SPAN_TO_TOP, 
			 MOVE_SPAN_TO_BOTTOM, FIXD,
			 ACTIVATE_ALL, DELETE_ALL, DELETE, DELETE_SPAN,
			 SELECTD, SET_THICKNESS};


class Draw : public InterModule {

  Q_OBJECT

 public:
  
  Draw(ImgSequence *is, QGLWidget *parent);
  
  void paintGL();

  //void interpolatePatchForwards(DrawPatch* A); 

  void propagatePathEverywhere(DrawPath* path);
  void propagateLook(DrawPath* path);

  void frameChange(int i);
  
  //void visualize(const DrawPath* A, const RotoPath* A_prime,
  //	 const DrawPath* B, const RotoPath* B_prime);

  void toolChange(DrawTool id);
  DrawTool getTool() const { return _tool; }

  void setColor(const Vec3f& col) { _currColor = col; }
  
  CmdMap* getPopupFunctions(const Vec2i& loc, DynArray<int,5>* enab);
  void callPopup(const int which);

  void reportPressure(const int p) { _pressure = float(p) * .15f; }
  
  void runCan();

  void setdAlpha(const float v) { _dAlpha=v; }
  
  void reportZoom(const float v) { _zoom=v; }

  void sensChanged(const float v) {
    _sens = v;
  }

  void corrPropAll();

  void deleteSpan(DrawPath* p);
  void moveSpanToTop(DrawPath* p);
  void moveSpanToBottom(DrawPath* p);

  void setRenderMode(const int s) { _renderMode = s; }
  void setFreeDraw(const bool s) { _freeDraw = s; }

  void copy();
  void paste();

  void optimizeDrawShape(DrawPath* aPath, DrawPath* bPath, int aFrame, int bFrame, bool pinLast);

  void makeDrawPath(const RotoPath* rp);

  const DrawPath* selected() const { return _selected; }

  QRgb getCurrColor();

  bool preRegister;
 
  /* for textures */ 
  //  int getCurrTexture() { return _currTexture ; }
  void setCurrTexture( int i ) ; 

 signals:
  void currColorChanged();
  void restoreCursor();

 protected:
  
  void keyPressEvent ( QKeyEvent * e );
  void mousePressEvent( QMouseEvent *e );
  void mouseReleaseEvent( QMouseEvent *e );
  void mouseMoveEvent( QMouseEvent *e );

 private:

  float mapPressure() const;
  //void keyframeInterpolatePatch(DrawPatch *dch);
  //void longTrackSpan(DrawPatch *dch);
  void spliceDrawCurve();
  void finishSplicedCurve();
  void biGenerateSpan(DrawPath* aPath, const int aFrame, DrawPath* bPath, const int bFrame);
  void forwGenerateSpan(DrawPath* aPath, int aFrame);
  void regeneratePath(DrawPath* path);
  void finishJustPasted();
  void rebuildStrokeSpan(DrawPath* path);
  void redoForwardSpan(DrawPath* path);
  void captureColor();
  bool QtpickColor(Vec3f& col);
  void loadTexture( char* fName, int width, int height) ; 

  int _w, _h;
  ImgSequence *_is;
  int _frame;

  int _mousePressed;
  DrawCurves* _dc;
  DrawPath* _currDPath;
  bool _corrShow;
  QGLWidget *_parent;
  DrawPath* _selected;
  //DrawPatch* _selectedPatch;
  DrawPath* _popupPicked;
  //DrawPatch* _popupPatchPicked;
  Vec3f _currColor;
  int _canNum;
  float _dAlpha;
  int _renderMode;
  bool _freeDraw;
  DrawTool _tool;

  //DrawPatch* _propPatch;
  int _propFrame;

  float _pressure;
  float _zoom;
  float _sens;

  int _spliceIndex, _spliceFlag; // -1 not active, 1 on curve, 0 off,  -2 for needs correspondences

  DrawPath* _copyPtr;
  int _copyFrame;
  Vec2f _dragLoc;
  bool _dragDirty, _justPasted;
  bool _captureColor;
  Vec2i _captureColorLoc;

  // patch capturing
  //bool _capturePatch;
  //A
  //DrawPatch* _captureDPatch;

  // patch tracking
  //KLT_FullPyramid* _newPyramid;

};

#endif

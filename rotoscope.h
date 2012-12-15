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



#ifndef ROTOSCOPE_H
#define ROTOSCOPE_H

#include <qgl.h>
#include <qimage.h>
#include <vector>
#include <list>
#include <qpainter.h>
#include "interModule.h"
#include "imageGraph.h"
#include "liveWire.h"
#include "dynArray.h"
#include "jl_vectors.h"
#include "llist.h"
#include "imgSequence.h"
#include "rotoCurves.h"
#include "transformer.h"

enum TrackTool { T_DRAW, T_MANUAL, T_TRACK, T_SELECT, T_TRANSFORM, T_OVERDRAW, T_NUDGE, T_ZOOM };

enum RotoPopupFunction { CLEAR_ACTIVES, TRACK_RIGHT, TRACK_LEFT, TRACK_BOTH, TRACK_EDGES,
			 DELETE_PATH, LONG_TRACK_RSPAN, DELETE_RSPAN, SELECTR, CORRFORW, 
			 RPRINT, REJIG, R_FIX, MAKE_JOINTS, BREAK_JOINTS, SPLIT_SEGMENT, 
                         SPLIT_CURVE, FORM_REGION, PROPAGATE_REGION, UNFIX_POINT, TRI_WIDTH, MAKE_DRAW_PATH};


class Rotoscope : public InterModule {

Q_OBJECT

 public: 

  //Rotoscope(QImage im, QWidget* parent);
  Rotoscope(ImgSequence* is, QGLWidget* parent);
  ~Rotoscope();
  
  void setImageGraph(ImageGraph* ig);
  
  void paintGL();
  
  void mousePressEvent( QMouseEvent *e );
  void mouseReleaseEvent( QMouseEvent *e );
  void mouseMoveEvent( QMouseEvent *e );
  void keyPressEvent ( QKeyEvent *e );

  void frameChange(int i);
 
  void killAnyRoto();

  void toolChange(TrackTool id);
  TrackTool getTool() const { return _tool; }

  // enab lists item menus that should be DISABLED
  CmdMap* getPopupFunctions(const Vec2i& loc, DynArray<int,5>* enab);

  void callPopup(const int which);

  void makeVisible();
  void makeInvisible();

  void runCan();

  void calcLerp();

  void deleteSpan(RotoPath* path);

  void clearCorrForw() {
    _corrForw=NULL; _corrForwFrame = -1;
  }

  int selectedSize() const { return _selected.size(); }
  
  bool modTrackPossible() const;
  void rejig(PathV* paths,  int frame, bool startOver);
  void rejigWrapper(bool startOver);
  
  void copySplinesAcrossTime();

  void timerEvent (QTimerEvent *);

  void formRegion();
  void propagateRegion();

  void setViewRegion(const bool s);
  void setVisEffort(const bool s);

  void trackForward();

  bool _showTrackPoints;
  bool _nudgeFreezesCurve;

  const RotoPath* getPopupPicked() const { return _popupPicked; }
 signals:
  void rotoToDrawPath();

 private:

  // Procedures:
  void addPathSegment(); 
  void finishRotoCurve(); 
  void finishManualRotoCurve();
  //void spliceRotoCurve();
  //void finishSplicedCurve();
  void tryToMakeJoints();
  void tryToBreakJoints();
  void splitCurve();
  void splitSegment();
  void setTrimapWidth();
  void makeDrawPath();
  
  void trackForwards(RotoPath* A);
  //void keyframeInterpolatePath(RotoPath* aPath, int aFrame, RotoPath *bPath, int bFrame);
  //void longTrackSpan(RotoPath* rp);
  void intrinsicBlend(Vec2f* out, Vec2f* lengths, Vec2f* angles, int n, double t);
  void intrinsicBlend2(Vec2f* out, const Vec2f* earlier, const Vec2f& length, const Vec2f& angle, double t);
  void performTracks(const int aFrame, const int bFrame, bool doInterp=true, bool useExistingInbetweens=false);
  void keyframeSedInterp(RotoPath* aPath, int aFrame, RotoPath *bPath, int bFrame);
  void addMasksToMulti(MultiSplineData* mts, const PathV& key0, const int frame0);

  void restartPausedTracking();
  void signalTracking();

  // Data:
  QGLWidget* _parent;
  ImageGraph* _graph;
  LiveWire* _wire;
  ImgSequence* _is;

  RotoCurves* _rc;
  RotoPath* _currPath;
  Vec2i _seedPoint, _freePoint;
  int _frame;
  RotoPath *_popupPicked;
  RotoRegion *_regionPicked;
  Vec2f _popupLoc;
  int _popupCtrl;
  vector<RotoPath*> _selected; // we maintain uniqueness
  int _corrShow; // 0 for not, 1 for yes!
  int _w, _h;
  bool _iscissors;
  bool _visible;
  //int _propFrame;
  int _canNum;
  //int _spliceIndex, _spliceFlag; // -1 not active, 1 on curve, 0 off, -2 needs bezier
  int _lastFrameTouched, _changeLastFrameTouched;
  RotoPath* _shadowSelected;
  bool _drawing;
  RotoPath* _corrForw; int _corrForwFrame;
  TrackTool _tool;
  PathV _toTrack;
  RotoPath* _ctrlDrag; int _dragCtrlNum; // indice of control point being dragged
  list<TrackGraph*> _ccompV;
  //list<MultiTrackData*> _trackDataV;
  list<MultiSplineData*> _trackDataV;
  list<KLT_TrackingContext*> _TCV;
  bool _tracking;
  int _trackingTimer;
  Vec2f _dragLoc;
  bool _viewRegions;
  bool _visEffort;

  bool _manualStage;
  Transformer _transformer;
};


#endif

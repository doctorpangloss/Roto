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



#ifndef IMGSEQUENCE_H
#define IMGSEQUENCE_H

#include <fstream>
#include <string.h>
#include <assert.h>
#include <qimage.h>
#include <qpixmap.h>
#include <qfile.h>
#include <qtextstream.h>
#include <string>

#include "imageGraph.h"
#include "rotoCurves.h"
#include "drawCurves.h"
#include "KLT/klt.h"

class RotoCurves;
class DrawCurves;
class RotoPath;
class DrawPath;
//class DrawPatch;

// current file version
#define _VERSION_ 7        



class ImgSequence {

 public:

  ImgSequence(char* imageRoot,int startF, int endF);
  ImgSequence() {} // better be followed by a load

  int width() const { return _w; }
  int height() const { return _h; }


  QImage getImage(const int a) const;

  ImageGraph* getImageGraph(const int a) const;

  QImage getRelativeImage(const int a) const;

  ImageGraph* getRelativeImageGraph(const int a) const;

  RotoCurves* getRotoCurves(const int a) const;  
  
  RotoCurves* getRelativeRotoCurves(const int a) const;
  
  DrawCurves* getDrawCurves(const int a) const;
  
  DrawCurves* getRelativeDrawCurves(const int a) const;

  KLT_FullPyramid* getFullPyramid(const int a);
  KLT_FullCPyramid* getFullCPyramid(const int a);
  KLT_FullPyramid* getFullEdgePyramid(const int a);
  
  void writeFullCPyramid(const int a); 
  void writeFullEdgePyramid(const int a); 
  bool loadFullCPyramid(const int a);  
  bool loadFullEdgePyramid(const int a);  
  void loadCalcFullCPyramids(const int a, const int b); // inclusive range
  void loadCalcFullEdgePyramids(const int a, const int b); // inclusive range

  void writeFullCPyramidImage(const int a);
  void writeFullEdgePyramidImage(const int a);

  void loadInitialData(FILE* fp, const char* root);
  void loadRestData(FILE* fp);
  void loadInitialDataqt(QFile* fp, const char* root);
  void loadRestDataqt(QFile* fp);

  void saveCurvesAs(const char* _root);
  void saveCurvesQt();
  void saveCurves();

  void saveRotoXML(std::ofstream& fp);

  void saveMattes();

  int getStart() const { return _startF; }
  int getEnd() const { return _endF; }
  bool validFrameNum(const int i) const {
    return (i>=_startF && i<=_endF); }

  RotoPath* resolveRotoWritePtr(int* ptr);
  DrawPath* resolveDrawWritePtr(int* ptr);
  //DrawPatch* resolveDrawPatchWritePtr(int* ptr);

  const char* getFileRoot() { return _fileRoot; }

  bool fileRootInited();
  
  void calculateStatistics();
  
  KLT_TrackingContext globalTC;

  void nlevelsChanged(int res) { globalTC.nPyramidLevels = res; }
  void smooth0DerivChanged(float res) { globalTC.smooth0Deriv = res; }
  void smooth2DerivChanged(float res) { globalTC.smooth2Deriv = res; }
  void smooth1DerivChanged(float res) { globalTC.smooth1Deriv = res; }
  void shape2DerivChanged(float res) { globalTC.shape2Deriv = res; }
  void edgeWeightChanged(float res) { globalTC.edgeWeight = res; }
 private:

  void mutualInit();

  RotoCurves* _rotoc;
  DrawCurves* _drawc;
  KLT_FullPyramid* _pyrms;
  KLT_FullCPyramid* _Cpyrms;
  KLT_FullPyramid* _pyrmsE;

  char _imageRoot[200];
  char _fileRoot[200];
  mutable char _fileName[200];
  int _startF, _endF, _version; // currently 0,1,2 versions (highest latest)
  int _w, _h;

};




#endif


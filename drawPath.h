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



#ifndef DRAWPATH_H
#define DRAWPATH_H


#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#include <qdatastream.h>
#include "geigerCorresponder.h"
#include "rotoCurves.h"
#include "hertzStrokes/Stroke.h"
#include "Texture.h"
#include "imgSequence.h"
#include "bBoxf2D.h"
#include "drawContCorr.h"
#include <set>

typedef Vec2f DrawPoint; // can later make this more sophisticated
class RotoCurves;
class RotoPath;
class DrawContCorr;

#define HCOLOR(a) (int(a*255.))

class DrawPath : public AbstractPath {

 public:
  
  DrawPath();
  ~DrawPath();

  DrawPath(const DrawPath& other);

  void addVertex(const float x, const float y, const float p);

  void setColor(const Vec3f& c);
  void setFillColor(const Vec3f& c);
  bool filled() const { return _filled;}
  void setFilled(bool s) { _filled=s; }
  Vec3f getColor() const { return _strokeColor; }

  static void setCurrTexture( int i ) ; 
  /* for textures */ 
  //int getCurrTexture() { return _currTexture ; }

  void loadTexture( char* fName, int width, int height) ; 



  //  void setTextureName( GLuint texture_name ) ;  // also sets _textured to 'true' 
  //  void setTextured( bool textured ) ; 
  //bool textured() { return _textured ; } 
  //  GLuint getTextureName() { return _texture_name ; } 
  
  

  //void setStrokeRadius(const float a);



  void nudge(const int i, const Vec2f& loc);

  void redoStroke();
  void fair();

  void render(const int rmode, bool useAlpha=true) const;
  void renderCorr() const;
  void renderSelected() const;
  void renderRecent() const;

  void copyLook(const DrawPath* a);

  //bool correspondRoto(RotoCurves* rc);
  bool correspondRoto2(RotoCurves* rc, const int w, const int h);

  void vacateRotoCorr();
  bool corrToRoto() const { if (_corrRoto) return true; else return false;}

  const set<RotoPath*>& getRotoSet() const;

  //RotoPath* getCorrespondRoto() { return _corrRoto; }
  //void setCorrespondRoto(RotoPath* rt) { _corrRoto = rt; }

  //int getCorrNum() const { return _corrNum; }
  //void setCorrNum(const int i)  { _corrNum=i; }

  void getP(const int i, Vec3f& res) const; // location in coordinate frame 
  // of corrRoto at i (which is index of this)

  /*
  RotoPath* getRotoCorr(const int i, int& outi) const {
    assert(_corrs && i>=0 && i<getNumElements());
    outi = _corrs[i];
    return _corrsPtr[i];
    }*/

  const DrawContCorr* getCorrRoto() const { return _corrRoto; }
  
  double distToLast2(const Vec2f& loc) const;

  int justDrawn() const { return _justDrawn; }

  void setJustDrawn(int val) { _justDrawn = val; }

  void initHertzStroke();

  Stroke* getStroke() { return _stroke; }

  //void renderCorrPS(FILE* fp) const;

  //void rigFrom(RotoPath* rp, const int numPoints);
  void offsetRig(const float rT, RotoPath* corrRoto, const Vec3f P);
  void empty();

  DrawPath* nextC() { return _nextFrame; }
  void setNextC(DrawPath* next) {
    _nextFrame = next; }
  DrawPath* prevC() { return _prevFrame; }
  void setPrevC(DrawPath* prev) {
    _prevFrame = prev; }

  void save(FILE* fp, int version) const;
  void load(FILE* fp, int version);
  void saveqt(QDataStream* fp, int version) const;
  void loadqt(QDataStream* fp, int version);

  void incrementThickness(const float a);

  void decrementThickness(const float a);

  float  getThick(const int i) const {
    assert(i>=0 && i<getNumElements());
    return _thick.getElement(i);
  }

  void setUniformThickness(const float t);

  void calculateStrokeDisplayList();
  void calculateFillDisplayList();
  void freshenAppearance(); 

  Bboxf2D calcBbox() const;

  void interpolateForwards(int frame);
  static void fillForwardInterpolatedCurve(DrawPath* B, const DrawPath* A);
  static void fillBackwardInterpolatedCurve(DrawPath* B, const DrawPath* A);
  static void fillBiInterpolatedCurve(DrawPath* D0, DrawPath* D1, DrawPath *Dnew, double t, 
				      int *D0_T0_T1_i, RotoPath** D0_T0_T1_p,
				      int *D0_T0_T_new_i, RotoPath** D0_T0_T_new_p, int *T0_D1);


  //static void fillInterpolatedCurve(DrawPath* B, const DrawPath* A, RotoPath** BR, int *BRi);
  //static void fillForwardInterpolatedCurve(DrawPath* B, const DrawPath* A, RotoPath* B_prime,
  //				   const RotoPath* A_prime);
  //static void fillBackwardInterpolatedCurve(DrawPath* B, const DrawPath* A, RotoPath* B_prime,
  //				    const RotoPath* A_prime);
  Stroke* giveStroke() { return _stroke; }
  void forgetStroke() { _stroke = NULL; }

  void spliceIn(const DrawPath *other, const int i1, const int i2);

  void setFixed(const bool s);
  bool fixed() const { return _fixed; }

  void resample(const int* numSamples) { AbstractPath::resample(numSamples, &_thick); }

  void fixLoc();

  const AbstractPath* getShouldBe() { return _shouldBe; }
  void repToShouldBe();

  void takeT1D1(int* a, int n) { _T1_D1 = a; _T1_D1_n = n; }
  int* getT1D1(int* n) { *n=_T1_D1_n; return _T1_D1; }

  DrawPath *prev,*next;

  static float _dAlpha;

  // global info
  // the info to use when the next stroke is created
  static bool _nxt_useTexture; 
  static Texture_Type _nxt_textureType; 
  static int _nxt_textureWidth; 
  static int _nxt_textureHeight; 
  static GLuint _nxt_texture_name;
  static int _nxt_texture_index; 
  static float _nxt_textureRadius; 
  static GLuint _texture_names [NUM_TEXTURES]; 


 private:


  //void sample5(Vec2f* here) const;

  //RotoPath* _corrRoto;
  //int _corrNum; // order in _corrRoto's list of attached drawPaths
  DrawPath *_prevFrame, *_nextFrame; 

  // correspondence data
  //double* _corrs; // maps this to corrRoto
  //RotoPath** _corrsPtr;
  DrawContCorr* _corrRoto;   // NEED TO SAVE/WRITE

  Vec3f _strokeColor; 
  bool _filled; 
  Vec3f _fillColor; 
  float _thicknessOffset; 
  GLuint _dFillList, _dStrokeList;
  DynArray<float,100> _thick; 


  // skq 
  bool _textured ; 
  int _texture ;  // index into draw.C texture array 
  int _texture_width , _texture_height ; 
  Texture_Type _type ;


  


  Stroke *_stroke;
  int _justDrawn;
  bool _fixed;

  AbstractPath* _shouldBe; // NEWDATA; must save, load
  int* _T1_D1; // NEWDATA
  int _T1_D1_n;
};

#endif

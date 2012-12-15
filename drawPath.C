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
#include "drawPath.h"
#include "drawCorresponder.h"


#define defaultThickness 2

bool DrawPath::_nxt_useTexture = false ;  
int  DrawPath::_nxt_textureWidth = 0 ;  
int  DrawPath::_nxt_textureHeight = 0 ;
float DrawPath::_nxt_textureRadius = 1.0 ; 
Texture_Type DrawPath::_nxt_textureType = none; 
GLuint DrawPath::_nxt_texture_name = 999; 
int DrawPath::_nxt_texture_index = -1; 
GLuint DrawPath::_texture_names [NUM_TEXTURES] = { 0, 0, 0, 0, 0, 0 } ; 


int modulo(const int a, const int b) {
  int tmp = a%b;
  return (tmp>=0 ? tmp : (b)+tmp);
}

DrawPath::DrawPath() : AbstractPath() {
  //_corrRoto = NULL; 
  //_corrs = NULL; 
  //_corrsPtr = NULL;
  _corrRoto = NULL;
  _justDrawn = 1; 
  _stroke = NULL;
  _nextFrame = _prevFrame = NULL;
  //_corrNum = -1;
  _filled = false;
  _dFillList = 0; _dStrokeList=0;
  _thicknessOffset = 0;
  _fixed = true;
  _shouldBe = NULL;
  _T1_D1 = NULL;
  _T1_D1_n = -1;

  _textured = false; 
  _texture  = -1 ; 
 
}

DrawPath::DrawPath(const DrawPath& other) : AbstractPath(other) {
  //_corrs = NULL; 
  //_corrsPtr = NULL;
  _corrRoto = NULL;
  _justDrawn = 1; 
  _stroke = other._stroke;
  _filled = other._filled;
  _dFillList = 0; _dStrokeList=0;
  _thicknessOffset = other._thicknessOffset;
  _nextFrame = _prevFrame = NULL;
  _fixed = other._fixed;
  _thick = other._thick;
  _strokeColor = other._strokeColor;
  _fillColor = other._fillColor;
  _shouldBe = NULL;
  _T1_D1 = NULL;
  _T1_D1_n = -1;

  // skq
  _textured =  other._textured; 
  _texture = other._texture ; 
  
  calculateStrokeDisplayList();
  calculateFillDisplayList();
}


void DrawPath::addVertex(const float x, const float y, const float p) {
  /*
    float myP;
  if (p==0)    // mouse stroke
    myP = defaultThickness;
  else
    myP = p;
  */

  DrawPoint dp(x,y);
  add(dp);
  _thick.add(p);
  if (_stroke)
    _stroke->addControlPoint(x,y,p); 
}

#define _alpha_ (useAlpha ? _dAlpha : 1.f)

void DrawPath::render(const int rmode, bool useAlpha) const {

  bool normal = (rmode==0 || rmode==1);

  if (_inMotion) {
    glPushMatrix();
    glTranslatef(_motion.x(), _motion.y(), 0);
  }
  
  if (_dFillList != 0 && _filled && normal) {  // fill display list 
    glColor4f(_fillColor.r(), _fillColor.g(), _fillColor.b(), _alpha_);
    glCallList(_dFillList);
  }

  // stroke color !
  if (_dStrokeList!=0 && normal) {  // stroke display list
    glColor4f(_strokeColor.r(), _strokeColor.g(), _strokeColor.b(), _alpha_);

    glCallList(_dStrokeList);
  }
  else  if (_stroke && normal) {      // hertzmann stroke
    glColor4f(_strokeColor.r(), _strokeColor.g(), _strokeColor.b(), _alpha_); 
     _stroke->drawWideLineCurve(); 
  }
  
  else {                  // broke-ass GL rendering
    glColor4f(_strokeColor.r(), _strokeColor.g(), _strokeColor.b(), _alpha_);
    glBegin(GL_POINTS);
    for (int i=0; i<getNumElements(); i++) {
      const DrawPoint *dp = getPointerToElement(i);
      glVertex2f(dp->x(),dp->y());
    }
    glEnd();

    if (_shouldBe) {
      glColor4f(1,0,0,_alpha_);
      glBegin(GL_POINTS);
      for (int i=0; i<getNumElements(); i++) {
	const Vec2f* dp = _shouldBe->getPointerToElement(i);
	glVertex2f(dp->x(),dp->y());
      }
      glEnd();
    }
  }

  if (_inMotion)
    glPopMatrix();

}
  
void DrawPath::renderRecent() const {
  int i = getNumElements();
  if ( i<2) return;
  glGetError();
  glBegin(GL_LINES);
  const DrawPoint *dp = getPointerToElement(i-1), *dp2 = getPointerToElement(i-2);
  glVertex2f(dp->x(),dp->y());
  printf("%f %f\n",dp->x(),dp->y());
  glVertex2f(dp2->x(),dp2->y());
  printf("%f %f\n",dp2->x(),dp2->y());
  glEnd();
  assert(!glGetError());
}



void DrawPath::renderCorr() const {
  
  if (!_corrRoto) return;
  Vec2f a,b;
  glLineWidth(1);
  glBegin(GL_LINES);
  for (int i=0; i<getNumElements(); i+=3) {
    a = getElement(i);
    glVertex2f(a.x(),a.y());
    b = _corrRoto->getSampleRoto(i)->getLoc(_corrRoto->getSampleRotoT(i));
    glVertex2f(b.x(),b.y());
  }
  glEnd();

  /*  if (!_corrRoto) return;
  glLineWidth(1);
  glBegin(GL_LINES);
  Vec2f a,b;
  for (int i=0; i<getNumElements(); i+=3) {
    a = _corrRoto->getElement(_corrs[i]);
    glVertex2f(a.x(),a.y());
    b= getElement(i);
    glVertex2f(b.x(),b.y());
  }
  glEnd();*/

  /*
  if (!(_corrs && _corrsPtr)) return;
  glLineWidth(1);
  glBegin(GL_LINES);
  Vec2f a,b;
  for (int i=0; i<getNumElements(); i+=1) {
    a = _corrsPtr[i]->getElement(_corrs[i]);
    glVertex2f(a.x(),a.y());
    b= getElement(i);
    glVertex2f(b.x(),b.y());

  }

  glEnd();
  */
}

void DrawPath::renderSelected() const {
  /*  Vec2f loc = getElement(0);
  glRectf(loc.x()-2, loc.y()-2, 
	  loc.x()+2, loc.y()+2);
  loc = getElement(getNumElements()-1);
  glRectf(loc.x()-2, loc.y()-2, 
  loc.x()+2, loc.y()+2);*/
  glColor3f(153./256., 51./256., 153./256.);

  //printf("render selected called\n") ;
  if (_textured ) { 
    setCurrTexture(_texture); 
  } 

  //printf("starting stroke render\n") ; 
  _stroke->render();
  //printf("make it out of stroke render\n") ; 
}

/*
void DrawPath::sample5(Vec2f* here) const {
 
  here[0] = getElement(0); 
  here[4] = getElement(getNumElements()-1);;
  here[1] = getElement(int(.25*float(getNumElements())));
  here[2] = getElement(int(.50*float(getNumElements())));
  here[3] = getElement(int(.75*float(getNumElements())));
}
*/

const set<RotoPath*>& DrawPath::getRotoSet() const {
  return _corrRoto->getRotoSet();
}

void DrawPath::vacateRotoCorr() {
  if (!_corrRoto) return;
  const set<RotoPath*>& rotoset = _corrRoto->getRotoSet();
  for (set<RotoPath*>::const_iterator c = rotoset.begin(); c != rotoset.end(); ++c)
    (*c)->notifyDrawDelete(this);  
  _corrRoto->vacateRotoCorr();

  /*
  set<RotoPath*>* rotos = getRotoSet();
  if (rotos) {
    for (set<RotoPath*>::iterator c = rotos->begin(); c != rotos->end(); ++c)
      (*c)->notifyDrawDelete(this);  
    delete rotos;
  }
  if (_corrsPtr) {
    delete[] _corrsPtr;  _corrsPtr = NULL;
  }
  if (_corrs) {
    delete[] _corrs; _corrs = NULL;
  }
  */
}


/*
TRULY DEAD
bool DrawPath::correspondRoto(RotoCurves* rc) {
  printf("correspond roto\n");
  // find RotoPath to rig self to
  assert(!_corrRoto); 
  rc->startCurvesIterator();
  RotoPath* currRoto;
  double dist=DBL_MAX,tmp;
  int start_i, end_i,ai,bi,i, shit;
  Vec2f drawLoc, samples[5];
  sample5(samples);
  while ((currRoto = rc->getNextCurve()) != NULL) {
    tmp = 0;
    tmp = currRoto->distanceTo2(samples[0],ai);
    tmp += currRoto->distanceTo2(samples[4],bi);
    for (i=1; i<4; i++)
      tmp += currRoto->distanceTo2(samples[i],shit);
    if (tmp < dist) {
      _corrRoto = currRoto;
      dist = tmp;
      start_i = ai;
      end_i = bi;
    }
  }
  if (!_corrRoto) {
    printf("distance too great 0\n");
    return 0;
  }

  int dCorrs, startCorrs;
  int corrCount = _corrRoto->getNumElements();
  if (dist==DBL_MAX || start_i==end_i) {
    _corrRoto = NULL;
    printf("distance too great 1\n"); assert(0); 
    return 0;
  }
  _corrNum = _corrRoto->addCorrDraw(this);
  
  bool rotoClosed = false;
  if (_corrRoto->getElement(corrCount-1) == _corrRoto->getElement(0))
    rotoClosed = true;

  // get subpath from _corrRoto
  AbstractPath subPath, subPathR;
  assert(end_i != start_i);
  if (end_i>start_i) dCorrs = 1;  //direction along roto curve...
  else dCorrs = -1;
  i=start_i;
  do {
    subPath.add(_corrRoto->getElement(i));
    i+=dCorrs;
  } while (i!=end_i);
  startCorrs = start_i;
  
  // other direction on closed rotoCurves
  if (rotoClosed) {
    i=start_i;
    do {
      printf("%d\n",i);
      subPathR.add(_corrRoto->getElement(i));
      i = modulo((i-dCorrs),corrCount);
    } while (i!=end_i);
  }
  
  double aCost, bCost;
  // corespond subpath
  assert(!_corrs); 
  _corrs = new int[getNumElements()];
  GeigerCorresponder *myCoor, *myCoorR=NULL;
  myCoor = new GeigerCorresponder(this,&subPath);
  aCost = myCoor->calculate();
  if (rotoClosed) {
    myCoorR = new GeigerCorresponder(this,&subPathR);
    bCost = myCoorR->calculate();
  }
  if (!rotoClosed || aCost < bCost) {
    myCoor->getSolution(_corrs);
    for (i=0; i<getNumElements(); i++) // convert _corr indices to actual
      _corrs[i] = startCorrs + dCorrs*_corrs[i];
  }
  else {
    printf("backwards\n");
    myCoorR->getSolution(_corrs);
    for (i=0; i<getNumElements(); i++) // convert _corr indices to actual
      _corrs[i] = modulo((startCorrs - dCorrs*_corrs[i]),corrCount);
  }
  
  //printf("roto: %d, drawn: %d\n", subPath.getNumElements(), getNumElements());
  
  for (int j=0; j<getNumElements(); j++)
    printf("%d ",_corrs[j]);
  printf("\n");
  

  delete myCoor;
  if (myCoorR) delete myCoorR;

  return 1;

}
*/


DrawPath::~DrawPath() {
  if (_corrRoto) delete _corrRoto;
  if (_stroke) delete _stroke;
  if (_shouldBe) delete _shouldBe;
  if (_T1_D1) delete[] _T1_D1;
}

void DrawPath::setColor(const Vec3f& c) { 
  _strokeColor = c; 
  //calculateDisplayList();
}


/*void DrawPath::rigFrom(RotoPath* rp, const int numPoints) {
  //_corrRoto = rp;
  //_corrNum = _corrRoto->addCorrDraw(this);
  _corrs = new int[numPoints];
  _corrsPtr = new RotoPath*[numPoints];
}
*/

double DrawPath::distToLast2(const Vec2f& loc) const {
  return getPointerToElement(getNumElements()-1)->distanceTo2(loc);
}


void DrawPath::empty() {
  resetElements();
  _thick.resetElements();
  if (_corrRoto) {
    delete _corrRoto; _corrRoto = NULL;
  }
  if (_shouldBe) {
    delete _shouldBe; _shouldBe = NULL;
  }
  if (_stroke) {
    delete _stroke;
    _stroke = NULL;
  }
}

void DrawPath::redoStroke() {
  
  if (!_stroke) {
    printf("No stroke to redo\n");
    return;
  }
  _stroke->clear();
  Vec2f loc;
  for (int i=0; i<getNumElements(); i++) {
    loc = getElement(i);
    _stroke->addControlPoint(loc.x(),loc.y(),_thick.getElement(i)); 
  }
}

void DrawPath::fair() {
  AbstractPath::fair();
  redoStroke();
}

void DrawPath::save(FILE* fp, int version) const {
  int dummy;
  
  dummy = getNumElements();
  fwrite(&dummy,sizeof(int),1,fp); // write # of points
  fwrite(getData(),sizeof(Vec2f),dummy,fp); // write points
  assert(_thick.getNumElements() == dummy);
  fwrite(_thick.getData(), sizeof(float), dummy, fp); // write thicknesses
  fwrite(&_fixed, sizeof(bool), 1, fp); // _fixed

  if (_prevFrame)    // _prevFrame
    _prevFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  if (_nextFrame)    // _nextFrame
    _nextFrame->writePtr(fp);
  else
    writeNullPtr(fp);

  /*
  dummy=0;   // _corrs intArray
  if (_corrs)
    writeIntArray(fp, _corrs, getNumElements());
  else
    fwrite(&dummy,sizeof(int),1,fp);

  if (_corrsPtr) {   // _corrsPtr array
    dummy = getNumElements();
    fwrite(&dummy,sizeof(int),1,fp);
    for (int i=0; i<dummy; i++)
      _corrsPtr[i]->writePtr(fp);
  }
  else {
    dummy=0;
    fwrite(&dummy,sizeof(int),1,fp);
  }
  */

  if (_shouldBe) { // _shouldBe
    dummy = _shouldBe->getNumElements();
    fwrite(&dummy,sizeof(int),1,fp);
    fwrite(_shouldBe->getData(),sizeof(Vec2f),dummy,fp);
  }
  else {
    dummy=0;
    fwrite(&dummy,sizeof(int),1,fp);
  }

  dummy=0;   // T1_D1
  if (_T1_D1)
    writeIntArray(fp, _T1_D1, _T1_D1_n);
  else
    fwrite(&dummy,sizeof(int),1,fp);

  fwrite(_strokeColor.data(), sizeof(Vec3f), 1, fp); // _strokeColor
  fwrite(_fillColor.data(), sizeof(Vec3f), 1, fp); //_fillColor
  fwrite(&_filled, sizeof(bool), 1, fp); // _filled

  printf("saving: textured %d, texture %d\n", _textured, _texture); 

  fwrite(&_thicknessOffset, sizeof(float),1,fp); // _thicknessOffset

  fwrite(&_textured, sizeof(bool),1,fp); // textured
  fwrite(&_texture,  sizeof(int),1,fp);  // texture index


}

void DrawPath::saveqt(QDataStream* fp, int version) const {
  int dummy,i;
  dummy = getNumElements();
  *fp << dummy;               // write # of points
  for (i=0; i<dummy; ++i)
    *fp << getData()[i].x() << getData()[i].y();     // write points
  writeBool(fp, _fixed);
  //*fp << (Q_INT8) _fixed;  // _fixed;

  assert(_thick.getNumElements() == dummy);
  for (i=0; i<dummy; ++i) {
    *fp << _thick.getData()[i];
  }

  if (_prevFrame)    // rotopath _prevFrame pointer
    _prevFrame->writePtr(fp);
  else
    writeNullPtr(fp);
  
  if (_nextFrame)    // _nextFrame pointer 
    _nextFrame->writePtr(fp);
  else
    writeNullPtr(fp);

  dummy=0;
  /*
  if (_corrs)     // _corrs array
    writeIntArray(fp,_corrs, getNumElements());
  else
    *fp << dummy;

  if (_corrsPtr) {   // _corrsPtr array
    dummy = getNumElements();
    *fp << dummy;
    for (int i=0; i<dummy; i++)
      _corrsPtr[i]->writePtr(fp);
  }
  else {
    dummy = 0; *fp << dummy;
  }
  */

  if (_shouldBe) { // _shouldBe
    dummy = _shouldBe->getNumElements();
    *fp << dummy;               // write # of points
    for (i=0; i<dummy; ++i)
      *fp << _shouldBe->getData()[i].x() << _shouldBe->getData()[i].y();     // write points  
  }
  else {
    dummy = 0; *fp << dummy;
  }

  dummy=0;
  if (_T1_D1)     // T1_D1
    writeIntArray(fp, _T1_D1, _T1_D1_n);
  else
    *fp << dummy;

  *fp << _strokeColor.r() << _strokeColor.g() << _strokeColor.b();
  *fp << _fillColor.r() << _fillColor.g() << _fillColor.b();
  //*fp << (Q_INT8)_filled;
  writeBool(fp, _filled);
  *fp << _thicknessOffset;
  //*fp << (Q_INT8)_textured;
  writeBool(fp,_textured);
  *fp << _texture;
}



void DrawPath::load(FILE* fp, int version) {
  int dummy, dptr[2], res;
  
  res = fread(&dummy,sizeof(int),1,fp); // read # of points
  

  assert(res == 1);
  addNElements(dummy);
  res = fread(getData(),sizeof(Vec2f),dummy,fp); // read points
  assert(res==dummy);
  _thick.addNElements(dummy);
  res = fread(_thick.getData(), sizeof(float), dummy, fp); // read thickness
  assert(res==dummy);
  res = fread(&_fixed, sizeof(bool), 1, fp); // _fixed
  assert(res==1);

  res = fread(dptr,sizeof(int),2,fp); // _prevFrame pointer
  assert(res == 2);
  if (dptr[0] == -1)
    _prevFrame = NULL;
  else
    _prevFrame = _is->resolveDrawWritePtr(dptr);

  res = fread(dptr,sizeof(int),2,fp); // _nextFrame pointer
  assert(res == 2);
  if (dptr[0] == -1)
    _nextFrame = NULL;
  else
  _nextFrame = _is->resolveDrawWritePtr(dptr);

  if (version < 4) {
    res = fread(&dummy,sizeof(int),1,fp);  // _corrs intArray
    assert(res == 1);
    fseek(fp,dummy*sizeof(int),1);
    /*if (dummy==0)
      _corrs = NULL;
    else {
      _corrs = new int[dummy];
      fread(_corrs,sizeof(int),dummy,fp);
      }*/
    
    
    res = fread(&dummy, sizeof(int),1,fp); // _corrsPtr array
    assert(res==1);
    fseek(fp,dummy*sizeof(RotoPath*), 1);
    /*if (dummy==0)
      _corrsPtr = NULL;
    else {
      _corrsPtr = new RotoPath*[dummy];
      for (int i=0; i<dummy; i++) {
	res = fread(dptr,sizeof(int),2,fp);
	_corrsPtr[i] = _is->resolveRotoWritePtr(dptr);
      }
      }*/
  }
  
  
  res = fread(&dummy, sizeof(int),1,fp); // _shouldBe
  assert(res == 1);
  if (dummy>0) {
    _shouldBe = new AbstractPath();
    _shouldBe->addNElements(dummy);
    res = fread(_shouldBe->getData(),sizeof(Vec2f),dummy,fp); // read points
    assert(res==dummy);
  }
  else
    _shouldBe = NULL;

  res = fread(&dummy,sizeof(int),1,fp);  // T1_D1
  assert(res == 1);
  if (dummy==0)
    _T1_D1 = NULL;
  else {
    _T1_D1 = new int[dummy];
    _T1_D1_n = dummy;
    fread(_T1_D1,sizeof(int),dummy,fp);
  }
  
  fread(_strokeColor.data(), sizeof(Vec3f), 1, fp); // _strokeColor
  fread(_fillColor.data(), sizeof(Vec3f), 1, fp); //_fillColor
  fread(&_filled, sizeof(bool), 1, fp); // _filled
  fread(&_thicknessOffset, sizeof(float),1,fp); // _thicknessOffset

  printf("version : %d\n", version); 

  // textures were added in version 2 by SKQ
  if ( version > 1 ) { 
    //skq

    // this is not reading right? 
    fread(&_textured, sizeof(bool),1,fp); // 
    fread(&_texture,  sizeof(int),1,fp);  // texture index

    printf("LOAD:_textured = %d\n", _textured); 
    printf("LOAD:_texture  = %d\n", _texture ) ; 
    
    // debug 
    
    if ( _texture < 0 || _texture >= NUM_TEXTURES ) { 
      _textured = false; 
      _texture = -1; 
      printf("LOAD:_texture  = %d\n", _texture ) ; 
    }

    if ( _textured ) { 
      
      setCurrTexture( _texture ) ; 
    } 
    
  } else {
    _textured = false; 
    _texture  = -1; 
  }

  
  initHertzStroke();


  Vec2f loc;
  for (int i=0; i<getNumElements()-1; i++) {
    loc = getElement(i);
    _stroke->addControlPoint(loc.x(),loc.y(),_thick.getElement(i)); 
  }

  calculateStrokeDisplayList();
  calculateFillDisplayList();
}



void DrawPath::loadqt(QDataStream* fp, int version) {
  int dummy,i, dptr[2];
  
  *fp >> dummy;
  addNElements(dummy);
  for (i=0; i<dummy; ++i)
    *fp >> getData()[i].data()[0] >> getData()[i].data()[1];  

  //*fp >> (Q_INT8&) _fixed;  // _fixed;
  readBool(fp, _fixed);

  _thick.addNElements(dummy); // thick
  //float* test = new float[dummy];
  for (i=0; i<dummy; ++i) {
    *fp >> _thick.getData()[i];
  }



  *fp >> dptr[0] >> dptr[1]; // _prevFrame
  if (dptr[0] == -1)
    _prevFrame = NULL;
  else
    _prevFrame = _is->resolveDrawWritePtr(dptr);
  
  *fp >> dptr[0] >> dptr[1]; // _nextFrame
  if (dptr[0] == -1)
    _nextFrame = NULL;
  else
    _nextFrame = _is->resolveDrawWritePtr(dptr);

  if (version < 4) {
    *fp >> dummy; // _corrs
    /*
      if (dummy==0)
      _corrs = NULL;
      else {
      _corrs = new int[dummy];
      readIntArray(fp,_corrs,dummy);
      }
    */
    readIntArray(fp,NULL,dummy);

    *fp >> dummy;  // _corrsPtr
    for (int i=0; i<dummy; i++)
      *fp >> dptr[0] >> dptr[1];
    /*if (dummy==0)
      _corrsPtr = NULL;
    else {
      _corrsPtr = new RotoPath*[dummy];
      for (int i=0; i<dummy; i++) {
	*fp >> dptr[0] >> dptr[1];
	_corrsPtr[i] = _is->resolveRotoWritePtr(dptr);
      }
      }*/
  }

   *fp >> dummy; // _shouldBe
   if (dummy>0) {
     _shouldBe = new AbstractPath();
     _shouldBe->addNElements(dummy);
     for (i=0; i<dummy; ++i)
       *fp >> _shouldBe->getData()[i].data()[0] >> _shouldBe->getData()[i].data()[1];  
   }
   else
    _shouldBe = NULL;

   *fp >> dummy; // T1_D1
   if (dummy==0)
     _T1_D1 = NULL;
   else {
     _T1_D1 = new int[dummy];
     _T1_D1_n = dummy;
     readIntArray(fp,_T1_D1,dummy);
   }

   *fp >> _strokeColor.data()[0] >> _strokeColor.data()[1] >> _strokeColor.data()[2];
  *fp >> _fillColor.data()[0] >> _fillColor.data()[1] >> _fillColor.data()[2];
  //*fp >> (Q_INT8&)_filled;
  readBool(fp, _filled);
  *fp >> _thicknessOffset;

  if (version > 1) {
    readBool(fp, _textured);
    //*fp >> (Q_INT8&)_textured;
    *fp >> _texture;
  }

  if ( _texture < 0 || _texture >= NUM_TEXTURES ) { 
    _textured = false; 
    _texture = -1; 
  }

  if ( _textured )
      setCurrTexture( _texture ) ; 

  else {
    _textured = false; 
    _texture  = -1; 
  }

  
  initHertzStroke();


  Vec2f loc;
  for (int i=0; i<getNumElements()-1; i++) {
    loc = getElement(i);
    _stroke->addControlPoint(loc.x(),loc.y(),_thick.getElement(i)); 
  }

  calculateStrokeDisplayList();
  calculateFillDisplayList();
}
  


void DrawPath::initHertzStroke() {
  assert(!_stroke);

  _stroke = new Stroke(defaultThickness + _thicknessOffset, 1,
  //	       makeColor(HCOLOR(_strokeColor.r()),
  //			 HCOLOR(_strokeColor.g()),
  //			 HCOLOR(_strokeColor.b())), //(GLubyte)0x80), 
		       0);
  _stroke->numLevels() = 0;
  

  /* figure out a better way to set these variables... */
  if ( _nxt_useTexture ) { 
    _stroke->useTexture() = _nxt_useTexture ; 
    _stroke->textureWidth() = _nxt_textureWidth; 
    _stroke->textureHeight() = _nxt_textureHeight ;
    _stroke->textureName() = _nxt_texture_name; 
    _stroke->textureType() = _nxt_textureType; 
    _stroke->textureRadius() = _nxt_textureRadius; 

    _textured = true ;
    _texture = _nxt_texture_index; 

  }else{ 

    // unset use texture 
    _stroke->useTexture() = false ; 
    _stroke->textureWidth() = 1 ; 
    _stroke->textureHeight() = 1 ; 
    _stroke->textureName() = 999 ; 
    _stroke->textureType() = none ; 
    _stroke->textureRadius() = 1.0; 
    
    _textured = false; 
    _texture = -1; 

  } 

}



void DrawPath::setFillColor(const Vec3f& c) {
  _fillColor = c;  
  _filled = true; 
  calculateFillDisplayList();
}

/*******
void DrawPath::setTextureName ( GLuint texture_name ) { 
  _texture_name = texture_name ; 
  _textured = true; 
}

void DrawPath::setTextured( bool textured ) { 
  _textured = textured ; 
  }********/

void fuck(GLenum errno) {
  const GLubyte *errString = gluErrorString(errno);
  printf("OpenGL error: %s\n",errString);
  exit(0);
}
void fuckCombine(GLdouble coords[3], GLdouble* vd[4], GLfloat weight[4], void **outData) {
  GLdouble *vertex = (GLdouble *) malloc(3*sizeof(GLdouble));
  memcpy(vertex, coords, 3*sizeof(GLdouble));
  *outData = vertex;

  //printf("vertex data\n");
  //for (int i=0; i<4; i++)
  //printf("%f %f %f\n",vd[i][0], vd[i][1], vd[i][2]);
  //assert(0);
}

void DrawPath::calculateStrokeDisplayList() {
  if (!_stroke) return;

  glColor3f(_strokeColor.r(), _strokeColor.g(), _strokeColor.b() ); 

  if (_dStrokeList == 0)
     _dStrokeList = glGenLists(1);
  assert(glGetError() == GL_NO_ERROR);
  assert(_dStrokeList != 0);
  
  glNewList(_dStrokeList, GL_COMPILE);
  _stroke->toffset() = _thicknessOffset;
  

  _stroke->render();
  glEndList();
  
  assert(glGetError() == GL_NO_ERROR);  
}

void DrawPath::calculateFillDisplayList() {
  if (!_filled) return;
  if (getNumElements() < 3) {
    _filled = false;
    return;
  }

  if (_dFillList == 0)
    _dFillList = glGenLists(1);
  assert(glGetError() == GL_NO_ERROR);
  assert(_dFillList != 0);

  glNewList(_dFillList, GL_COMPILE);
  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // copy our data into 3 doubles for GLU
  GLdouble *v = new GLdouble[getNumElements()*3];
  Vec2f ourLoc;
  for (int i=0; i<getNumElements(); i++) {
    ourLoc = getElement(i);
    v[i*3] = ourLoc.x(); v[i*3+1] = ourLoc.y(); v[i*3+2] = 0;
    //printf("%f %f\n",ourLoc.x(), ourLoc.y());
  }
  
  GLUtesselator* tobj = gluNewTess();
  gluTessCallback(tobj, GLU_TESS_VERTEX, 
		  (GLvoid (*) ()) &glVertex3dv);
  gluTessCallback(tobj, GLU_TESS_BEGIN, 
		  (GLvoid (*) ()) &glBegin);
  gluTessCallback(tobj, GLU_TESS_END, 
		  (GLvoid (*) ()) &glEnd);
  gluTessCallback(tobj, GLU_TESS_ERROR, 
		  (GLvoid (*) ()) &fuck);
  gluTessCallback(tobj, GLU_TESS_COMBINE, 
		  (GLvoid (*) ()) &fuckCombine);
  
  gluTessBeginPolygon(tobj, NULL);
  gluTessBeginContour(tobj);
  
  for (int i=0; i<getNumElements(); i++)
    gluTessVertex(tobj, v+i*3, v+i*3);
  
  gluTessEndContour(tobj);
  gluTessEndPolygon(tobj);
  gluDeleteTess(tobj);

  
  assert(glGetError() == GL_NO_ERROR);
    
  delete[] v;
  glEndList();
}

void DrawPath::freshenAppearance() {
  if (_stroke) {
    //_stroke->radius() = defaultThickness + _thicknessOffset;
    _stroke->toffset() = _thicknessOffset;
    /*
    _stroke->color() = makeColor(HCOLOR(_strokeColor.r()),
				 HCOLOR(_strokeColor.g()),
				 HCOLOR(_strokeColor.b())); 
    */
  }
  calculateStrokeDisplayList();
  calculateFillDisplayList();
}

/*
void DrawPath::setStrokeRadius(const float a) {
  if (!_stroke) return;
  _stroke->radius() = MAX(0, a);
  calculateStrokeDisplayList();
}
*/

void DrawPath::setUniformThickness(const float t) {
  //_stroke->setUniformRadius(t);
  //calculateStrokeDisplayList();
  for (int i=0; i<getNumElements(); i++)
    _thick.setElement(i,t);
  redoStroke();
  calculateStrokeDisplayList();
}

void DrawPath::incrementThickness(const float a) {
  _thicknessOffset += a;
  if (_stroke)
    _stroke->toffset() = _thicknessOffset;
  //_stroke->radius() = defaultThickness + _thicknessOffset;
  calculateStrokeDisplayList();
}

void DrawPath::decrementThickness(const float a) {
  _thicknessOffset -= a;
  if (_stroke)
    _stroke->toffset() = _thicknessOffset;
  //_stroke->radius() = MAX(0, defaultThickness + _thicknessOffset);
  calculateStrokeDisplayList();
}


void DrawPath::copyLook(const DrawPath* a) {
  _strokeColor = a->_strokeColor;
  _filled = a->_filled;
  _fillColor = a->_fillColor;
  _thicknessOffset = a->_thicknessOffset;
}


Bboxf2D DrawPath::calcBbox() const {
  Bboxf2D bbox;
  for (int i=0; i<getNumElements(); i++) {
    Bboxf2D currBox;
    Vec2f currLoc = getElement(i), extLoc, 
      thick(_thick.getElement(i), _thick.getElement(i));
    
    Vec2f_Sub(extLoc,currLoc,thick);
    currBox.includePoint(extLoc);
    
    Vec2f_Add(extLoc,currLoc,thick);
    currBox.includePoint(extLoc);
    
    bbox.setToUnion(currBox);
  }

  return bbox;
}


void DrawPath::interpolateForwards(int frame) {
  /*
    NOT SO IMPORTANT

  DrawPath *A=this;

  RotoPath** nextRotos = new RotoPath*[getNumElements()];
  int i;
  for (i=0; i<getNumElements(); i++) {
    nextRotos[i] = _corrsPtr[i]->nextC();
    if (!nextRotos[i]) {
      delete[] nextRotos;
      return;
    }
  }

  // if later drawPath exists, clear it, otherwise create new one
  DrawPath* B;
  if ((B=nextC())) {
    B->empty();
    assert(B->prevC() == this);
  }
  else {
    B = new DrawPath();
    _is->getDrawCurves(frame)->addPath(B);
    setNextC(B);
    B->setPrevC(this);
  }

  // create interpolated drawPath
  B->copyLook(this);
  fillForwardInterpolatedCurve(B,A, nextRotos);
  
  delete[] nextRotos;
  */
}

/*
  PROBABLY TRULY DEAD
void DrawPath::fillInterpolatedCurve(DrawPath* B, const DrawPath* A, RotoPath** BR, int *BRi) {
  B->initHertzStroke();
  int A_prime_index, B_prime_index;
  RotoPath *A_prime, *B_prime;
  Vec3f P;
  int i=0;
  int shit = A->getNumElements();
  
  //B->rigFrom(B_prime, A->getNumElements());
  B->_corrs = new int[A->getNumElements()];
  B->_corrsPtr = new RotoPath*[A->getNumElements()];

  //const RotoPath* curr;
  while (i < A->getNumElements()) {
    A_prime_index = A->_corrs[i];
    A_prime = A->_corrsPtr[i];  assert(A_prime);
    A->getP(i,P);
    B_prime_index = BRi[i];
    B_prime = BR[i];
    B->offsetRig(B_prime_index,  B_prime, P);
    i++;
  }  
  
  B->calculateStrokeDisplayList();
  B->calculateFillDisplayList();
  //visualize(A,A_prime,B,B_prime);  
}
*/

//------------------------------

bool DrawPath::correspondRoto2(RotoCurves* rc, const int w, const int h) { 

  if (rc->getNumCurves() == 0) return false;
  //if (_corrs || _corrsPtr) return false;

  DrawCorresponder* dcorer = new DrawCorresponder(this, rc, w, h);
  double score = dcorer->calculate();
  printf("correspondence score %f\n",score);
  //assert(!_corrs); 
  int* corrs = new int[getNumElements()];  
  RotoPath** corrsPtr = new RotoPath*[getNumElements()];
  dcorer->getSolution(corrs, corrsPtr);
  _corrRoto = new DrawContCorr(corrs,corrsPtr, getNumElements());
  delete[] corrs;
  delete[] corrsPtr;


  const set<RotoPath*>& rotos = getRotoSet();

  int numRotos = 0;
  for (set<RotoPath*>::const_iterator c = rotos.begin(); c != rotos.end(); ++c) {
    (*c)->addCorrDraw(this);
    ++numRotos;
  }
  printf("attached to %d rotoCurves in all\n",numRotos);
  //delete rotos;

  delete dcorer;
  return true;
}



void DrawPath::offsetRig(const float rT, RotoPath* corrRoto, const Vec3f P) {
  
  Vec2f U,V,blocf;
  blocf = corrRoto->getLoc(rT);  //corrRoto->getElement(index);
  corrRoto->getUV(rT,U,V);
  Vec2f_AddScale(blocf,blocf,U,P.x());
  Vec2f_AddScale(blocf,blocf,V,P.y());
  assert(finite(blocf.x()) && finite(blocf.y()));
  addVertex(blocf.x(),blocf.y(), P.z()); 
  _shouldBe->add(Vec2f(blocf.x(),blocf.y()));

  //_corrs[getNumElements()-1] = index;
  //_corrsPtr[getNumElements()-1] = corrRoto;
  
}

void DrawPath::getP(const int i, Vec3f& res) const { 
  Vec2f U,V;
  assert(_corrRoto);
  float rT = _corrRoto->getSampleRotoT(i);
  const RotoPath* rp = _corrRoto->getSampleRoto(i);
  rp->getUV(rT,U,V);
  Vec2f XO( getElement(i) , rp->getLoc(rT));

  res.set_x(XO.x()*U.x() + XO.y()*U.y());
  res.set_y(XO.x()*V.x() + XO.y()*V.y());
  res.set_z(_thick.getElement(i));
  assert(finite(res.x()) && finite(res.y()) && finite(res.z()));

}




void DrawPath::fillForwardInterpolatedCurve(DrawPath* B, const DrawPath* A) {
  B->initHertzStroke();
  float A_prime_t, B_prime_t;
  RotoPath *A_prime, *B_prime;
  Vec3f P;
  int i=0;
  assert(!B->getShouldBe());
  B->_shouldBe = new AbstractPath();
  
  B->_corrRoto = new DrawContCorr();
  A->_corrRoto->fillForward(B->_corrRoto);

  //const RotoPath* curr;
  while (i < A->getNumElements()) {

    A_prime = A->_corrRoto->getSampleRoto(i);
    A_prime_t = A->_corrRoto->getSampleRotoT(i);
    A->getP(i,P);
    B_prime_t = A_prime->getNextT(A_prime_t);  //A->_corrRoto->getSampleRotoT(i);
    B_prime = A_prime->nextC();  
    assert(B_prime);
    B->offsetRig(B_prime_t, B_prime, P);

    i++;
  }  
  
}

void DrawPath::fillBackwardInterpolatedCurve(DrawPath* B, const DrawPath* A) {
  B->initHertzStroke();
  float A_prime_t, B_prime_t;
  RotoPath *A_prime, *B_prime;
  Vec3f P;
  int i=0;
  assert(!B->getShouldBe());
  B->_shouldBe = new AbstractPath();
  
  B->_corrRoto = new DrawContCorr();
  A->_corrRoto->fillForward(B->_corrRoto);
  
  //const RotoPath* curr;
  while (i < A->getNumElements()) {
    
    A_prime = A->_corrRoto->getSampleRoto(i);
    A_prime_t = A->_corrRoto->getSampleRotoT(i);
    A->getP(i,P);
    B_prime_t = A_prime->getPrevT(A_prime_t);  //A->_corrRoto->getSampleRotoT(i);
    B_prime = A_prime->prevC();  
    assert(B_prime);
    B->offsetRig(B_prime_t, B_prime, P);
    
    i++;
  }  
  
}

// STUB
void DrawPath::fillBiInterpolatedCurve(DrawPath* D0, DrawPath* D1, DrawPath *Dnew, double t, 
				       int *D0_T0_T1_i, RotoPath** D0_T0_T1_p,
				       int *D0_T0_T_new_i, RotoPath** D0_T0_T_new_p, int *T1_D1) {
  /*Dnew->empty();
  Dnew->initHertzStroke();
  int numP = D0->getNumElements();
  Dnew->_corrs = new int[numP];
  Dnew->_corrsPtr = new RotoPath*[numP];
  assert(!Dnew->getShouldBe());
  Dnew->_shouldBe = new AbstractPath();
  
  int j;
  Vec3f P1, P2, P3;
  Vec2f XO, U, V;
  for (j=0; j<numP; j++) {
    D0->getP(j,P1);
    
    Vec2f_Sub(XO, D1->getElement(T1_D1[j]), 
	      D0_T0_T1_p[j]->getElement(D0_T0_T1_i[j]));
    D0_T0_T1_p[j]->getUV(D0_T0_T1_i[j], U, V);
    P2.set_x(XO.x()*U.x() + XO.y()*U.y());
    P2.set_y(XO.x()*V.x() + XO.y()*V.y());
    P2.set_z(D1->getThick(T1_D1[j]));

    Vec3f_Lerp(P3, P1, P2, t);
    Dnew->offsetRig(D0_T0_T_new_i[j], D0_T0_T_new_p[j], P3);
    D0_T0_T_new_p[j]->addCorrDraw(Dnew); // set should handle multiple inserts
  }
  
  Dnew->calculateStrokeDisplayList();
  Dnew->calculateFillDisplayList();
  */
}


/*
  NOT IMPORTANT
void DrawPath::interpolateForwards(int frame) {
  
  DrawPath *A=this;

  RotoPath *A_prime, *B_prime;
  // find earliest example this
  int df = 0, i=0;
  while (A->prevC() != NULL) {
    A = A->prevC();
    df++;
  }
  // find rotoCurves
  A_prime = A->getCorrespondRoto();
  if (!A_prime) return;
  B_prime = A_prime;
  ++df;              // now number of time steps from A' to newFrame B'
  printf("df: %d\n",df);
  for (i=0; i<df; i++) {
    B_prime = B_prime->nextC();
    if (!B_prime) return;
  }
  assert(B_prime->prevC() == getCorrespondRoto());

  printf("interpolating forward\n");
  
  // if later drawPath exists, clear it, otherwise create new one
  DrawPath* B;
  if ((B=nextC())) {
    B->empty();
    assert(B->prevC() == this);
  }
  else {
    B = new DrawPath();
    _is->getDrawCurves(frame)->addPath(B);
    setNextC(B);
    B->setPrevC(this);
  }

  // create interpolated drawPath
  B->copyLook(this);
  fillForwardInterpolatedCurve(B,A,B_prime,A_prime);
}
*/

/*
 PROBABLY DEAD
void DrawPath::fillForwardInterpolatedCurve(DrawPath* B, const DrawPath* A, RotoPath* B_prime,
			   const RotoPath* A_prime) {
  B->initHertzStroke();
  int A_prime_index, B_prime_index;
  Vec3f P;
  int i=0;
  
  B->rigFrom(B_prime, A->getNumElements());
  
  const RotoPath* curr;
  while (i < A->getNumElements()) {
    A_prime_index = A->getRotoCorr(i);
    A->getP(i,P);
    //B_prime_index = A_prime->getNextCorrIndex(A_prime_index);
    curr = A_prime; 
    B_prime_index = A_prime_index;
    do {
      B_prime_index = curr->getNextCorrIndex(B_prime_index);
      curr = curr->nextC();
      assert(curr);
    } while (curr != B_prime);
    B->offsetRig(B_prime_index, P);
    i++;
  }  
  
  B->calculateDisplayList();
  //visualize(A,A_prime,B,B_prime);
  
}
*/
/*
  PROBABLY DEAD
void DrawPath::fillBackwardInterpolatedCurve(DrawPath* B, const DrawPath* A, RotoPath* B_prime,
					      const RotoPath* A_prime) {
  B->initHertzStroke();
  int A_prime_index, B_prime_index;
  Vec3f P;
  int i=0;
  
  B->rigFrom(B_prime, A->getNumElements());
  
  const RotoPath* curr;
  while (i < A->getNumElements()) {
    A_prime_index = A->getRotoCorr(i);
    A->getP(i,P);
    //B_prime_index = A_prime->getNextCorrIndex(A_prime_index);
    curr = A_prime; 
    B_prime_index = A_prime_index;
    do {
      B_prime_index = curr->getPrevCorrIndex(B_prime_index);
      curr = curr->prevC();
      assert(curr);
    } while (curr != B_prime);
    B->offsetRig(B_prime_index, P);
    i++;
  }  
  
  B->calculateDisplayList();
  //visualize(A,A_prime,B,B_prime);
  
}
*/


void DrawPath::spliceIn(const DrawPath *other, const int i1, const int i2) {

  printf("DrawPAth splicing %d %d\n",i1,i2);
  int oldSize = getNumElements();
  int newsize = i1 + other->getNumElements() + oldSize-i2-1, i,j, k;
  ensureCapacity(newsize);
  _thick.ensureCapacity(newsize);
  count=newsize;
  _thick.setCount(newsize);

  Vec2f* oldPts = new Vec2f[oldSize];
  float* oldThick = new float[oldSize];
  memcpy(oldPts, data, oldSize*sizeof(Vec2f));
  memcpy(oldThick, _thick.getData(), oldSize*sizeof(float));

  /*
  for (i=newsize-1,j=oldSize-1; j>i2; i--,j--) {
    setElement(i,oldPts[j]);
    _thick.setElement(i,oldThick[j]);
  }
  for (i=i1,j=0; j<other->getNumElements(); i++,j++) {
    setElement(i,other->getElement(j));
    _thick.setElement(i, other->_thick.getElement(j));
    }
  */
  float fac, blT;
  int spliceL = other->getNumElements();

  for (i=i1, j=0; j<spliceL; i++, j++) {
    setElement(i,other->getElement(j));
    if (j<9 && (k=i)<oldSize) {
      fac = float(j)/8.f;
      blT = oldThick[k] + fac*(other->_thick.getElement(j) - oldThick[k]);
      _thick.setElement(i, blT);
    }
    else if (j>=spliceL-9 && (k= i2 - (spliceL-j-1) ) < oldSize) {
      fac = float(spliceL-j-1)/8.f; assert(fac>=0);
      printf("end %f\n",fac);
      blT = oldThick[k] + fac*(other->_thick.getElement(j) - oldThick[k]);
      _thick.setElement(i, blT);
    }
    else
      _thick.setElement(i, other->_thick.getElement(j));
  }
  
  for (k=i2+1; i<newsize; i++, k++) {
    setElement(i,oldPts[k]);
    _thick.setElement(i,oldThick[k]);
  }
  
  //i: new curve, j: spliceIn, k: old curve
  /*
  for (i=i1-lB,m=0,j=0, k=i1-lB; m<2*lB; m++,i++,j++,k++) {
    fac = float(m)/float(lB);
    Vec2f_LinInterp(blend, oldPts[k], other->getElement(j), fac);
    blT = oldThick[k] + fac*(other->_thick.getElement(j) - oldThick[k]);
    setElement(i,blend);
    _thick.setElement(i, blT);
  }
  */

  /*  for (m=0, k=i2+1-lB; m<hB; m++,i++,j++,k++) {
    fac = float(m)/float(lB);
    Vec2f_LinInterp(blend, oldPts[k], other->getElement(j), fac);
    setElement(i,blend);
    blT = oldThick[j] + fac*(other->_thick.getElement(j) - oldThick[k]);
    _thick.setElement(i, blT);
    }*/

 
  
  


  delete[] oldPts; delete[] oldThick;

  //for (i=0; i<newsize; i++)
  //printf("%f %f\n",getElement(i).x(), getElement(i).y());
  

  delete _corrRoto; _corrRoto = NULL;
  //delete[] _corrs; _corrs = NULL;
  //delete[] _corrsPtr; _corrsPtr = NULL;

  redoStroke();
  calculateStrokeDisplayList();
  calculateFillDisplayList();
}


void DrawPath::fixLoc() {
  AbstractPath::fixLoc();
  redoStroke();
  calculateStrokeDisplayList();
  calculateFillDisplayList();
}


void DrawPath::setFixed(const bool s) { 
  if (!_fixed && _shouldBe && s==true) {
    delete _shouldBe;
    _shouldBe=false;
  }
  _fixed = s; 
}

void DrawPath::nudge(const int i, const Vec2f& loc) {
  AbstractPath::nudge(i, loc);
}


// useless
void DrawPath::repToShouldBe() {
  if (_shouldBe || _fixed) return;
  _shouldBe = new AbstractPath();
  _shouldBe->addNElements(getNumElements());
  for (int j=0; j<getNumElements(); j++)
    _shouldBe->setElement(j,getElement(j));
}

float DrawPath::_dAlpha;


void DrawPath::setCurrTexture ( int i ) { 

  //printf("set current texture\n") ; 


  if ( i < 0 || i > NUM_TEXTURES ) { 
    printf("i = %d\n", i ) ; 
    assert ( i >= 0 && i <= NUM_TEXTURES ) ; 
  } 



  // if the chosen color is not loaded, load it
  if ( DrawPath::_texture_names[i] == 0 ) {

    printf("LOADING TEXTURE\n"); 
    GLuint new_texture ; 
    
    if ( i == PEN ) { 

      new_texture = 999; 

    }else{ 

      glGenTextures(1, &new_texture);
      Texture* canvas = new Texture (PNG, new_texture, Texture::Texture_FNames[i], 
				       Texture::Texture_Widths[i], Texture::Texture_Heights[i]) ;      
      delete canvas ; 
    }

    DrawPath::_texture_names[i] = new_texture; 
  }
  
  if ( i == PEN ) { 

      DrawPath::_nxt_textureWidth = 0 ; 
      DrawPath::_nxt_textureHeight = 0 ; 
      DrawPath::_nxt_textureType = none ; 
      DrawPath::_nxt_useTexture = false; 

  } else { 

    DrawPath::_nxt_textureWidth = Texture::Texture_Widths[i]; 
    DrawPath::_nxt_textureHeight = Texture::Texture_Heights[i];  
    DrawPath::_nxt_useTexture = true; 
    DrawPath::_nxt_textureType = Texture::Texture_Types[i]; 
  } 

  DrawPath::_nxt_textureRadius = Texture::Texture_Radii[i]; 

  // bind to the current texture

  DrawPath::_nxt_texture_index = i ; 
  DrawPath::_nxt_texture_name = DrawPath::_texture_names[i] ; 
} 


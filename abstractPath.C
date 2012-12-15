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
#include "abstractPath.h"
//#include "GraphicsGems.h"


AbstractPath::AbstractPath() {
  _tangents = NULL;
  _inMotion = false;
  //_centroidCalculated = false;
}

AbstractPath::AbstractPath(const AbstractPath& other) : DynArray<Vec2f,100>(other) {
  _tangents = NULL;
  _inMotion = false;
  //_centroidCalculated = false;
}

AbstractPath::~AbstractPath() {
  if (_tangents) 
    delete[] _tangents;
}


#define PHI 0.25f
#define MU -0.2564f

void AbstractPath::smooth() {
  int ul = getNumElements()-1;
  Vec2f offset, newVal;
  for (int i=1; i<ul; i++) {
    filterOffset(i,offset);
    offset *= PHI;
    newVal = getElement(i);
    newVal += offset;
    setElement(i,newVal);
  }
}

void AbstractPath::sharpen() {
  int ul = getNumElements()-1;
  Vec2f offset, newVal;
  for (int i=1; i<ul; i++) {
    filterOffset(i,offset);
    offset *= MU;
    newVal = getElement(i);
    newVal += offset;
    setElement(i,newVal);
  }
}

inline void AbstractPath::filterOffset(const int i, Vec2f& here) const {
  Vec2f P1(getElement(i-1), getElement(i));
  Vec2f P2(getElement(i+1), getElement(i));
  
  float p1len = P1.Len(), p2len = P2.Len(), sum = p1len+p2len;
  P1 *= p2len / sum;
  P2 *= p1len / sum;

  here = P1;
  here += P2;
}

// CORRECTNESS: not calculating tangents across closed rotoPaths
Vec2f AbstractPath::tangent(const int i) const { 
  assert(i>=0 && i<getNumElements());
  if (_tangents)
    return _tangents[i];
  Vec2f U;
  if (i==0)
    Vec2f_Sub(U,getElement(1), getElement(0));
  else if (i==getNumElements()-1) 
    Vec2f_Sub(U,getElement(getNumElements()-1), 
	      getElement(getNumElements()-2));
  else {
    Vec2f_Sub(U,getElement(i+1), getElement(i-1));
  }
  U.Normalize();
  return U;
}

void AbstractPath::renderPS(FILE* fp) const {
  Vec2f location;
  for (int i=0; i<getNumElements(); i++) {
    location = getElement(i);
    fprintf(fp,"%f %f .5 .5 rectfill\n",location.x()-.25, 
	    _globalh - location.y()-.25);
  }

}


double AbstractPath::distanceToEnd(const Vec2f pt) const {
  return getPointerToElement(getNumElements()-1)->distanceTo2(pt);
}
double AbstractPath::distanceToStart(const Vec2f pt) const {
  return getPointerToElement(0)->distanceTo2(pt);
}


double AbstractPath::distanceTo2(const Vec2f pt, int& index) const { 
  int i, max = getNumElements();
  float res=FLT_MAX,tmp;
  for (i=0; i<max; i++) {
    if ((tmp=getPointerToElement(i)->distanceTo2(pt)) < res) {
      res = tmp;
      index = i;
    }
  }
  return ((double)res);
}

double AbstractPath::distanceTo2(const Vec2f pt) const { 
  int shit;
  return distanceTo2(pt,shit);
}

void AbstractPath::writePtr(FILE* fp) const {
  fwrite(_ptrDoc,sizeof(int),2,fp);
}

void AbstractPath::writePtr(QDataStream* fp) const {
  *fp << _ptrDoc[0] << _ptrDoc[1];
}

void AbstractPath::writeNullPtr(FILE* fp) {
  int dummy[2] = {-1,-1};
  fwrite(&dummy,sizeof(int),2,fp);
}

void AbstractPath::writeNullPtr(QDataStream* fp) {
  int dummy[2] = {-1,-1};
  *fp << dummy[0] << dummy[1];
}

void AbstractPath::writeIntArray(FILE* fp, int* array, int size) {
  assert(array && size>0);
  fwrite(&size,sizeof(int),1,fp);
  fwrite(array,sizeof(int),size,fp);
}
void AbstractPath::writeIntArray(QDataStream* fp, int* array, int size) {
  assert(array && size>0);
  *fp << size;
  for (int i=0; i<size; ++i)
    *fp << array[i];
}

void AbstractPath::readIntArray(QDataStream* fp, int* array, int size) {
  if (array != NULL) {
    for (int i=0; i<size; ++i)
      *fp >> array[i];
  }
  else {
    int j;
    for (int i=0; i<size; ++i)
      *fp >> j;
  }
}

void AbstractPath::writeBool(QDataStream *fp, bool b) {
  int i = b;
  *fp << i;
}

void AbstractPath::readBool(QDataStream *fp, bool& b) {
  int i;
  *fp >> i;
  assert(i==0 || i==1);
  b = (bool)i;
}

void AbstractPath::print() const {
  for (int i=0; i<getNumElements(); i++) {
    Vec2f pt = getElement(i);
    printf("%.4f %.4f\n",pt.x(),pt.y());
  }
  printf("\n");

}

#define _SSPACING_ 2.
void AbstractPath::resample(const int* numSamples, DynArray<float,100>* interpMe) {
  int oldCount = getNumElements();
  Vec2f* oldSamples = new Vec2f[oldCount];
  float* oldVals = new float[oldCount];
  memcpy(oldSamples, getData(), oldCount*sizeof(Vec2f));
  memcpy(oldVals, interpMe->getData(), oldCount*sizeof(float));
  resetElements();
  interpMe->resetElements();
  
  float m=0, totSegLen=0, L=0, n, s, segLen, st=0, rSegLen, tmpF;
  int i;
  Vec2f tmp;
  for (i=1; i<oldCount; i++) {
    Vec2f_Sub(tmp,oldSamples[i], oldSamples[i-1]);
    L +=tmp.Len();
  }

  add(oldSamples[0]);
  interpMe->add(oldVals[0]);
  
  if (numSamples==NULL)
    n = ceil(L/_SSPACING_);
  else
    n = float(*numSamples-1);
  s = L / n; 	      assert(finite(s)); 
  for (i=1; i<oldCount; i++) {
    Vec2f_Sub(tmp,oldSamples[i], oldSamples[i-1]);
    segLen = tmp.Len();
    totSegLen += segLen;
    if (s-st <= segLen) {
      m = (s-st)/segLen; //assert(m>=0 && m<=1);

      Vec2f_LinInterp(tmp,oldSamples[i-1], oldSamples[i],m);
      add(tmp);
      tmpF = oldVals[i-1] + m*(oldVals[i]-oldVals[i-1]);
      interpMe->add(tmpF);

      rSegLen = segLen;
      segLen -= (s-st);

      while (segLen>s) {
	m += s/rSegLen; //assert(m>=0 && m<=1);
	//printf("in here\n");
	Vec2f_LinInterp(tmp,oldSamples[i-1], oldSamples[i],m);
	add(tmp);
	tmpF = oldVals[i-1] + m*(oldVals[i]-oldVals[i-1]);
	interpMe->add(tmpF);

	segLen -= s;
      } // endwhile
      st = segLen;
    } 
    else
      st += segLen;
  } // endfor

  add(oldSamples[oldCount-1]);
  interpMe->add(oldVals[oldCount-1]);
  
}


/*
void AbstractPath::centroidBbox(Vec2f* cent, Vec2f* dims) const {
  cent->Set(0,0);
  Vec2f curr;
  Bboxf2D bbox;
  for (int i=0; i<getNumElements(); i++) {
    curr = getElement(i);
    *cent += curr;
    bbox.includePoint(curr);
  }
  *cent /= float(getNumElements());

  _centroid = *cent;
  _centroidCalculated = true;

  dims->Set(bbox.width(), bbox.height()); // not completely right
  // also, centroid may not be on curve
}

void AbstractPath::getCentroid(Vec2f& here) const {
  if (!_centroidCalculated) {
    here.Set(0,0);
    Vec2f curr;
    for (int i=0; i<getNumElements(); i++) {
      curr = getElement(i);
      here += curr;
    }
    here /= float(getNumElements());
    _centroid = here;
    _centroidCalculated = true;
  }
  else 
    here = _centroid;
  
}
*/

void AbstractPath::translate(const Vec2f& delta) {
  for (int i=0; i<getNumElements(); i++)
    (*getPointerToElement(i)) += delta;
  _inMotion = true;
  _motion += delta;
}

void AbstractPath::fixLoc() {
  _motion.Set(0,0);
  _inMotion = false;
}


void AbstractPath::renderDirection() const {
  glBegin(GL_LINES);
  for (int i=0; i<getNumElements(); i+=15) {
    Vec2f tan = tangent(i);
    Vec2f perp(-tan.y(), tan.x());
    Vec2f ourLoc = getElement(i);
    perp *= 5.f;
    tan *= 5.f;
    Vec2f left(ourLoc, perp);
    left -= tan;
    glVertex2f(left.x(), left.y());
    glVertex2f(ourLoc.x(), ourLoc.y());
    glVertex2f(ourLoc.x(), ourLoc.y());
    ourLoc += perp;
    ourLoc -= tan;
    glVertex2f(ourLoc.x(), ourLoc.y());
  }
  glEnd();
}

void AbstractPath::nudge(const int i, const Vec2f& loc) {
  printf("Nudge %d, %f %f\n",i,loc.x(), loc.y());
  //int low = MIN(0,i-_nudgeRadius), high = MAX(size-1
}


void AbstractPath::modifyEnd(const int which, const Vec2f newloc) {
  if (which==0)
    setElement(0, newloc);
  else if (which==1)
    setElement(getNumElements()-1, newloc);
  else
    assert(0);
}

int AbstractPath::_globalh;
ImgSequence* AbstractPath::_is;
TalkFitCurve AbstractPath::_talkFit;
int AbstractPath::_nudgeRadius;


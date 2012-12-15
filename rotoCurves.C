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



#include <stdio.h>
#include <float.h>
#include "rotoCurves.h"

RotoCurves::RotoCurves() {
  //_active = NULL;
  
}

/*
void RotoCurves::save(const char* fileName) const {
  int numPaths = _paths.getSize();
  if (numPaths == 0) return;

  FILE* fp = fopen(fileName,"w");
  assert(fp);
  
  fwrite(&numPaths,sizeof(int),1,fp);
  int p;
  const RotoPath* curr;
  _paths.SetIterationHead();
  for (p=0; p<numPaths; p++) {
    curr = _paths.IterateNext();
    int numPts = curr->getNumElements();
    fwrite(&numPts,sizeof(int),1,fp);
    fwrite(curr->getData(),sizeof(Vec2f),numPts,fp);
  }
  fclose(fp);
}
*/

void RotoCurves::save(FILE* fp, int version) const { 
  
  /*_paths.SetIterationHead();
  const RotoPath* curr;
  //while ((curr = _paths.IterateNext()) != NULL) */
    
  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) 
    (*c)->save(fp, version);  
}


void RotoCurves::saveqt(QDataStream* fp, int version) { 
  checkAllJoints(); 
  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c)
    (*c)->saveqt(fp, version);  
}

void RotoCurves::saveMatte(int w, int h) const {
  if (_regions.empty()) return;
  QImage im(w,h,32);
  im.fill(qRgb(255,255,255));
  
  std::deque<RotoRegion*>::const_iterator c;
  for (c=_regions.begin(); c!=_regions.end(); ++c)
    (*c)->renderMask(im.bits());

  char name[200];
  sprintf(name, "matte%.3d.png",_frame);
  im.save(name, "PNG");
}

void RotoCurves::saveRotoXML(std::ofstream& fp, int frame) const {
  fp << "<RotoFrame: " << frame << std::endl;
  fp << _paths.size() << std::endl;
  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c)
    (*c)->saveXML(fp);  
  fp << ">" << std::endl;
}


void RotoCurves::clearXMLLabels() {
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c)
    (*c)->clearXMLLabels();
}

void RotoCurves::createXMLLabels(int& labelNum) {

  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c)
    (*c)->labelSelfXML(labelNum);
}


/*
void RotoCurves::load(const char* fileName) {
  FILE* fp = fopen(fileName,"r");
  if (!fp) return;

  int numPaths,numPts,p, res;
  res = fread(&numPaths,sizeof(int),1,fp); assert(res==1);
  for (p=0; p<numPaths; p++) {
    printf("loading %s %d\n",fileName,p);
    RotoPath* newPath = new RotoPath();
    res = fread(&numPts,sizeof(int),1,fp);
    newPath->ensureCapacity(numPts);
    res = fread(newPath->getData(),sizeof(Vec2f),numPts,fp);
    newPath->setCount(numPts);
    assert(res==numPts);
    _paths.AddToTail(newPath);
  }

  fclose(fp);
  
}

*/


void RotoCurves::load(FILE* fp, int version) { 
  // get data in
  /*
  RotoPath* curr;
  LList<RotoPath> copy;
  copy.copyPointers(&_paths); // Hack to allow nested loops
  copy.SetIterationHead();
  while ((curr = copy.IterateNext()) != NULL) 
    curr->load(fp, version);
    copy.forgetContents();*/

  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c)
    (*c)->load(fp, version);
}

void RotoCurves::loadqt(QDataStream* fp, int version) { 
  // get data in
  /*
  RotoPath* curr;
  LList<RotoPath> copy;
  copy.copyPointers(&_paths); // Hack to allow nested loops
  copy.SetIterationHead();
  while ((curr = copy.IterateNext()) != NULL) */


  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c)
    (*c)->loadqt(fp, version);
 
  //copy.forgetContents();
  checkAllJoints();
}

void RotoCurves::checkAllJoints() {
  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) 
    assert((*c)->allJointsOk());
  //printf("checked joints\n");
}



void RotoCurves::deletePath(RotoPath* rp) {
  rp->vacateRotoCorr();

  const set<DrawPath*>* corrDraws = rp->getCorrDraws();
  for (set<DrawPath*>::const_iterator c = corrDraws->begin(); c!=corrDraws->end(); ++c)
    (*c)->vacateRotoCorr();

  rp->destroyJoints();

  //for (int i=0; i < rp->getNumCorrDraws(); i++) {
  //DrawPath* dp = rp->getCorrDraw(i);
  //dp->vacateRotoCorr();
    //dp->setCorrNum(-1);
    //dp->setCorrespondRoto(NULL);
  //}

  //if (rp == _active)
  //_active = NULL;
  //_paths.RemoveNode(rp,1);  
  _paths.remove(rp);

  for (short s=0; s<2; ++s) {
    if (rp->getRegion(s)) {
      std::deque<RotoRegion*>::iterator c = std::remove(_regions.begin(), _regions.end(), rp->getRegion(s)), c2;
      c2=c;
      for ( ; c!= _regions.end(); ++c)
	delete (*c);
      _regions.erase(c2, _regions.end());  
    }
  }
  

  delete rp;
}


void RotoCurves::documentSelf() { 
  int i=0;

  /*_paths.SetIterationHead();
  RotoPath* curr;
  while ((curr=_paths.IterateNext()) != NULL) {*/

  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    (*c)->setPtrDoc(_frame,i);
    i++;
  }
}


void RotoCurves::addNPaths(const int n) {
  for (int i=0; i<n; i++)
    _paths.push_back(new RotoPath());
    //_paths.AddToTail(new RotoPath());
}

RotoPath* RotoCurves::getCurve(const int i) {
  assert(i>=0 && i<(int) _paths.size());
  RotoPathList::iterator c = _paths.begin();
  for (int j=0; j<i; ++j, ++c) {}
  return *c;
    
  //return (  *(_paths.begin()+i) );
  /*
  _paths.SetIterationHead();
  RotoPath* res = _paths.IterateForward(i);
  assert(res);
  return res;*/
}



void RotoCurves::renderAllPaths(bool showTrackPoints, bool visEffort) {
  //startCurvesIterator();
  //const RotoPath* curr;
  glColor3f(1,0,1);
  glPointSize(2.0);
  int name=0;
  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    //while ((curr=getNextCurve()) != NULL) {
    glPushName(name++);
    if (!(*c)->fixed())//G!
      glColor3f(1,1,0);
    else
	glColor3f(1,0,0);
    //glColor3f(.5,.5,.5);
    (*c)->render(showTrackPoints, visEffort);
    glPopName();
  }
}

void RotoCurves::renderNextCorrLines() {
  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    if (!(*c)->nextC()) continue;
    (*c)->renderNextCorrLines();
  }
}

void RotoCurves::renderFreePaths(bool showTrackPoints) {
  //startCurvesIterator();
  //const RotoPath* curr;
  int name=0;

  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    //while ((curr=getNextCurve()) != NULL) {
    glPushName(name++);
    
    if (!(*c)->nextC()) {
      (*c)->render(showTrackPoints);
    }

    glPopName();
  }
}


RotoPath* RotoCurves::pickPrimitive(const int x, const int y) {
  //printf("picking\n");
  GLint viewport[4];
  GLdouble projection[16];
  GLuint selectBuf[100];
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glSelectBuffer(100, selectBuf);
  glRenderMode(GL_SELECT); 
  glInitNames();
  
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPickMatrix(x, AbstractPath::_globalh-y, 10., 10., viewport);
  glMultMatrixd(projection);
  glMatrixMode(GL_MODELVIEW);

  renderAllPaths(false);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glFlush();

  int hits = glRenderMode(GL_RENDER);
  //printf("%d hits\n",hits);
  if (hits==0) return NULL;
  int maxHit = selectBuf[3];
  for (int i=1; i<hits; i++) {
    if (selectBuf[3+4*i] > (unsigned) maxHit) {
      maxHit = selectBuf[3+4*i];
    }
  }
  
  return getCurve(maxHit);
  //_paths.SetIterationHead();
  //printf("selected path\n");
  //return _paths.IterateForward(maxHit);
}

RotoPath* RotoCurves::pickFreePrimitive(const int x, const int y) {
  printf("picking\n");
  GLint viewport[4];
  GLdouble projection[16];
  GLuint selectBuf[100];
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glSelectBuffer(100, selectBuf);
  glRenderMode(GL_SELECT); 
  glInitNames();
  
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPickMatrix(x, AbstractPath::_globalh-y, 10., 10., viewport);
  glMultMatrixd(projection);
  glMatrixMode(GL_MODELVIEW);

  renderFreePaths(false);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glFlush();

  int hits = glRenderMode(GL_RENDER);
  printf("%d hits\n",hits);
  if (hits==0) return NULL;
  unsigned int maxHit = selectBuf[3];
  for (int i=1; i<hits; i++) {
    if (selectBuf[3+4*i] > maxHit) {
      maxHit = selectBuf[3+4*i];
    }
  }
  
  //_paths.SetIterationHead();
  //printf("selected path\n");
  //return _paths.IterateForward(maxHit);
  return getCurve(maxHit);
}

Vec2i RotoCurves::calcStat() {
  Vec2i res;
  //RotoPath* curr;
  //_paths.SetIterationHead();
  //while ((curr=_paths.IterateNext()) != NULL) {

  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    res += (*c)->getStats();
    /*int n = (*c)->getNumControls();
    for (int i=0; i<n; ++i) {
      res.Inc(0,1);
      if ((*c)->touched(i)) 
	res.Inc(1,0);
	}*/
  }
  return res;
}


RotoPath* RotoCurves::cycle(RotoPath* p) {
  RotoPathList::iterator c = find(_paths.begin(), _paths.end(), p);
  RotoPathList::iterator last = _paths.end();
  --last;
  if (c == _paths.end() || c == last)
    return _paths.front();
  else {
    ++c;
    return *c;
  }
  /*if (p->next) 
    return p->next;
  else
  return _paths.front();*/
}



void RotoCurves::assertNoDraws() {
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    //while ((curr=_paths.IterateNext()) != NULL) {
    //assert(curr->getNumCorrDraws() == 0);
    (*c)->emptyCorrDraws();
  }
  printf("assertion success\n");
}


void RotoCurves::printJoints() {
  int i=0;
  RotoPathList::const_iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    printf("%d: \n", i++);
    RotoPath::printJoints((*c)->bJoints());
    RotoPath::printJoints((*c)->eJoints());
  }
}

#define _J_DIST_ 25
void RotoCurves::setupJoints(RotoPath* rp) {
  //setupBegJoints(rp);
  //setupEndJoints(rp);
  setupJointsOneSide(rp,0);
  setupJointsOneSide(rp,1);
  makeJointsConsistent(rp->joints(0));
  makeJointsConsistent(rp->joints(1));
  checkAllJoints();
}

void RotoCurves::setupJointsOneSide(RotoPath* rp, int side) {

  double tmp;
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    if ((*c)==rp) continue; // CORRECTNESS: allow looping curves
    if ((tmp=(*c)->distanceToSide2(rp->getEnd(side), 0)) < _J_DIST_) {
      printf("Joining %d side to beginning of other curve\n",side);
      makeJoint(rp,(*c),side,0);
    }
    if ((tmp=(*c)->distanceToSide2(rp->getEnd(side), 1)) < _J_DIST_)  {
      printf("Joining %d side to end of other curve\n",side);
      makeJoint(rp,(*c),side,1);
    }
    assert((*c)->okEnds());
  }
}

void makeJointsConsistent(Joints& jnts) {
  Joints::iterator c;
  
  // build a list of all rotopaths, sides involved.
  // make sure each rotopath knows about all the others
  JointSet all; 
  for (c=jnts.begin(); c!=jnts.end(); ++c)
    c->_rp->addToJointSet(all, c->_side);
  

  JointSet::iterator c1,c2;
  for (c1=all.begin(); c1!=all.end(); ++c1) {
    for (c2=all.begin(); c2!=all.end(); ++c2) {
      if (c1->_rp == c2->_rp) continue; // no joint to self

      // duplicates taken care of by insertJoint
      c1->_rp->insertJoint(*c2, c1->_side); 

    }
  }
}

void RotoCurves::makeJoint(RotoPath* rp1, RotoPath* rp2, int w1, int w2) {
  Joint j1, j2;
  j1._rp = rp2; j1._side = w2;
  j2._rp = rp1; j2._side = w1;
  rp1->insertJoint(j1,w1);
  rp2->insertJoint(j2,w2);

  //Vec2f loc1 = rp1->getElement(w1 * (rp1->getNumElements()-1));
  //Vec2f loc;
  //Vec2f_Average(loc, loc1, loc2);

  //Vec2f loc2 = rp2->getElement(w2 * (rp2->getNumElements()-1));

  rp1->setEnd(rp2->getEnd(w2), w1);

  /*
  int rp1i, rp1i1;
  if (w1==0) {
    rp1i = 0; rp1i1 = 1;
  }
  else {
    rp1i = rp1->getNumElements()-1;
    rp1i1 = rp1i - 1;
  }

  rp1->setElement(rp1i, loc2);
  if (rp1->getElement(rp1i) == rp1->getElement(rp1i1))
    rp1->removeElement(rp1i1);
  
    assert(rp1->okEnds());*/
}


RotoPath* RotoCurves::distanceToEnds2(const Vec2f& loc, int& whichEnd) {  

  double minDist = DBL_MAX, tmp;
  RotoPath* res; int resSide;
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    //while ((curr=getNextCurve()) != NULL) {
    tmp = (*c)->distanceToSide2(loc,1);
    if (tmp < minDist) {
      res = (*c); resSide = 1; minDist = tmp;
    }
    tmp = (*c)->distanceToSide2(loc, 0);
    if (tmp < minDist) {
      res = (*c); resSide = 0; minDist = tmp;
    }
  }

  if (minDist < 10) {
    whichEnd = resSide;
    return res;
  }
  else 
    return NULL;
}


void RotoCurves::setupBegJoints(RotoPath* rp) {
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    //startCurvesIterator();
    //RotoPath* curr;
    //while (((*c)=getNextCurve()) != NULL) {
    if ((*c)==rp) continue; // CORRECTNESS: allow looping curves
    if ((*c)->distanceToStart(rp->getElement(0)) < _J_DIST_)
      makeJoint(rp,(*c),0,0);
    if ((*c)->distanceToEnd(rp->getElement(0)) < _J_DIST_)
      makeJoint(rp,(*c),0,1);
    assert((*c)->okEnds());
  }
}

void RotoCurves::setupEndJoints(RotoPath* rp) {
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    //  while ((curr=getNextCurve()) != NULL) {
    if ((*c)==rp) continue; // CORRECTNESS: allow looping curves
    if ((*c)->distanceToStart(rp->getElement(rp->getNumElements()-1)) < _J_DIST_)
      makeJoint(rp,(*c),1,0);
    if ((*c)->distanceToEnd(rp->getElement(rp->getNumElements()-1)) < _J_DIST_)
      makeJoint(rp,(*c),1,1);
    assert((*c)->okEnds());
  }

}




RotoPath* RotoCurves::distanceToCtrls2(const Vec2f& loc, int* whichControl) {
  double minDist = DBL_MAX, tmp;
  RotoPath* res; 
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    //  while ((curr=getNextCurve()) != NULL) {
    const std::vector<Vec2f>& ctrls = (*c)->getBez()->getControls();
    std::vector<Vec2f>::const_iterator i = ctrls.begin();
    for (; i!= ctrls.end(); ++i) {
      tmp = i->distanceTo2(loc);
      if (tmp<minDist) {
	minDist = tmp;
	res = (*c);
	*whichControl = i - ctrls.begin(); // is this right?
      }      
    }
  }

  if (minDist < 10)
    return res;
  else
    return NULL;
}

/*
void RotoCurves::reconcileJoints() {
  startCurvesIterator();
  RotoPath* curr;
  while ((curr = getNextCurve()) != NULL) {
    reconcileJoint(curr->bJoints(), curr->getPointerToElement(0));
    reconcileJoint(curr->eJoints(), curr->getPointerToElement(curr->getNumElements()-1));
  }
}
*/

/*
RotoPath* RotoCurves::getNextCurve() {
  return _paths.IterateNext();
}


const RotoPath* RotoCurves::getNextCurve() const {
  return _paths.IterateNext();
}
*/
void RotoCurves::addPath(RotoPath* newPath) {
  _paths.push_back(newPath);
}

void RotoCurves::deleteTail() { 
  delete *(_paths.end());
  _paths.pop_back(); 
}


void RotoCurves::splitCurve(RotoPath* rp, const int ctrl) {
  

  // remove rp
  //_paths.RemoveNode(rp, 0); 
  _paths.remove(rp);
  // TODO: should elegantly handle roto correspondences (split them),
  // and or propagate change to other moments in time


  // make into two
  RotoPathPair two = rp->split(ctrl);
  
  // handle joints
  //rp->replaceMeInJoints(two.first, 0);
  //rp->replaceMeInJoints(two.second, 1);


  // add them
  addPath(two.first);
  addPath(two.second);
  delete rp;

  checkAllJoints(); 
}

RotoPath* RotoCurves::pickSplineSegment(const Vec2f loc, float& t) {
  // STUB
  double minDist = DBL_MAX, dist;
  RotoPath* res; 
  RotoPathList::const_iterator c;
  DSample sample;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    dist = (*c)->distToSplineSamples2(loc, sample);
    if (dist < minDist) {
      minDist = dist;
      t = sample.totT();
      res = *c;
    }
  }
  
  //printf("pickSplineSegment: minDist %f, t %f, RotoPath %x\n",minDist,t,res);

  if (minDist < 25)
    return res;
  else
    return NULL;
}

void RotoCurves::addRotoRegion(RotoRegion* r) {
  if (find(_regions.begin(), _regions.end(), r) == _regions.end())
    _regions.push_back(r);
}

RotoRegion* RotoCurves::pickRegion(const float x, const float y) {
  std::deque<RotoRegion*>::reverse_iterator c;
  for (c=_regions.rbegin(); c!=_regions.rend(); ++c) {
    if ((*c)->withinRegion(x,y))
      return *c;
  }

  return NULL;
}

void RotoCurves::renderRegions() {
  std::deque<RotoRegion*>::const_iterator c;
  for (c=_regions.begin(); c!=_regions.end(); ++c) {
    (*c)->render();
  }
}


// makes mask by finding lowest region in deck that these paths are
// part of, and then AND'ing together the masks of regions above
// it which would occlude them.
// Returns null if no such occluding regions
unsigned char* RotoCurves::makeCummMask(const PathV& paths) {
    std::deque<RotoRegion*>::const_iterator c;
    PathV::const_iterator c2;
    bool found=false;
    for (c=_regions.begin(); c!=_regions.end() && !found; ++c) {
      for (c2 = paths.begin(); c2!=paths.end() && !found; ++c2) {
	if ( ((*c2)->getRegion(RR_RIGHT) == *c) ||
	     ((*c2)->getRegion(RR_LEFT) == *c))
	  found = true;
      }
    }

    if (!found)
      c = _regions.begin(); // un-associated curves are below rest of deck
    else
      ++c; // advance to first possible occluding contour

    if (c == _regions.end()) // no occluding regions
      return NULL;

    printf("Occluding region found\n");
    unsigned char* res = (*c)->copyMask();
    ++c;
    for ( ; c!=_regions.end(); ++c)
      (*c)->orMask(res);
    return res;
}

void RotoCurves::calcLerp() {
  RotoPathList::iterator c;
  for (c=_paths.begin(); c!= _paths.end(); ++c) {
    (*c)->calcLerp();
  }
}

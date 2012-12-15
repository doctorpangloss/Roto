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



#include "imgSequence.h"

ImgSequence::ImgSequence(char* imageRoot,int startF, int endF) : 
  _startF(startF), _endF(endF) {
  strcpy(_imageRoot, imageRoot);
  //strcpy(_fileRoot, fileRoot);
  strcpy(_fileRoot,"");
  mutualInit();
}

void ImgSequence::mutualInit() {
  sprintf(_fileName,"%s%.3d.png",_imageRoot,_startF);
  QImage nothing;
  printf("%s\n",_fileName);
  printf("%s\n",_imageRoot);
  int res = nothing.load(_fileName);
  if (!res) {
    printf("Failed to find image\n");
    exit(0);
  }
  _w = nothing.width();
  _h = nothing.height();

  /*** remove ****/ 
  printf("WIDTH = %d\n", _w ) ;
  printf("HEIGHT = %d\n", _h) ; 

  AbstractPath::_globalh = _h;
  
  _rotoc = new RotoCurves[_endF-_startF+1];
  _drawc = new DrawCurves[_endF-_startF+1];
  _pyrms = new KLT_FullPyramid[_endF-_startF+1];
  _Cpyrms = new KLT_FullCPyramid[_endF-_startF+1];
  _pyrmsE = new KLT_FullPyramid[_endF-_startF+1];
  for (int i= _startF; i<=_endF; i++) {
    _rotoc[i-_startF].setFrame(i);
    _drawc[i-_startF].setFrame(i);
  }
}

bool ImgSequence::fileRootInited() { 
  return (strcmp(_fileRoot, "")!=0);
}

QImage ImgSequence::getImage(const int a) const {
  assert(a >= _startF && a <= _endF);
  sprintf(_fileName,"%s%.3d.png",_imageRoot,a);
  QImage px(_fileName);
  return px;
}


ImageGraph* ImgSequence::getImageGraph(const int a) const {
  assert(a >= _startF && a <= _endF);
  sprintf(_fileName,"%s%.3d.dat",_imageRoot,a);
  FILE *fp = fopen(_fileName,"r");
  if (!fp) {
    printf("Need to run imageChewer first\n");
    exit(0);
    //return NULL;
  }
  ImageGraph* ig = new ImageGraph(fp);
  //QImage im = getImage(a);
  //ImageGraph* ig = new ImageGraph(im);
  return ig;
}


QImage ImgSequence::getRelativeImage(const int a) const {
  return getImage(_startF+a);
}

ImageGraph* ImgSequence::getRelativeImageGraph(const int a) const {
  return getImageGraph(_startF+a);
}

RotoCurves* ImgSequence::getRotoCurves(const int a) const {
  if (!(a >= _startF && a <= _endF))
    return NULL;
  return (_rotoc+a-_startF);
}

RotoCurves* ImgSequence::getRelativeRotoCurves(const int a) const {
  if (!(a >= 0 && a <= _endF-_startF))
      return NULL;
  return (_rotoc+a);
}

DrawCurves* ImgSequence::getDrawCurves(const int a) const {
  if (!(a >= _startF && a <= _endF))
    return NULL;
  return (_drawc+a-_startF);
}

DrawCurves* ImgSequence::getRelativeDrawCurves(const int a) const {
  if (!(a >= 0 && a <= _endF-_startF))
    return NULL;
  return (_drawc+a);
}

KLT_FullPyramid* ImgSequence::getFullPyramid(const int a) {
  if (!(a >= _startF && a <= _endF))
    return NULL;
  KLT_FullPyramid *fp = _pyrms+a-_startF;
  if (!fp->img)
    fp->initMe(getImage(a), &globalTC);
  return (_pyrms+a-_startF);
}


RotoPath* ImgSequence::resolveRotoWritePtr(int* ptr) { 
  //assert(ptr[0] >= _startF && ptr[0] <= _endF);
  if (!(ptr[0] >= _startF && ptr[0] <= _endF)) { // G
    printf("Big problem 1\n");
    return NULL;
  }
  RotoCurves* rc = _rotoc+ ptr[0]-_startF;
  //assert(ptr[1] >=0 ) ; //G!
  //assert( ptr[1] < rc->getNumCurves());
  if (!(ptr[1] >=0 && ptr[1] < rc->getNumCurves())) {
  printf("Big problem 2\n");
  return NULL;
  }
  return rc->getCurve(ptr[1]);
}

DrawPath* ImgSequence::resolveDrawWritePtr(int* ptr) {
  //assert(ptr[0] >= _startF && ptr[0] <= _endF);
  if (!(ptr[0] >= _startF && ptr[0] <= _endF)) { // G
  printf("D Big problem 1\n ");
  return NULL;
  }
  DrawCurves* dc = _drawc+ ptr[0]-_startF;
  //assert(ptr[1] >=0 && ptr[1] < dc->getNumCurves());
  if (!(ptr[1] >=0 && ptr[1] < dc->getNumCurves())) {
  printf("D Big problem 2\n");
  return NULL;
  }
  return dc->getCurve(ptr[1]);
}

/*
DrawPatch* ImgSequence::resolveDrawPatchWritePtr(int* ptr) {
  assert(ptr[0] >= _startF && ptr[0] <= _endF);
  DrawCurves* dc = _drawc+ ptr[0]-_startF;
  assert(ptr[1] >=0 && ptr[1] < dc->getNumPatches());
  return dc->getPatch(ptr[1]);
}
*/

void ImgSequence::saveMattes() {
  int i;
  for (i = _startF; i<=_endF; i++) 
    _rotoc[i-_startF].saveMatte(_w, _h);
}


void ImgSequence::saveRotoXML(std::ofstream& fp) {
  int i;
  for (i = _startF; i<=_endF; i++) 
    _rotoc[i-_startF].clearXMLLabels();
  int labelNum=0;
  for (i = _startF; i<=_endF; i++) 
    _rotoc[i-_startF].createXMLLabels(labelNum);
  
  fp << "<Movie: {name.mov}>" << std::endl;
  for (i = _startF; i<=_endF; i++) 
    _rotoc[i-_startF].saveRotoXML(fp, i);
  //fp << ">" << std::endl;
}

void ImgSequence::saveCurvesAs(const char* root) {


  strcpy(_fileRoot,root);
  printf("file root is now %s\n",_fileRoot);
  
  saveCurvesQt();
}

void ImgSequence::saveCurves() {
  char name[300];
  FILE* fp;
  int i;

  // write param file
  sprintf(name,"%s.param",_fileRoot);
  fp = fopen(name,"w");
  if (!fp) {
    printf("Failed to save.  Perhaps you need to create the directory?\n");
    return;
  }
  fprintf(fp,"version: %d\n",_VERSION_);   // version number 
  //fprintf(fp,"%s\n",_fileRoot);
  fprintf(fp,"%s\n",_imageRoot);
  fprintf(fp,"%d %d\n\n", _startF, _endF);

  for (i= _startF; i<=_endF; i++) {
    fprintf(fp,"%d %d %d\n", _rotoc[i-_startF].getNumCurves(),
	    _drawc[i-_startF].getNumCurves(),0);
  }
  fclose(fp);
  
  // prepare for writing
  for (i = _startF; i<=_endF; i++) {
    _rotoc[i-_startF].documentSelf(); // roto then draw, important
    _drawc[i-_startF].documentSelf();
  }

  // actually write
  for (i = _startF; i<=_endF; i++) {
    if (_rotoc[i-_startF].getNumCurves() > 0 || _drawc[i-_startF].getNumCurves() > 0) {
      sprintf(name,"%s%.3d.rd", _fileRoot,i);
      fp = fopen(name,"w");
      assert(fp);
      _rotoc[i-_startF].save(fp, _VERSION_);
      _drawc[i-_startF].save(fp, _VERSION_);
      fclose(fp);
    }
  }

  printf("Wrote files to %s\n", _fileRoot);
}



void ImgSequence::saveCurvesQt() {
  char name[300];
  int i;

  // write param file
  sprintf(name,"%s.param",_fileRoot);
  QFile qf(name);
  bool res = qf.open(IO_WriteOnly);
  if (!res) {
    printf("Failed to save.  Perhaps you need to create the directory?\n");
    return;
  }
  QTextStream qts(&qf);
  qts << "version: " << _VERSION_ << endl;  // version number 
  qts << _imageRoot << endl;
  qts << _startF << " " << _endF << endl << endl;

  for (i= _startF; i<=_endF; i++) {
    qts << _rotoc[i-_startF].getNumCurves() << " "
	<< _drawc[i-_startF].getNumCurves() << " 0"<< endl;
  }
  qf.close();
  
  // prepare for writing
  for (i = _startF; i<=_endF; i++) {
    _rotoc[i-_startF].documentSelf(); // roto then draw, important
    _drawc[i-_startF].documentSelf();
  }

  // actually write
  for (i = _startF; i<=_endF; i++) {
    if (_rotoc[i-_startF].getNumCurves() > 0 || _drawc[i-_startF].getNumCurves() > 0) {
      sprintf(name,"%s%.3d.rd", _fileRoot,i);
      QFile qf(name); 
      bool res = qf.open(IO_WriteOnly);
      assert(res);
      QDataStream qd(&qf);
      _rotoc[i-_startF].saveqt(&qd, _VERSION_);
      _drawc[i-_startF].saveqt(&qd, _VERSION_);
      qf.close();
    }
  }

  printf("Wrote files to %s\n", _fileRoot);
}


void ImgSequence::loadInitialData(FILE* fp, const char* root) {

  int res;
  char s_version[20];
  // load param file
  assert(fp);
  strcpy(_fileRoot, root);

  res = fscanf(fp,"%s",s_version);
  assert(res==1);  assert(strcmp(s_version,"version:")==0);
  res = fscanf(fp,"%d",&_version); assert(res==1);

  //res = fscanf(fp,"%s", _fileRoot); assert(res==1);
  res = fscanf(fp,"%s", _imageRoot); assert(res==1);
  res = fscanf(fp,"%d %d", &_startF, &_endF); assert(res==2);
  mutualInit();  
}
void ImgSequence::loadRestData(FILE* fp) {


  char name[300];
  int i, res;
  bool files=false;
  
  if (!feof(fp)) {
    files = true;
    int rnum, dnum, rpnum;
    for (i= _startF; i<=_endF; i++) {
      res = fscanf(fp,"%d %d %d",&rnum, &dnum, &rpnum); assert(res==3);
      _rotoc[i-_startF].addNPaths(rnum);
      _drawc[i-_startF].addNPaths(dnum);
      //_drawc[i-_startF].addNPatches(rpnum);
    }
  }
  
  if (files) {
    //printf("files\n");
    for (i=_startF; i<=_endF; i++) {
      if (_rotoc[i-_startF].getNumCurves() > 0 || _drawc[i-_startF].getNumCurves() > 0) {
	sprintf(name,"%s%.3d.rd", _fileRoot,i);
	FILE* outfp = fopen(name,"r");
	assert(outfp);
	_rotoc[i-_startF].load(outfp, _version);
	_drawc[i-_startF].load(outfp, _version);
	fclose(outfp);
      }
    }

  }
}


void ImgSequence::loadInitialDataqt(QFile* fp, const char* root) {


  char s_version[20];
  // load param file
  assert(fp);
  strcpy(_fileRoot, root);
  QTextStream qtm(fp);

  qtm >> s_version;
  assert(strcmp(s_version,"version:")==0);
  qtm >> _version;  assert (_version <= _VERSION_);

  qtm >> _imageRoot;
  qtm >> _startF >> _endF;
  mutualInit();  
}

void ImgSequence::loadRestDataqt(QFile* fp) {


  char name[300];
  int i;
  bool files=false;
  QTextStream qtm(fp);
  
  if (!qtm.atEnd()) {
    files = true;
    int rnum, dnum, rpnum;
    for (i= _startF; i<=_endF; i++) {
      qtm >> rnum >> dnum >> rpnum;
      _rotoc[i-_startF].addNPaths(rnum);
      _drawc[i-_startF].addNPaths(dnum);
      //_drawc[i-_startF].addNPatches(rpnum);
    }
  }
  
  if (files) {
    //printf("files\n");
    for (i=_startF; i<=_endF; i++) {
      if (_rotoc[i-_startF].getNumCurves() > 0 || _drawc[i-_startF].getNumCurves() > 0) {
	sprintf(name,"%s%.3d.rd", _fileRoot,i);
	QFile subqfp(name);
	subqfp.open(IO_ReadOnly);
	QDataStream subqds(&subqfp);
	_rotoc[i-_startF].loadqt(&subqds, _version);
	_drawc[i-_startF].loadqt(&subqds, _version);
	subqfp.close();
      }
    }

  }
}


KLT_FullCPyramid* ImgSequence::getFullCPyramid(const int a) {
  if (!(a >= _startF && a <= _endF))
    return NULL;
  KLT_FullCPyramid *fp = _Cpyrms+a-_startF;
  if (!fp->img)
    if (!loadFullCPyramid(a)) {
      fp->initMe(getImage(a), &globalTC);
      printf("Calculated frame %d\n",a);
    }
  return (_Cpyrms+a-_startF);
}

void ImgSequence::loadCalcFullCPyramids(const int a, const int b) {
  for (int i=a; i<=b; i++) {
    KLT_FullCPyramid *fp = _Cpyrms+i-_startF;
    if (!fp->img)
      if (!loadFullCPyramid(i))
	fp->initMe(getImage(i), &globalTC);
    assert(fp->img);
  }
}


void ImgSequence::writeFullCPyramid(const int a) { 
  if (!(a >= _startF && a <= _endF)) {
    printf("Failed to write (no such pyramid)");
    return;
  }

  KLT_FullCPyramid *fp = _Cpyrms+a-_startF;  
  if (!fp->img)
    fp->initMe(getImage(a), &globalTC);

  sprintf(_fileName,"%s%.3d.fcp",_imageRoot,a);
  FILE* fileptr = fopen(_fileName,"w");
  assert(fileptr);
  fp->write(fileptr);
  fclose(fileptr);
  printf("Wrote file %s\n",_fileName);
}


void ImgSequence::writeFullCPyramidImage(const int a) { 
  if (!(a >= _startF && a <= _endF)) {
    printf("Failed to write (no such pyramid)");
    return;
  }
  KLT_FullCPyramid *fp = _Cpyrms+a-_startF;  
  if (!fp->img)
    fp->initMe(getImage(a), &globalTC);

  char name1[200], name2[200], name3[200];
  sprintf(name1,"%s%.3d",_imageRoot,a);
  sprintf(name2,"%s%.3d_dx",_imageRoot,a);
  sprintf(name3,"%s%.3d_dy",_imageRoot,a);
  fp->writeImages(name1, name2, name3);
}

bool ImgSequence::loadFullCPyramid(const int a) { 
  if (!(a >= _startF && a <= _endF))
    return false;
  KLT_FullCPyramid *fp = _Cpyrms+a-_startF;
  assert(!fp->img); // other procedures should check for this

  sprintf(_fileName,"%s%.3d.fcp",_imageRoot,a);
  FILE* fileptr = fopen(_fileName,"r");
  if (!fileptr) return false;
  printf("Found file for %d\n",a);
  bool res = fp->load(fileptr, &globalTC);
  fclose(fileptr);
  if (res)
    return true;
  else {
    printf("Found file different number of levels\n");
    return false;
  }
}


//-------------
KLT_FullPyramid* ImgSequence::getFullEdgePyramid(const int a) {
  if (!(a >= _startF && a <= _endF))
    return NULL;
  KLT_FullPyramid *fp = _pyrmsE+a-_startF;
  if (!fp->img)
    if (!loadFullEdgePyramid(a)) {
      fp->initMeFromEdges(getImage(a), &globalTC);
      printf("Calculated edge frame %d\n",a);
    }
  return (_pyrmsE + a-_startF);
}

void ImgSequence::loadCalcFullEdgePyramids(const int a, const int b) {
  for (int i=a; i<=b; i++) {
    KLT_FullPyramid *fp = _pyrmsE + i-_startF;
    if (!fp->img)
      if (!loadFullEdgePyramid(i))
	fp->initMeFromEdges(getImage(i), &globalTC);
    assert(fp->img);
  }
}


void ImgSequence::writeFullEdgePyramid(const int a) { 
  if (!(a >= _startF && a <= _endF)) {
    printf("Failed to write (no such pyramid)");
    return;
  }
  KLT_FullPyramid *fp = _pyrmsE + a-_startF;  
  if (!fp->img)
    fp->initMeFromEdges(getImage(a), &globalTC);

  sprintf(_fileName,"%s%.3d.ep",_imageRoot,a);
  FILE* fileptr = fopen(_fileName,"w");
  assert(fileptr);
  fp->write(fileptr);
  fclose(fileptr);
  printf("Wrote file %s\n",_fileName);
}

bool ImgSequence::loadFullEdgePyramid(const int a) { 
  if (!(a >= _startF && a <= _endF))
    return false;
  KLT_FullPyramid *fp = _pyrmsE + a-_startF;
  assert(!fp->img); // other procedures should check for this

  sprintf(_fileName,"%s%.3d.ep",_imageRoot,a);
  FILE* fileptr = fopen(_fileName,"r");
  if (!fileptr) return false;
  printf("Found edge file for %d\n",a);
  bool res = fp->load(fileptr, &globalTC);
  fclose(fileptr);
  if (res)
    return true;
  else {
    printf("Found edge file different number of levels\n");
    return false;
  }
}

void ImgSequence::writeFullEdgePyramidImage(const int a) { 
  if (!(a >= _startF && a <= _endF)) {
    printf("Failed to write (no such pyramid)");
    return;
  }
  KLT_FullPyramid *fp = _pyrmsE+a-_startF;  
  if (!fp->img)
    fp->initMeFromEdges(getImage(a), &globalTC);

  char name1[200], name2[200], name3[200];
  sprintf(name1,"%s%.3d",_imageRoot,a);
  sprintf(name2,"%s%.3d_dx",_imageRoot,a);
  sprintf(name3,"%s%.3d_dy",_imageRoot,a);
  fp->writeImages(name1, name2, name3);
}



//--------------


void ImgSequence::calculateStatistics() {
  printf("Statistics:\n");
  Vec2i d, t;
  for (int i= _startF; i<=_endF; i++) {
    t += _rotoc[i-_startF].calcStat();
    d += _drawc[i-_startF].calcStat();
  }

  printf("Draw: %d touched out of %d\n",d.x(), d.y());
  printf("Roto: %d touched out of %d\n",t.x(), t.y());
}

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



#include <math.h>
#include <float.h>
#include "imageGraph.h"


ImageGraph::ImageGraph(const QImage im) {  
  _w = im.width();
  _h = im.height();
  _maxI = _w*_h;
  _graph = new GraphNode[_maxI];
  float maxD=-1;

  int i,j, k,  pix;
  for (j=0; j<_h; j++) {
    for (i=0; i<_w; i++) {
      pix = j*_w+i;
      _graph[pix]._loc.set_x(i);
      _graph[pix]._loc.set_y(j);
      _graph[pix]._totalCost = FLT_MAX;

      // horizontal and vertical links
      //_graph[pix]._linkCost[0] = horizLink(im,
      //				   i,j-1, i+1,j-1,
      //				   i,j+1, i+1,j+1);
      _graph[pix]._linkCost[4] = horizLink(im,
					   i,j-1, i-1,j-1,
					   i,j+1, i-1,j+1);
      if (i!=0) _graph[pix-1]._linkCost[0] = _graph[pix]._linkCost[4];
      //_graph[pix]._linkCost[2] = horizLink(im,
      //				   i-1,j, i-1,j-1, 
      //				   i+1,j, i+1,j-1);
      _graph[pix]._linkCost[6] = horizLink(im,
					   i-1,j, i-1,j+1,
					   i+1,j, i+1,j+1);
      if (j!=_h-1) 
	_graph[(j+1)*_w+i]._linkCost[2] =  _graph[pix]._linkCost[6];

      // diagonal links
      //_graph[pix]._linkCost[1] = diagLink(im,
      //				  i+1,j, i,j-1);
      //_graph[pix]._linkCost[3] = diagLink(im,
      //				  i-1,j, i,j-1);
      _graph[pix]._linkCost[5] = diagLink(im,
					  i-1,j, i,j+1);
      if (i!=0 && j!=_h-1)
	_graph[(j+1)*_w+i-1]._linkCost[1] =  _graph[pix]._linkCost[5];
      _graph[pix]._linkCost[7] = diagLink(im,
					  i,j+1, i+1,j);
      if (i!=_w-1 && j!=_h-1)
	_graph[(j+1)*_w+i+1]._linkCost[3] =  _graph[pix]._linkCost[7];

      // record maximum derivative strength
      for (k=4; k<8; k++) {
	assert(_graph[pix]._linkCost[k] >= 0);
	if (_graph[pix]._linkCost[k] > maxD)
	  maxD = _graph[pix]._linkCost[k];
      }
    }
  }

  // compute final costs
  for (j=0; j<_h; j++) {
    for (i=0; i<_w; i++) {
      pix = j*_w+i;
      // final divide-by-2 is just to bring down value to under 255
      for (k=0; k<8; k+=2)
	_graph[pix]._linkCost[k] = (maxD-_graph[pix]._linkCost[k]) / 2.;
      for (k=1; k<8; k+=2)
	_graph[pix]._linkCost[k] = (maxD-_graph[pix]._linkCost[k])*M_SQRT2 / 2.;
    }
  }
  
  printf("finished initing graph\n");
}


// SPEED: should store floats, cut file size and file load time in half
ImageGraph::ImageGraph(FILE* fp) {
  if (!fp) {
    printf("Couldn't open image graph\n");
    exit(0);
  }

  int res,i,j,pix;
  res = fread(&_w,sizeof(int),1,fp); assert(res==1);
  res = fread(&_h,sizeof(int),1,fp); assert(res==1);
  _maxI = _w*_h;

  _graph = new GraphNode[_maxI];
  for  (j=0; j<_h; j++) 
    for (i=0; i<_w; i++) {
      pix = j*_w+i;
      _graph[pix]._loc.set_x(i);
      _graph[pix]._loc.set_y(j);
      _graph[pix]._totalCost = FLT_MAX;

      res = fread(_graph[pix]._linkCost+4,sizeof(float),4,fp);
      //if (j==0) printf("%f %f %f %f\n", _graph[pix]._linkCost[4],_graph[pix]._linkCost[5],
      //	       _graph[pix]._linkCost[6],_graph[pix]._linkCost[7]);
      assert(res==4);

      // reflect values
      if (i!=0) _graph[pix-1]._linkCost[0] = _graph[pix]._linkCost[4];
      if (j!=_h-1) 
	_graph[pix+_w]._linkCost[2] =  _graph[pix]._linkCost[6];
      if (i!=0 && j!=_h-1)
	_graph[pix+_w-1]._linkCost[1] =  _graph[pix]._linkCost[5];
      if (i!=_w-1 && j!=_h-1)
	_graph[pix+_w+1]._linkCost[3] =  _graph[pix]._linkCost[7];
    }
  fclose(fp);
}


float ImageGraph::horizLink(const QImage im, const int x1, const int y1, 
		 const int x2, const int y2,
		 const int x3, const int y3, 
		 const int x4, const int y4) {

  Vec3i a = paddedImageQuery(im,x1,y1), 
    b = paddedImageQuery(im,x2,y2),
    c = paddedImageQuery(im,x3,y3),
    d = paddedImageQuery(im,x4,y4);
  Vec3f res;
  res.set_r(fabs(float(a.r() + b.r() - (c.r()+d.r())))/4.0);
  res.set_g(fabs(float(a.g() + b.g() - (c.g()+d.g())))/4.0);
  res.set_b(fabs(float(a.b() + b.b() - (c.b()+d.b())))/4.0);
  
  return sqrt(res.Len2()/3.0); // should be 3, avoiding square root
}

float ImageGraph::diagLink(const QImage im, const int x1, const int y1, 
		const int x2, const int y2) {
  Vec3i a = paddedImageQuery(im,x1,y1), 
    b = paddedImageQuery(im,x2,y2);
  
  Vec3f res;
  res.set_r(fabs(a.r() - b.r())/M_SQRT2);
  res.set_g(fabs(a.g() - b.g())/M_SQRT2);
  res.set_b(fabs(a.b() - b.b())/M_SQRT2);

  return sqrt(res.Len2()/3.0); // should be 3, avoiding square root
}

Vec3i ImageGraph::paddedImageQuery(const QImage im, int i, int j) {
  int ii,jj;
  if (i<0) 
    ii=0;
  else if (i>= im.width())
    ii = im.width()-1;
  else
    ii = i;
  if (j<0) 
    jj=0;
  else if (j>= im.height())
    jj = im.height()-1;
  else
    jj = j;
  QRgb rgb = im.pixel(ii,jj);
  return Vec3i(qRed(rgb), qGreen(rgb), qBlue(rgb));
}


void ImageGraph::init() {
  int i,j;
  for (j=0; j<_h; j++) {
    for (i=0; i<_w; i++) {
      int pix = j*_w+i;
      _graph[pix]._state = INITIAL;
      _graph[pix]._prevNode = NULL;
    }
  }
}

void ImageGraph::save(const char* fileName) const {
  FILE* fp = fopen(fileName,"w");
  assert(fp);
  int i,j,res;
  res = fwrite(&_w,sizeof(int),1,fp); assert(res==1);
  res = fwrite(&_h,sizeof(int),1,fp); assert(res==1);
  //  printf("saving to %s, %d %d\n",fileName, _w, _h);
  for  (j=0; j<_h; j++) 
    for (i=0; i<_w; i++) { 
      res = fwrite(_graph[j*_w+i]._linkCost+4,sizeof(float),4,fp);

      assert(res==4);
    }
  fclose(fp);
}

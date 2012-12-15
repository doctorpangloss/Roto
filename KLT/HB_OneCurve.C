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



#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <qdatetime.h>
#include "HB_OneCurve.h"


#define _SCA(x,y) (2*((y)*_w + (x)))
// takes advantage that only space endpoints are fixed
#define _fixed(x,y) (((x)==0 || (x)==_w-1) && (find(_fixedLocs.begin(),_fixedLocs.end(), 2*((y)*_w + (x))) != _fixedLocs.end()))
#define _PL_ 3

//HB_OneCurve::HB_OneCurve(int xlevels, int ylevels, int span) : 
HB_OneCurve::HB_OneCurve(int w, int h, vector<int>& fixer, int offset) {
  init(w,h,fixer,offset);
}

void HB_OneCurve::init(int w, int h, vector<int>& fixer, int offset) {
  _w = w; _h=h;
  _maxdim = MAX(_w,_h);
  _qu = new ushort[2*_maxdim*_maxdim];  assert(_qu);
  memset(_qu,0,2*_maxdim*_maxdim*sizeof(ushort));
  _totalVar = _w*_h*2;
  _fixedLocs = fixer;
  for (unsigned int i=0; i<_fixedLocs.size(); ++i) {
    _fixedLocs[i] -= offset;
    assert(_fixedLocs[i] < _totalVar-1);
    //printf("%d,%d /2 space fixed\n", (_fixedLocs[i]/2)%_w, (_fixedLocs[i]/2)/_w);
  }
  memset(_dummy,0,2*sizeof(double));
}

// multiply by S'
void HB_OneCurve::sweepUp(double* r) const {
  
  //ushort* debugGrid = new ushort[_w*_h]; // DEBUG
  //memset(debugGrid, 0, _w*_h*sizeof(ushort));
  //assertFixed(r); 
  QTime time;
  time.start();

  int qu_h= 0; // treat as stack, now
  ushort curr[5] = {0,0,_w-1,_h-1, 0}; // last number is # children traversed (0,1,2,3,4). 5 means dead, pop or terminate
  ushort mx,my,dx,dy,i,j;
  double k,m;
  double *ll,*hl,*lh,*hh;

  do {   // drive to leaf


    dx =  curr[2] - curr[0];
    dy =  curr[3] - curr[1];
    mx =  curr[0] + (dx>>1);
    my =  curr[1] + (dy>>1);
    //assert(dx>=0 && dy>=0);
    //l = MAX(dx,dy);    

    ll = _fixed(curr[0], curr[1]) ? _dummy : r+_SCA(curr[0], curr[1]);
    hl = _fixed(curr[2], curr[1]) ? _dummy : r+_SCA(curr[2], curr[1]);
    lh = _fixed(curr[0], curr[3]) ? _dummy : r+_SCA(curr[0], curr[3]);
    hh = _fixed(curr[2], curr[3]) ? _dummy : r+_SCA(curr[2], curr[3]);

    
  
    if (curr[4]==4 && dx>1 && dy>1) {  // done with this 4-interpolant, process & continue

      // do interpolation
      k = double(mx-curr[0]) / double(dx); 
      m = double(my-curr[1]) / double(dy); 
      if (curr[1]==0) { // bottom
	ll[0] += (1.-k)*r[_SCA(mx,0)];      ll[1] += (1.-k)*r[_SCA(mx,0)+1]; 
	hl[0] += k*r[_SCA(mx,0)];           hl[1] += k*r[_SCA(mx,0)+1];   
      }
      hl[0] += (1.-m)*r[_SCA(curr[2],my)];  hl[1] += (1.-m)*r[_SCA(curr[2],my)+1]; // right
      hh[0] += m*r[_SCA(curr[2],my)];       hh[1] += m*r[_SCA(curr[2],my)+1];
      if (curr[0]==0) {
	ll[0] += (1.-m)*r[_SCA(0,my)];      ll[1] += (1.-m)*r[_SCA(0,my)+1];  //left
	lh[0] += m*r[_SCA(0,my)];           lh[1] += m*r[_SCA(0,my)+1];
      }
      lh[0] += (1.-k)*r[_SCA(mx,curr[3])]; lh[1] += (1.-k)*r[_SCA(mx,curr[3])+1]; // top
      hh[0] += k*r[_SCA(mx,curr[3])];      hh[1] += k*r[_SCA(mx,curr[3])+1];
      
      // center
      ll[0] += (1.-k)*(1.-m)*r[_SCA(mx,my)]; ll[1] += (1.-k)*(1.-m)*r[_SCA(mx,my)+1];
      hl[0] += k*(1.-m)*r[_SCA(mx,my)];      hl[1] += k*(1.-m)*r[_SCA(mx,my)+1];
      lh[0] += (1.-k)*m*r[_SCA(mx,my)];      lh[1] += (1.-k)*m*r[_SCA(mx,my)+1];
      hh[0] += k*m*r[_SCA(mx,my)];           hh[1] += k*m*r[_SCA(mx,my)+1];
      
      /*
	debugGrid[my*_w + mx] += l;
	if (curr[1]==0)      
	debugGrid[mx] += l;                //bottom
	debugGrid[my*_w + curr[2]] += l; // right
	if (curr[0]==0)
	debugGrid[my*_w] +=l;              // left
	debugGrid[curr[3]*_w + mx] += l; // top*/
      
      curr[4]=5; // induce pop
    }
    
    else if (curr[4]==2 && (dx==1 || dy==1)) { // done with this 2-interpolant, process & continue
      // do interpolation
      if (dx==1) {
	m = double(my-curr[1]) / double(dy); 
	if (curr[0]==0) {  // left
	  ll[0] += (1.-m)* r[_SCA(0,my)];       ll[1] += (1.-m)* r[_SCA(0,my)+1];
	  lh[0] += m * r[_SCA(0,my)];           lh[1] += m * r[_SCA(0,my)+1];
	}
	hl[0] += (1.-m)* r[_SCA(curr[2],my)];   hl[1] += (1.-m)* r[_SCA(curr[2],my)+1]; // right
	hh[0] += m * r[_SCA(curr[2],my)];       hh[1] += m * r[_SCA(curr[2],my)+1];
      }
      else {
	k = double(mx-curr[0]) / double(dx); 
	if (_qu[qu_h+1]==0) { // bottom
	  ll[0] += (1.-k)* r[_SCA(mx,0)];       ll[1] += (1.-k)* r[_SCA(mx,0)+1];
	  hl[0] +=k* r[_SCA(mx,0)];             hl[1] += k* r[_SCA(mx,0)+1]; 
	}
	lh[0] += (1.-k)*r[_SCA(mx,curr[3])];    lh[1] += (1.-k)*r[_SCA(mx,curr[3])+1];
	hh[0] += k * r[_SCA(mx,curr[3])];       hh[1] += k * r[_SCA(mx,curr[3])+1];	 
      }
      
      curr[4]=5; // induce pop
    }
    
    else if ((dx>_PL_  && dy>1) || (dy>_PL_ && dx>1)) {  // subdivide, if we have children left add to stack 
      
      assert(curr[4] < 4);
      memcpy(_qu+qu_h, curr, 5*sizeof(ushort));
      qu_h += 5;
      _qu[qu_h-1]++; // we're about to traverse one
      
      if (curr[4]==0) {       // hh
	curr[0] = mx;     curr[1] = my;
      }
      else if (curr[4]==1) {  // lh
	curr[1] = my;   curr[2] = mx;    
      }
      else if (curr[4]==2) { // hl
	curr[0] = mx;     curr[3] = my;
      }
      else if (curr[4]==3) { // hh
	curr[2] = mx;   curr[3] = my;
      }
      curr[4]=0;      
    }
    
    else if (dx<=_PL_ && dy<=_PL_) { // leaf node, handle 
      // handle 1x2,1x3,2x2,2x3,3x3 or flipped (3 more for total 8)
      // separate axes

      for (i=curr[0]+1; i<curr[2]; ++i) { // horizontal edges
	k = double(i-curr[0]) / double(dx);  
	if (curr[1]==0) {  // bottom
	  ll[0] += (1.-k)*r[_SCA(i,0)];   ll[1] += (1.-k)*r[_SCA(i,0)+1];
	  hl[0] += k*r[_SCA(i,0)];        hl[1] += k*r[_SCA(i,0)+1];
	}
	lh[0] += (1.-k)*r[_SCA(i,curr[3])];  lh[1] += (1.-k)*r[_SCA(i,curr[3])+1]; // top
	hh[0] += k*r[_SCA(i,curr[3])];       hh[1] += k*r[_SCA(i,curr[3])+1];	
      }
      
      for (j=curr[1]+1; j<curr[3]; ++j) { // vertical edges
	k = double(j-curr[1]) / double(dy); 
	if (curr[0]==0) {        // left
	  ll[0] += (1.-k)*r[_SCA(0,j)];  ll[1] += (1.-k)*r[_SCA(0,j)+1]; 
	  lh[0] += k*r[_SCA(0,j)];       lh[1] += k*r[_SCA(0,j)+1]; 
	}
	hl[0] += (1.-k)*r[_SCA(curr[2], j)];  hl[1] += (1.-k)*r[_SCA(curr[2], j)+1];  // right
	hh[0] += k*r[_SCA(curr[2], j)];       hh[1] += k*r[_SCA(curr[2], j)+1];  
      }
      
      for (j=curr[1]+1; j<curr[3]; ++j) 
	for (i=curr[0]+1; i<curr[2]; ++i) {
	  k = double(i-curr[0]) / double(dx); 
	  m = double(j-curr[1]) / double(dy); 
	  ll[0] += (1.-k)*(1.-m)*r[_SCA(i,j)]; ll[1] += (1.-k)*(1.-m)*r[_SCA(i,j)+1];
	  hl[0] += k*(1.-m)*r[_SCA(i,j)];      hl[1] += k*(1.-m)*r[_SCA(i,j)+1];
	  lh[0] += (1.-k)*m*r[_SCA(i,j)];      lh[1] += (1.-k)*m*r[_SCA(i,j)+1];
	  hh[0] += k*m*r[_SCA(i,j)];           hh[1] += k*m*r[_SCA(i,j)+1];
	}

      /*
      for (i=curr[0]+1; i<curr[2]; ++i) { // horizontal edges
	if (curr[1]==0)   
	  debugGrid[i] += l;  // bottom
	debugGrid[curr[3]*_w + i] += l;  // top
      }
      
      for (j=curr[1]+1; j<curr[3]; ++j) { // vertical edges
	if (curr[0]==0)         // left
	  debugGrid[j*_w] += l;
	debugGrid[j*_w + curr[2]] += l; // right
      }
      
      for (j=curr[1]+1; j<curr[3]; ++j) 
	for (i=curr[0]+1; i<curr[2]; ++i) {
	  debugGrid[j*_w + i] += l;
	}
      */

      curr[4] = 5;
    }

    else if (dx==1) {

      memcpy(_qu+qu_h, curr, 5*sizeof(ushort));
      qu_h += 5;
      _qu[qu_h-1]++; // we're about to traverse one

      if (curr[4]==0)      // h
	curr[1]=my;
      else if (curr[4]==1) // l
	curr[3]=my;
      curr[4]=0;      
    }

    else if (dy==1) {

      memcpy(_qu+qu_h, curr, 5*sizeof(ushort));
      qu_h += 5;
      _qu[qu_h-1]++; // we're about to traverse one

      if (curr[4]==0)      // h
	curr[0]=mx;
      else if (curr[4]==1) // l
	curr[2]=mx;
      curr[4]=0;
    }


    // pop
    if (curr[4]==5) {  // leaf handled, need to pop stack
      if (qu_h==0) break; // stack empty, done, terminate
      qu_h-=5;
      memcpy(curr,_qu+qu_h, 5*sizeof(ushort));
    }
    
    /*int index=0,ii,jj;
    for (jj=0; jj<_h; ++jj) {
      for (ii=0; ii<_w; ++ii, index+=2)
	printf("%.2f,%.2f   ", r[index],r[index+1]);
      printf("\n");
    }
    printf("\n");*/

    
  } while (curr[4]!=5);  // end drive to leaf

  // assertFixed(r); 

  /*
  int index=0;
  for (j=0; j<_h; ++j) {
    for (i=0; i<_w; ++i, ++index)
      printf("%.2d ", debugGrid[index]);
    printf("\n");
  }
  printf("\n");
  */

  //printf("%d time taken to sweep\n",time.elapsed());

  //delete[] debugGrid; // DEBUG

} 


// multipy by S
void HB_OneCurve::sweepDown(double* r) const {
  
  //ushort* debugGrid = new ushort[_w*_h]; // DEBUG
  //memset(debugGrid, 0, _w*_h*sizeof(ushort));
  QTime time;
  time.start();

  // representation if lower-left corner, upper right corner
  int qu_h=2*_maxdim*_maxdim - 4;
  int qu_l = qu_h - 4;
  assert(qu_l>=0);
  _qu[qu_h] = _qu[qu_h+1] = 0;
  _qu[qu_h+2] = _w-1;   _qu[qu_h+3] = _h-1; 
  ushort mx,my,dx,dy,i,j,l;
  double llx,lly,lhx,lhy,hlx,hly,hhx,hhy, k,y1x,y1y,y2x,y2y,m;
  int tmp;

  while (qu_h>qu_l) {
    dx =  _qu[qu_h+2] - _qu[qu_h]; 
    dy =  _qu[qu_h+3] - _qu[qu_h+1]; 
    mx =  _qu[qu_h] + (dx>>1);
    my =  _qu[qu_h+1] + (dy>>1);
    //assert(dx>0 && dy>0);
    l = MAX(dx,dy);
    
    if (!_fixed(_qu[qu_h],_qu[qu_h+1])) {
      llx = r[tmp=_SCA(_qu[qu_h],_qu[qu_h+1])];  lly = r[tmp+1];  }
    else {
      llx = lly = 0;  }
    if (!_fixed(_qu[qu_h+2],_qu[qu_h+1])) {
      hlx = r[tmp=_SCA(_qu[qu_h+2],_qu[qu_h+1])];  hly = r[tmp+1]; }
    else {
      hlx = hly = 0;  }
    if (!_fixed(_qu[qu_h],_qu[qu_h+3])) {
      lhx = r[tmp=_SCA(_qu[qu_h],_qu[qu_h+3])];  lhy = r[tmp+1]; }
    else {
      lhx = lhy = 0;  }
    if (!_fixed(_qu[qu_h+2],_qu[qu_h+3])) {
      hhx = r[tmp=_SCA(_qu[qu_h+2],_qu[qu_h+3])];  hhy = r[tmp+1]; }
    else {
      hhx = hhy = 0;  }
    

    // SPEED: create offset pointers, remove 1 add per access

    if ((dx>_PL_  && dy>1) || (dy>_PL_ && dx>1)) {  //enque subdivides 

      // ll
      _qu[qu_l] = _qu[qu_h];        _qu[qu_l+1] = _qu[qu_h+1];
      _qu[qu_l+2] = mx;             _qu[qu_l+3] = my;
      qu_l -=4;
      //printf("added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      // hl
      _qu[qu_l] = mx;               _qu[qu_l+1] =  _qu[qu_h+1];
      _qu[qu_l+2] =  _qu[qu_h+2];   _qu[qu_l+3] = my;
      qu_l -=4;
      //printf("added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      // hl
      _qu[qu_l] = _qu[qu_h];        _qu[qu_l+1] = my;
      _qu[qu_l+2] = mx;             _qu[qu_l+3] = _qu[qu_h+3];
      qu_l -=4;
      //printf("added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      // hh
      _qu[qu_l] = mx;               _qu[qu_l+1] = my;
      _qu[qu_l+2] = _qu[qu_h+2];    _qu[qu_l+3] = _qu[qu_h+3];
      qu_l -=4;      
      //printf("added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      k = double(mx-_qu[qu_h]) / double(dx); 
      m = double(my-_qu[qu_h+1]) / double(dy); 
      // do actual interpolation
      if (_qu[qu_h+1]==0) { // bottom
	r[_SCA(mx,0)]   += llx + k * (hlx - llx);
	r[_SCA(mx,0)+1] += lly + k * (hly - lly);
      }
      r[_SCA(_qu[qu_h+2],my)]   += hlx + m * (hhx - hlx);   // right				    
      r[_SCA(_qu[qu_h+2],my)+1] += hly + m * (hhy - hly);

      if (_qu[qu_h]==0) {    // left
	r[_SCA(0,my)]   += llx + m * (lhx - llx);
	r[_SCA(0,my)+1] += lly + m * (lhy - lly);
      }
      r[_SCA(mx,_qu[qu_h+3])]   += lhx + k * (hhx - lhx); //top
      r[_SCA(mx,_qu[qu_h+3])+1] += lhy + k * (hhy - lhy);

      y1x = llx + k*(hlx-llx);   y2x = lhx + k*(hhx-lhx); // center
      y1y = lly + k*(hly-lly);   y2y = lhy + k*(hhy-lhy);
      r[_SCA(mx,my)] += y1x + m * (y2x-y1x);
      r[_SCA(mx,my)+1] += y1y + m * (y2y-y1y);

      /*
      debugGrid[my*_w + mx] += l;
      if (_qu[qu_h+1]==0)      
	debugGrid[mx] += l;                //bottom
      debugGrid[my*_w + _qu[qu_h+2]] += l; // right
      if (_qu[qu_h]==0)
	debugGrid[my*_w] +=l;              // left
      debugGrid[_qu[qu_h+3]*_w + mx] += l; // top
      */

    }
    else if (dx<=_PL_ && dy<=_PL_) { // handle 

 
      // handle 1x2,1x3,2x2,2x3,3x3 or flipped (3 more for total 8)
      // separate axes
      // do actual interpolation
      for (i=_qu[qu_h]+1; i<_qu[qu_h+2]; ++i) { // horizontal edges
	k = double(i-_qu[qu_h]) / double(dx);
	if (_qu[qu_h+1]==0) {  // bottom
	  r[_SCA(i,0)] += llx + k * (hlx-llx);
	  r[_SCA(i,0)+1] += lly + k * (hly-lly);
	}
	r[_SCA(i,_qu[qu_h+3])] += lhx + k * (hhx-lhx);   // top
	r[_SCA(i,_qu[qu_h+3])+1] += lhy + k * (hhy-lhy);
      }

      for (j=_qu[qu_h+1]+1; j<_qu[qu_h+3]; ++j) { // vertical edges
	k = double(j-_qu[qu_h+1]) / double(dy);
	if (_qu[qu_h]==0) {         // left
	  r[_SCA(0,j)] += llx + k * (lhx-llx);
	  r[_SCA(0,j)+1] += lly + k * (lhy-lly);
	}
	r[_SCA(_qu[qu_h+2], j)] += hlx + k * (hhx-hlx); // right
	r[_SCA(_qu[qu_h+2], j)+1] += hly + k * (hhy-hly);
      }
      
      for (j=_qu[qu_h+1]+1; j<_qu[qu_h+3]; ++j) 
	for (i=_qu[qu_h]+1; i<_qu[qu_h+2]; ++i) {
	  k = double(i-_qu[qu_h]) / double(dx);
	  m = double(j-_qu[qu_h+1]) / double(dy);
	  y1x = llx + k*(hlx-llx);   y2x = lhx + k*(hhx-lhx);
	  y1y = lly + k*(hly-lly);   y2y = lhy + k*(hhy-lhy);
	  r[_SCA(i,j)] += y1x + m * (y2x-y1x);
	  r[_SCA(i,j)+1] += y1y + m * (y2y-y1y);
	}
      
      /*
      for (i=_qu[qu_h]+1; i<_qu[qu_h+2]; ++i) { // horizontal edges
	if (_qu[qu_h+1]==0)   
	  debugGrid[i] += l;  // bottom
	debugGrid[_qu[qu_h+3]*_w + i] += l;  // top
      }

      for (j=_qu[qu_h+1]+1; j<_qu[qu_h+3]; ++j) { // vertical edges
	if (_qu[qu_h]==0)         // left
	  debugGrid[j*_w] += l;
	debugGrid[j*_w + _qu[qu_h+2]] += l; // right
      }
      
      for (j=_qu[qu_h+1]+1; j<_qu[qu_h+3]; ++j) 
	for (i=_qu[qu_h]+1; i<_qu[qu_h+2]; ++i) {
	  debugGrid[j*_w + i] += l;
	}
      */
    }
    else if (dx==1) {
      //assert(dy>3); 
      // l
      _qu[qu_l] = _qu[qu_h];        _qu[qu_l+1] = _qu[qu_h+1];
      _qu[qu_l+2] = _qu[qu_h+2];    _qu[qu_l+3] = my;
      qu_l -=4;
      //printf("1 added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      // h
      _qu[qu_l] =  _qu[qu_h];       _qu[qu_l+1] = my;
      _qu[qu_l+2] = _qu[qu_h+2];    _qu[qu_l+3] = _qu[qu_h+3];
      qu_l -=4;
      //printf("1 added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      // do actual interpolation
      m = double(my-_qu[qu_h+1]) / double(dy); 
      if (_qu[qu_h]==0) {  // left
	r[_SCA(0,my)]   += llx + m * (lhx - llx);
	r[_SCA(0,my)+1] += lly + m * (lhy - lly);
      }
      r[_SCA(_qu[qu_h+2],my)]   += hlx + m * (hhx - hlx);   // right				    
      r[_SCA(_qu[qu_h+2],my)+1] += hly + m * (hhy - hly);

      //if (_qu[qu_h]==0)
      //debugGrid[my*_w] += l; // left
      //debugGrid[my*_w + _qu[qu_h+2]] += l; // right
      
    }
    else if (dy==1) {
      //assert(dx>3);
      // l
      _qu[qu_l] = _qu[qu_h];        _qu[qu_l+1] = _qu[qu_h+1];
      _qu[qu_l+2] = mx;             _qu[qu_l+3] = _qu[qu_h+3];
      qu_l -=4;
      //printf("1 added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      // h
      _qu[qu_l] = mx;               _qu[qu_l+1] = _qu[qu_h+1];
      _qu[qu_l+2] = _qu[qu_h+2];    _qu[qu_l+3] = _qu[qu_h+3];
      qu_l -=4;
      //printf("1 added %d %d, %d %d\n", _qu[qu_l+4],_qu[qu_l+5],_qu[qu_l+6], _qu[qu_l+7]);

      // do actual interpolation
      k = double(mx-_qu[qu_h]) / double(dx); 
      if (_qu[qu_h+1]==0) { // bottom
	r[_SCA(mx,0)]   += llx + k * (hlx - llx);
	r[_SCA(mx,0)+1] += lly + k * (hly - lly);
      }
      r[_SCA(mx,_qu[qu_h+3])]   += lhx + k * (hhx - lhx); //top
      r[_SCA(mx,_qu[qu_h+3])+1] += lhy + k * (hhy - lhy);
      
      //if (_qu[qu_h+1]==0)
      //debugGrid[mx] += l;
      //debugGrid[_qu[qu_h+3]*_w + mx] += l; // top

    }

    qu_h -=4;
    assert(qu_l >=0);
    
    /*int index=0,ii,jj;
    for (jj=0; jj<_h; ++jj) {
      for (ii=0; ii<_w; ++ii, index+=2)
	printf("%.2f,%.2f   ", r[index],r[index+1]);
      printf("\n");
    }
    printf("\n");*/

  }

  zeroFixed(r);

  //printf("%d time taken to sweep\n",time.elapsed());

  /*int index=0;
  for (j=0; j<_h; ++j) {
    for (i=0; i<_w; ++i, ++index)
      printf("%.2d ", debugGrid[index]);
    printf("\n");
  }
  printf("\n");
  
  delete[] debugGrid; // DEBUG
  */
  //printf("Final q_l %d\n",qu_l);
}


void HB_OneCurve::assertFixed(const double* r) const { 
  int j;
  for (unsigned int i=0; i<_fixedLocs.size(); ++i) {
    j = _fixedLocs[i];
    assert(r[j] == 0 && r[j+1] == 0);
  }
}

void HB_OneCurve::zeroFixed(double* r) const { 
  for (unsigned int i=0; i<_fixedLocs.size(); ++i) 
    r[_fixedLocs[i]] = r[_fixedLocs[i]+1] = 0;
}


double* HB_OneCurve::extractMatrix(int which) const {
  // 0 for S, 1 for S'
  double* out = new double[_totalVar*_totalVar];
  double* bV = new double[_totalVar]; // basis vector
  memset(out,0,_totalVar*_totalVar*sizeof(double));
  int i,j,c;
  assert(which==0 || which==1);

  for (i=0; i<_totalVar; ++i) {
    memset(bV,0,_totalVar*sizeof(double));
    bV[i]=1.;
    zeroFixed(bV);
    if (which==0) 
      sweepDown(bV);
    else
      sweepUp(bV);
    
    for (j=0,c=0; j<_totalVar; ++j, c+=_totalVar)
      out[c+i] = bV[j];
  }


  delete[] bV;
  return out;
}


void HB_OneCurve::printMatrix(double* m, FILE* fp) const {
  int i,j,index=0;
  for (j=0; j<_totalVar; ++j) {
    for (i=0;i<_totalVar; ++i,++index) {
      fprintf(fp,"%.2f ",m[index]);
    }
    fprintf(fp,"\n");
  }
}


/*
// multipy by S
void HB_EvenBlock::sweepDown(double* r) const {
  
  int* debugGrid = new int[_w*_h]; // DEBUG
  memset(debugGrid, 0, _w*_h*sizeof(int));

  int l, sx, sy;
  int xl=_xlevels-1, yl=_ylevels-1; // When these hit 0, don't do
  int xstride = pow2(_xlevels-2), ystride=pow2(_ylevels-2);
  int xstart,ystart,xstep,ystep;

  for (l=_maxLevels-1; l>0; --l, --xl, --yl, xstride>>=1, ystride>>=1) {
    // one pass collapses level l to level l-1
    
    if (xl <=0) {
      xstart=0; xstep=1;
    }
    else {
      xstart=pow2(xl-1); xstep = xstart << 1;
    }
    if (yl <=0) {
      ystart=0; ystep=1;
    }
    else {
      ystart=pow2(yl-1); ystep = ystart << 1;
    }

    for (sy=ystart; sy<_h; sy += ystep) {
      for (sx=xstart; sx<_w; sx += xstep) {
	 
	 
	
	if (xstep==1) {  // interpolate y's only
	  // sx, sy-ystride, up to sx,sy+ystride
	  assert(ystep!=1);
	  debugGrid[sy*_w + sx] = l; // DEBUG
	  
	}
	else if (ystep==1) { // interpolate x's only
	  // sx-xstride, sy, up to sx+xstride,sy
	  assert(xstep!=1);
	  debugGrid[sy*_w + sx] = l; // DEBUG
	} 
	else {  // full subdivide
	  // bbox sx-xstride,sy-ystride, up to sx+xstride, sy+ystride
	  debugGrid[sy*_w + sx+xstride] = l; // DEBUG
	  debugGrid[sy*_w + sx-xstride] = l; // DEBUG
	  debugGrid[(sy-ystride)*_w + sx] = l; // DEBUG
	  debugGrid[(sy+ystride)*_w + sx] = l; // DEBUG
	  debugGrid[sy*_w + sx] = l;
	}

	 // remember to go to full-2 space, etc...
       }  // sx
     } // sy

    
    int index=0;
    for (int j=0; j<_h; ++j) {
      for (int i=0; i<_w; ++i, ++index)
	printf("%.2d ", debugGrid[index]);
      printf("\n");
    }
    printf("\n");
    
    
  } // end over levels

  int index=0;
  for (int j=0; j<_h; ++j) {
    for (int i=0; i<_w; ++i, ++index)
      printf("%.2d ", debugGrid[index]);
    printf("\n");
  }
  printf("\n");

  delete[] debugGrid; // DEBUG
} 
*/

/*
// multiply by S'
void HB_OneCurve::sweepUp(double* r) const {

} 

// multipy by S
void HB_OneCurve::sweepDown(double* r) const {

}
*/



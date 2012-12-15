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



#include "keeper.h"


double Keeper::calculateRo(const double* p, const double newTheta, const double oldTheta) const {
  double num = oldTheta - newTheta;
  double denom = vecDot(numVar(), g(), p);
  double *Bp = new double[numVar()];
  matVecMult(p, Bp);
  denom += .5 * vecDot(numVar(), p, Bp);
  denom *= -1.;
  assert(denom > 0);
  
  double result = num/denom;
  delete[] Bp;
  return result; 
}


void Keeper::refresh() {
  memset(_main,0,_numVar * _WIDTH * sizeof(double));
  memset(_left,0,_offNumVar * _WIDTH * sizeof(double));
  memset(_right,0,_offNumVar * _WIDTH * sizeof(double));
  memset(_g,0,_numVar * sizeof(double));
}


Keeper::Keeper(int num_blocks, int block_size, int curveSampling) {
  init(num_blocks, block_size, curveSampling, 1);
}

void Keeper::init(Keeper* other) {
  _cs = other->_cs;
  _num_blocks = other->_num_blocks;

  _block_size = other->_block_size;
  _o_bs = other->_o_bs;
  _subs = new SubRcd[_o_bs/2];
  memcpy(_subs, other->_subs, (_o_bs/2)*sizeof(SubRcd));

  _numVar = other->_numVar;
  _offNumVar = other->_offNumVar;
   _main = new double[_numVar * _WIDTH];
  memset(_main,0,_numVar * _WIDTH * sizeof(double));
  _left = new double[_offNumVar * _WIDTH];
  memset(_left,0,_offNumVar * _WIDTH * sizeof(double));
  _right = new double[_offNumVar * _WIDTH];
  memset(_right,0,_offNumVar * _WIDTH * sizeof(double));
 
  _g = new double[_numVar];
  memset(_g,0,_numVar * sizeof(double));
}



void Keeper::init(int num_blocks, int block_size, int curveSampling, int hold) {

  _cs = curveSampling;
  _num_blocks = num_blocks;

  // use curveSampling to decide _block_size, _o_bs
  _o_bs = block_size;
  block_size /= 2;                         // in half-coordinates               
  _subs = new SubRcd[block_size];
   hold = MAX(MIN(hold-1, block_size/2-1), 0);
  block_size = MAX(0, block_size - 2*hold); 
  int npb = (block_size-1) / MAX(1,(block_size-1) / _cs), i;
  _block_size = (block_size-1) / npb + 1;

  for (i=0; i<hold; i++)
    _subs[i].set(1.,i);

  int left = (block_size-1) % npb;         // these pts must be distributed into last buckets
  _subs[hold].set(1.,hold);
  int bucket=0, ct=1;                      // which bucket we're dumping into, # points in current bucket
  double k=0;
  for (i=1; i<block_size; i++) {
    if (bucket >= (_block_size-left-1)) {  // into last 'left' buckets
      if (ct >= npb+1) {
	++bucket;
	ct=0;
	k=1.;
      }
      else
	k = 1. - double(ct) / double(npb+1);
    }
    else {
      if (ct >= npb) {
	++bucket;
	ct=0;
	k=1.;
      }
      else
	k = 1. - double(ct) / double(npb);
    }
    
    _subs[i+hold].set(k, bucket+hold);
    assert (bucket < _block_size);
	
    ++ct;
  }
  
   _block_size += 2*hold;
   for (i=0; i<hold; i++)
     _subs[i+block_size+hold].set(1., i+bucket+hold+1);

  assert(k==1);
  assert(bucket = _block_size-1);

  _block_size *= 2;    // back to double coordinates
  printf("Keeper: %d blocks, %d block size\n",_num_blocks, _block_size);

  _numVar = _num_blocks * _block_size;
  _offNumVar = _numVar - _block_size; 
  _main = new double[_numVar * _WIDTH];
  memset(_main,0,_numVar * _WIDTH * sizeof(double));
  _left = new double[_offNumVar * _WIDTH];
  memset(_left,0,_offNumVar * _WIDTH * sizeof(double));
  _right = new double[_offNumVar * _WIDTH];
  memset(_right,0,_offNumVar * _WIDTH * sizeof(double));
  
  _g = new double[_numVar];
  memset(_g,0,_numVar * sizeof(double));
}


void Keeper::add(Keeper* other, double c) {
  assert(other);
  assert(other->_numVar == _numVar);
  assert(other->_offNumVar == _offNumVar);
  int i;
  int ml = _numVar*_WIDTH, ol = _offNumVar*_WIDTH;
  double *optr1, *ptr1;
  double *optr2, *ptr2;

  optr1 = other->_main; ptr1 = _main;
  for (i=0; i<ml; ++i, ++optr1, ++ptr1)
    *ptr1 += c * *optr1;

  optr1 = other->_left;  ptr1 = _left;
  optr2 = other->_right; ptr2 =_right;
  for (i=0; i<ol; ++i, ++optr1, ++ptr1,++optr2, ++ptr2) {
    *ptr1 += c * *optr1;
    *ptr2 += c * *optr2;
  }

  for (i=0; i<_numVar; ++i)
    _g[i] += c * other->_g[i];

  

}

// SPEED: improve this code

void Keeper::matVecMult(const double x[], double res[]) const {  
  memset(res,0,_numVar*sizeof(double));
  double *resptr=res, *ptr=_main;
  int r,c, i, br, bw=0;

  // main strip
  for (r=0,br=0; r<_numVar; ++r, ++resptr,++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	*resptr += *ptr * x[i+bw];
    }
  }
  
  // right strip
  ptr = _right;
  resptr = res;   bw=0;
  for (r=0,br=0; r<_offNumVar; ++r, ++resptr,++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	*resptr += *ptr * x[i+bw+_block_size];
    }
  }

  // left strip
  ptr = _left;   bw=0;
  resptr = res + _block_size;
  for (r=0,br=0; r<_offNumVar; ++r, ++resptr,++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	*resptr += *ptr * x[i+bw];
    }
  }
  
  
  //printf("result:\n");
  /*
  FILE *fp = fopen("/homes/gws/aseem/r.txt","w");
  for (r=0; r<_numVar; r++) 
    fprintf(fp,"%.5f\n",res[r]);
  //printf("\n");
  fclose(fp);

  fp = fopen("/homes/gws/aseem/x.txt","w");
  for (r=0; r<_numVar; r++) 
    fprintf(fp,"%.5f\n",x[r]);
  //printf("\n");
  fclose(fp);

  fp = fopen("/homes/gws/aseem/Z.txt","w");
  printFull(fp);
  fclose(fp);
  
  exit(0);
  */
}


void Keeper::printSparse(int off, FILE* fp) {
  ++off;
  int resptr=0; 
  double *ptr=_main;
  int r,c,i,br,bw=0;
  
  // main strip
  for (r=0,br=0; r<_numVar; ++r, ++resptr,++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	fprintf(fp,"%d %d %.10f\n",resptr+off, i+bw+off, *ptr);
    }
  }
  
  // right strip
  ptr = _right;
  resptr = 0;   bw=0;
  for (r=0,br=0; r<_offNumVar; ++r, ++resptr,++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	fprintf(fp,"%d %d %.10f\n",resptr+off, i+bw+off+_block_size, *ptr);
    }
  }
  
  // left strip
  ptr = _left;   bw=0;
  resptr = _block_size;
  for (r=0,br=0; r<_offNumVar; ++r, ++resptr,++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	fprintf(fp,"%d %d %.10f\n",resptr+off, i+bw+off, *ptr);
    }
  }
  
}


double Keeper::createTestSol(const Vec2f* Z, const double* sol, Vec2f* Z2) const {
  int obs2 = _o_bs/2;
  memcpy(Z2, Z, _num_blocks*obs2*sizeof(Vec2f));

  double maxStep = 0, tmp;
  int Z_off = 0, sol_off=0, v, sol_i, b;
  for (b=0; b<_num_blocks; ++b, sol_off+=_block_size) {  // iterate over blocks
    for (v=0; v < obs2; ++v, ++Z_off) {                    // iterate over variables in a block
      sol_i = sol_off + 2*_subs[v]._P;
      Z2[Z_off].Inc(_subs[v]._k  *  sol[sol_i],
		    _subs[v]._k  *  sol[sol_i+1]);
      if ( (tmp=MAX(fabs(sol[sol_i]), fabs(sol[sol_i+1]) )) > maxStep) {
	maxStep = tmp;
	assert(maxStep < 1000);
      }
      if (_subs[v]._k != 1.) {
	Z2[Z_off].Inc((1.-_subs[v]._k)  *  sol[sol_i+2],
		      (1.-_subs[v]._k)  *  sol[sol_i+3]);
	if ( (tmp=MAX(fabs(sol[sol_i+2]), fabs(sol[sol_i+3]) )) > maxStep) {
	  maxStep = tmp;
	  assert(maxStep < 1000);
	}
      }
    }
  }
  
  //printf("max Step %.5f\n",maxStep);
  /*
  int on = _num_blocks*obs2;
  FILE* fp;
  fp = fopen("Z.txt","w");
  for (b=0; b<on; b++)
    fprintf(fp,"%.5f\n%.5f\n",Z[b].x(), Z[b].y());
  fclose(fp);

  fp = fopen("Z2.txt","w");
  for (b=0; b<on; b++)
    fprintf(fp,"%.5f\n%.5f\n",Z2[b].x(), Z2[b].y());
  fclose(fp);
  exit(0);*/

  return maxStep;
}

void Keeper::printSubs2(FILE* fp, int w, int h) {
  int r;
  int onv = _o_bs / 2;
  for (r=0; r<onv; r++) {
    fprintf(fp, "%d %d %.10f\n",w+2*r+1, h + 2*_subs[r]._P+1, _subs[r]._k);
    fprintf(fp, "%d %d %.10f\n",w+2*r+2, h + 2*_subs[r]._P+2, _subs[r]._k);
    if (_subs[r]._k != 1.) {
      fprintf(fp, "%d %d %.10f\n",w+2*r+1, h + 2*_subs[r]._P + 3, 1. -  _subs[r]._k);
      fprintf(fp, "%d %d %.10f\n",w+2*r + 2, h + 2*_subs[r]._P + 4, 1. -  _subs[r]._k);
    }
  }
}

void Keeper::printSubs(FILE* fp) {
  int r,c;
  int onv = _o_bs / 2;
  for (r=0; r<onv; r++) {
    for (c=0; c<_subs[r]._P; c++)
      fprintf(fp,"0 0 ");
    fprintf(fp,"%.5f 0 ", _subs[r]._k);
    ++c;
    if (_subs[r]._k != 1.) {
       fprintf(fp,"%.5f 0 ", 1. - _subs[r]._k);
       ++c;
    }    
    for (; c<_block_size/2; c++)
      fprintf(fp,"0 0 ");      
    fprintf(fp,"\n");

    for (c=0; c<_subs[r]._P; c++)
      fprintf(fp,"0 0 ");
    fprintf(fp,"0 %.5f ", _subs[r]._k);
    ++c;
    if (_subs[r]._k != 1.) {
       fprintf(fp,"0 %.5f ", 1. - _subs[r]._k);
       ++c;
    }    
    for (; c<_block_size/2; c++)
      fprintf(fp,"0 0 ");      
    fprintf(fp,"\n");
  }
  
}

/*
void Keeper::fill_e(double* e) { 
  memcpy(e,_g, _numVar * sizeof(double));
}
*/

// SPEED: get rid of g1,g2, use _G1, _G2

void Keeper::take2Key0(int ib, const double* g2, double res, const double c) {
  memset(_G1,0,4*sizeof(double));
  int nib, maxj;
  maxj = convert1(ib,1, nib, g2);
  
  if (maxj==4)
    _take4Key0(nib,_G1,res,c);
  else if (maxj==2)
    _take2Key0(nib,_G1,res,c);
  else (assert(0));
}
void Keeper::take2(int b, int ib, const double* g1, const double* g2, double res, const double c) {
  memset(_G1,0,4*sizeof(double));
  memset(_G2,0,4*sizeof(double));

  int nib, maxj;
  maxj = convert2(ib,1, nib, g1, g2);
  
  if (maxj==4)
    _take4(b,nib,_G1, _G2, res, c);
  else if (maxj==2)
    _take2(b,nib,_G1, _G2, res, c);
  else (assert(0));
}
void Keeper::take2Key1(int ib, const double* g1, double res, const double c) {
  memset(_G1,0,4*sizeof(double));
  int nib, maxj;
  maxj = convert1(ib,1, nib, g1);
  
  if (maxj==4)
    _take4Key1(nib,_G1,res,c);
  else if (maxj==2)
    _take2Key1(nib,_G1,res,c);
  else (assert(0));
}
  
void Keeper::take4Key0(int ib, const double* g2, double res, const double c) {
  memset(_G1,0,4*sizeof(double));
  int nib, maxj;
  maxj = convert1(ib,2, nib, g2);
  
  if (maxj==4)
    _take4Key0(nib,_G1,res,c);
  else if (maxj==2)
    _take2Key0(nib,_G1,res,c);
  else (assert(0));
}
/*
void Keeper::take4Keyb(int b, int ib, const double* g2, const double ch, const double cg) {
  memset(_G1,0,4*sizeof(double));
  int nib, maxj;
  maxj = convert1(ib,2, nib, g2);
  
  if (maxj==4)
    _take4Keyb(b, nib,_G1,ch,cg);
  else if (maxj==2)
    _take2Keyb(b, nib,_G1,ch,cg);
  else (assert(0));
  }*/


void Keeper::take4(int b, int ib, const double* g1, const double* g2, double res, const double c) {
  memset(_G1,0,4*sizeof(double));
  memset(_G2,0,4*sizeof(double));

  int nib, maxj;
  maxj = convert2(ib,2, nib, g1, g2);
  
  if (maxj==4)
    _take4(b,nib,_G1, _G2, res, c);
  else if (maxj==2)
    _take2(b,nib,_G1, _G2, res, c);
  else (assert(0));
}


void Keeper::take4Key1(int ib, const double* g1, double res, const double c) {
  memset(_G1,0,4*sizeof(double));
  int nib, maxj;
  maxj = convert1(ib,2, nib, g1);
  
  if (maxj==4)
    _take4Key1(nib,_G1,res,c);
  else if (maxj==2)
    _take2Key1(nib,_G1,res,c);
  else (assert(0));
}

void Keeper::take6Key0(int ib, const double* g2, const double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  int nib, maxj;
  maxj = convert1(ib,3, nib, g2);
  
  if (maxj==6)
    _take6Key0(nib,_G1,res,c);
  else if (maxj==4)
    _take4Key0(nib,_G1,res,c);
  else if (maxj==2)
     _take2Key0(nib,_G1,res,c);
  else (assert(0));
}



void Keeper::take6(int b, int ib, const double* g1, const double* g2, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  memset(_G2,0,6*sizeof(double));

  int nib, maxj;
  maxj = convert2(ib,3, nib, g1, g2);
  
  if (maxj==6)
    _take6(b,nib,_G1, _G2, res, c);
  else if (maxj==4)
    _take4(b,nib,_G1, _G2, res, c);
  else if (maxj==2)
    _take2(b,nib,_G1, _G2, res, c);
  else (assert(0));
}


void Keeper::take6Key1(const int ib, const double* g1, double res, const double c) {
  memset(_G1,0,6*sizeof(double));
  int nib, maxj;
  maxj = convert1(ib,3, nib, g1);
  
  if (maxj==6)
    _take6Key1(nib,_G1,res,c);
  else if (maxj==4)
    _take4Key1(nib,_G1,res,c);
  else if (maxj==2)
    _take2Key1(nib,_G1,res,c);
  else (assert(0));
}


// ib in full 2 space, hm a count in /2 (this is a legacy, should /2 space ib)
// returns full 2 space
int Keeper::convert2(int ib, const int hm, int& nib, const double* g1, const double* g2) {
  return convert2(ib,hm,nib,g1,g2,_G1, _G2);
}

int Keeper::convert2(int ib, const int hm, int& nib, const double* g1, const double* g2, double* G1, double* G2) {
  ib /= 2;
  
  nib = 2 * _subs[ib]._P;
  int i,j, maxj=0;
  for (i=0; i<hm; i++) {
    j = 2 * _subs[ib + i]._P;
    maxj = MAX(maxj,j-nib+2);
    assert (j-nib+1 < 6);
    G1[j-nib] +=   _subs[ib + i]._k * g1[2*i];
    G1[j-nib+1] += _subs[ib + i]._k * g1[2*i+1];
    G2[j-nib] +=   _subs[ib + i]._k * g2[2*i];
    G2[j-nib+1] += _subs[ib + i]._k * g2[2*i+1];
    if (_subs[ib + i]._k != 1.) {
      assert (j-nib+3 < 6);
      maxj = MAX(maxj,j-nib+4);
      G1[j-nib+2] += (1.-_subs[ib + i]._k) * g1[2*i];
      G1[j-nib+3] += (1.-_subs[ib + i]._k) * g1[2*i+1];
      G2[j-nib+2] += (1.-_subs[ib + i]._k) * g2[2*i];
      G2[j-nib+3] += (1.-_subs[ib + i]._k) * g2[2*i+1];
    }
  }
  
  return maxj;
}

int Keeper::convert1(int ib, const int hm, int& nib, const double* g1) {
  return convert1(ib,hm,nib,g1, _G1);
}

int Keeper::convert1(int ib, const int hm, int& nib, const double* g1, double* G1) {
  ib /= 2;
  

  nib = 2 * _subs[ib]._P;
  int i,j, maxj=0;
  for (i=0; i<hm; i++) {
    j = 2 * _subs[ib + i]._P;
    maxj = MAX(maxj,j-nib+2);
    assert (j-nib+1 < 6);
    G1[j-nib] +=   _subs[ib + i]._k * g1[2*i];
    G1[j-nib+1] += _subs[ib + i]._k * g1[2*i+1];
    if (_subs[ib + i]._k != 1.) {
      assert (j-nib+3 < 6);
      maxj = MAX(maxj,j-nib+4);
      G1[j-nib+2] += (1.-_subs[ib + i]._k) * g1[2*i];
      G1[j-nib+3] += (1.-_subs[ib + i]._k) * g1[2*i+1];
    }
  }
  
  return maxj;
}



//---------------------------------------------------------------------------------------
// SPEED: not taking advantage of symmetry anywhere!!

void Keeper::_take2Key0(const int ib, const double* g2, double res, const double c) {
  assert(ib <= _block_size-2);
  
  int row = ib, i;
  
  make2Block(g2,g2,c);
  add2Block(_main, row);
  
  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<2; ++i, ++gPtr)
    *gPtr += g2[i]*res;
}

void Keeper::_take4Key0(const int ib, const double* g2, double res, const double c) {
  assert(ib <= _block_size-4);

  int row = ib, i;

  make4Block(g2,g2,c);
  add4Block(_main, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<4; ++i, ++gPtr)
    *gPtr += g2[i]*res;
}

void Keeper::_take6Key0(const int ib, const double* g2, double res, const double c) {
  assert(ib <= _block_size-6);

  int row = ib, i;

  make6Block(g2,g2,c);
  add6Block(_main, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<6; ++i, ++gPtr)
    *gPtr += g2[i]*res;
}

void Keeper::_take2Key1(const int ib, const double* g1, double res, const double c) {
  assert(ib <= _block_size-2);
  int row = (_num_blocks-1) * _block_size + ib, i;

  make2Block(g1,g1,c);
  add2Block(_main, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<2; ++i, ++gPtr)
    *gPtr += g1[i]*res;
}



void Keeper::_take4Key1(const int ib, const double* g1, double res, const double c) {
  assert(ib <= _block_size-4);
  int row = (_num_blocks-1) * _block_size + ib, i;

  make4Block(g1,g1,c);
  add4Block(_main, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<4; ++i, ++gPtr)
    *gPtr += g1[i]*res; 
}

void Keeper::_take6Key1(const int ib, const double* g1, double res, const double c) {
  assert(ib <= _block_size-6);
  int row = (_num_blocks-1) * _block_size + ib, i;

  make6Block(g1,g1,c);
  add6Block(_main, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<6; ++i, ++gPtr)
    *gPtr += g1[i]*res;
}

void Keeper::_take2Keyb(const int b, const int ib, const double* g1, double res, const double c) {
  assert(ib <= _block_size-2);
  int row = b * _block_size + ib, i;

  make2Block(g1,g1,c);
  add2Block(_main, row);

  // add c * r to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<2; ++i, ++gPtr)
    *gPtr += g1[i]*res; 
}


void Keeper::_take4Keyb(const int b, const int ib, const double* g1, double res, const double c) {
  assert(ib <= _block_size-4);
  int row = b * _block_size + ib, i;

  make4Block(g1,g1,c);
  add4Block(_main, row);

  // add c * r to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<4; ++i, ++gPtr)
    *gPtr += g1[i]*res; 
}


void Keeper::_take2(const int b, const int ib, const double* g1, const double* g2, double res, const double c) {
  assert(b < _num_blocks-1);
  assert(ib <= _block_size-2);

  int row = b*_block_size + ib, i;

  make2Block(g1,g1,c);
  add2Block(_main, row);
  make2Block(g2,g2,c);
  add2Block(_main, row + _block_size);

  make2Block(g1,g2,c);
  add2Block(_right, row);
  make2Block(g2,g1,c);
  add2Block(_left, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<2; ++i, ++gPtr)
    *gPtr += g1[i]*res;
  gPtr =  _g + row +_block_size;
  for (i=0; i<2; ++i, ++gPtr) 
    *gPtr += g2[i]*res;
}

void Keeper::_take4(const int b, const int ib, const double* g1, const double* g2, double res, const double c) {
  assert(b < _num_blocks-1);
  assert(ib <= _block_size-4);

  int row = b*_block_size + ib, i;

  make4Block(g1,g1,c);
  add4Block(_main, row);
  make4Block(g2,g2,c);
  add4Block(_main, row + _block_size);

  make4Block(g1,g2,c);
  add4Block(_right, row);
  make4Block(g2,g1,c);
  add4Block(_left, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<4; ++i, ++gPtr)
    *gPtr += g1[i]*res;
  gPtr =  _g + row +_block_size;
  for (i=0; i<4; ++i, ++gPtr) 
    *gPtr += g2[i]*res;
}



void Keeper::_take6(const int b, const int ib, const double* g1, const double* g2, double res, const double c) {
  assert(b < _num_blocks-1);
  assert(ib <= _block_size-6);

  int row = b*_block_size + ib, i;

  make6Block(g1,g1,c);
  add6Block(_main, row);
  make6Block(g2,g2,c);
  add6Block(_main, row + _block_size);

  make6Block(g1,g2,c); // SPEED: these two blocks are the same! Ditto take4, take2
  add6Block(_right, row);
  make6Block(g2,g1,c);
  add6Block(_left, row);

  // add r * grad to right 
  res *= c;
  double* gPtr = _g + row;
  for (i=0; i<6; ++i, ++gPtr)
    *gPtr += g1[i]*res;

  gPtr =  _g + row +_block_size;
  for (i=0; i<6; ++i, ++gPtr)
    *gPtr += g2[i]*res;
}

void Keeper::_takejKey0(int j, const int ib, const double* g2, double res, const double c) {
  if (j==2)
    _take2Key0(ib,g2,res,c);
  else if (j==4)
    _take4Key0(ib,g2,res,c);
  else if (j==6)
    _take6Key0(ib,g2,res,c);
  else(assert(0));
}

void Keeper::_takejKey1(int j, const int ib, const double* g1, double res, const double c) {
  if (j==2)
    _take2Key1(ib,g1,res,c);
  else if (j==4)
    _take4Key1(ib,g1,res,c);
  else if (j==6)
    _take6Key1(ib,g1,res,c);
  else(assert(0));
}

void Keeper::_takej(int j, const int b, const int ib, const double* g1, const double* g2, double res, const double c) {
  if (j==2)
    _take2(b,ib,g1,g2,res,c);
  else if (j==4)
    _take4(b,ib,g1,g2,res,c);
  else if (j==6)
    _take6(b,ib,g1,g2,res,c);
  else(assert(0));
}

void Keeper::_takejKeyb( const int b, int j, const int ib, const double* g1, double res, const double c) {
  if (j==2)
    _take2Keyb(b,ib,g1,res,c);
  else if (j==4)
    _take4Keyb(b,ib,g1,res,c);
  else if (j==6)
    assert(0);
  else (assert(0));
    }


// SPEED: should also use symmetry for these operations

void Keeper::make2Block(const double* a, const double* b, const double co) {
  int r,c,off=0;
  for (r=0; r<2; ++r)
    for (c=0; c<2; ++c,++off)
      _D[off] = co * a[r]*b[c];      
}

void Keeper::make4Block(const double* a, const double* b, const double co) {
  int r,c,off=0;
  for (r=0; r<4; ++r)
    for (c=0; c<4; ++c,++off)
      _D[off] = co * a[r]*b[c];      
}

void Keeper::make6Block(const double* a, const double* b, const double co) {
  int r,c,off=0;
  for (r=0; r<6; ++r)
    for (c=0; c<6; ++c,++off)
      _D[off] = co * a[r]*b[c];      
}

void Keeper::add2Block(double* strip, const int row) {
  int dOff=0, r,c, sOff=row*_WIDTH + _DIAG;
  for (r=0; r<2; ++r, sOff += _WIDTH-1)
    for (c=0; c<2; ++c, ++dOff)
      strip[sOff+c] += _D[dOff];
}

void Keeper::add4Block(double* strip, const int row) {
  int dOff=0, r,c, sOff=row*_WIDTH + _DIAG;
  for (r=0; r<4; ++r, sOff += _WIDTH-1)
    for (c=0; c<4; ++c, ++dOff)
      strip[sOff+c] += _D[dOff];
}

void Keeper::add6Block(double* strip, const int row) {
  int dOff=0, r,c, sOff=row*_WIDTH + _DIAG;
  for (r=0; r<6; ++r, sOff += _WIDTH-1)
    for (c=0; c<6; ++c, ++dOff)
      strip[sOff+c] += _D[dOff];
}
  

Keeper::~Keeper() {
  delete[] _main; delete[] _left; delete[] _right; delete[] _g;
  delete[] _subs;
}


void Keeper::print() {
  int r,c,off;
  off=0;
  printf("main:\n");
  for (r=0; r<_numVar; r++) {
    for (c=0; c<_WIDTH; c++, ++off)
      printf("%.2f ",_main[off]);
    printf("\n");
  }
  printf("\n");

  off=0;
  printf("left:\n");
  for (r=0; r<_offNumVar; r++) {
    for (c=0; c<_WIDTH; c++, ++off)
      printf("%.2f ",_left[off]);
    printf("\n");
  }
  printf("\n");

  off=0;
  printf("right:\n");
  for (r=0; r<_offNumVar; r++) {
    for (c=0; c<_WIDTH; c++, ++off)
      printf("%.2f ",_right[off]);
    printf("\n");
  }
  printf("\n");

  printf("g:\n");
  for (r=0; r<_numVar; r++) 
      printf("%.2f\n",_g[r]);
  printf("\n");
}


void Keeper::printFull(FILE* fp) {
  int rb,cb,r,c,ir,ic;
  for (rb=0; rb<_num_blocks; rb++) {     // iterate over rows of blocks
    for (r=0; r<_block_size; r++) {      // iterate over rows
      for (cb=0; cb<_num_blocks; cb++) { // iterate over column blocks
	for (c=0; c<_block_size; ++c) {    // print one row of one block
	  ic = c - r + _DIAG;
	  ir = rb*_block_size + r;
	  if (rb==cb) {               // main strip
	    if (ic<0 || ic >= _WIDTH)
	      fprintf(fp,"0 ");
	    else
	      fprintf(fp,"%.10f ", _main[ir*_WIDTH + ic]);
	  }
	  else if (cb-rb==1) {        // right strip
	    if (ic<0 || ic >= _WIDTH)
	      fprintf(fp,"0 ");
	    else
	      fprintf(fp,"%.10f ", _right[ir*_WIDTH + ic]);
	  }
	  else if (rb-cb==1) {        // left strip
	    if (ic<0 || ic >= _WIDTH)
	      fprintf(fp,"0 ");
	    else {
	      ir -= _block_size;
	      fprintf(fp,"%.10f ", _left[ir*_WIDTH + ic]);
	    }
	  }
	  else {                      // zeros
	    fprintf(fp,"0 ");
	  }
	}
      }
      fprintf(fp,"\n");
    }
  }

  fprintf(fp,"\n");
  
}

void Keeper::fillMatrix(double* mat, int stride) const {
  double *ptr=_main;
  int r,c, i, br, bw=0;

  // main strip
  for (r=0,br=0; r<_numVar; ++r, ++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	mat[r*stride + i+bw] = *ptr;
    }
  }

  // right strip
  ptr = _right;
  bw=0;
  for (r=0,br=0; r<_offNumVar; ++r, ++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	mat[r*stride + i + bw + _block_size] = *ptr;
    }
  }

  // left strip
  ptr = _left;   bw=0;
  for (r=0,br=0; r<_offNumVar; ++r, ++br) {
    if (br==_block_size) {
      bw+=br; br=0;
    }
    for (c=0; c<_WIDTH; ++c, ++ptr) {
      i = (r-bw) + c - _DIAG;
      if (i>=0 && i<_block_size)
	mat[(r+_block_size)*stride + i+bw] = *ptr;
    }
  }
  
  
}


void Keeper::printE(FILE* fp) const {
  for (int r=0; r<_numVar; r++) 
    fprintf(fp,"%.5f\n",_g[r]);
}

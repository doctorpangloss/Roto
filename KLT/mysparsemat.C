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



#include "mysparsemat.h"

MySparseMat::~MySparseMat() {
  if (_A) delete[] _A;
  delete[] _ind;
  delete[] _pckdA; 
  if (_pckdAT)
    delete[] _pckdAT;
}

void MySparseMat::init(const int nv, const bool T) {
  assert(_A==NULL);
  _T = T;
  _nv = nv;
  if (_A) delete[] _A;
  _A = new double[nv*nv];
  assert(_A);
  memset(_A,0,nv*nv*sizeof(double));
  _pckdA = new double[nv*10];
  memset(_pckdA, 0, 10*nv*sizeof(double)); // SPEED?
  if (T) {
    _pckdAT = new double[nv*10];
    memset(_pckdAT, 0, 10*nv*sizeof(double)); // SPEED?
  }
  else 
    _pckdAT = NULL;

  _ind = new int[nv];
  memset(_ind,0,6*sizeof(int));
  int ind = 0;
  for (int i=6; i<_nv; ++i) {
    if (i%2==0)
      ind+=2;
    _ind[i]=ind;
  }
      
}

void MySparseMat::pack() {
  _packed = true;
  if (_T) {
    double *_Aptr = _A, *_pAptr = _pckdA;
    int i,j,count, up, tmp;
    for (i=0; i<_nv; ++i,_pAptr+=10, _Aptr+=_nv) {
      count = MIN(10, _nv-_ind[i]);
      memcpy(_pAptr, _Aptr+_ind[i], count * sizeof(double));
      tmp = _ind[i]+count;
      for (j=_ind[i],up=_ind[i]*10; j<tmp; j++, up+=10) {  // SPEED!  get rid of multiply...
	if (i>=_ind[j]) {
	  _pckdAT[up +i -_ind[j]] = *(_Aptr+j);
	}
      }
    }
  }
    
  else {
    double *_Aptr = _A, *_pAptr = _pckdA;
    int i,count;
    for (i=0; i<_nv; ++i,_pAptr+=10, _Aptr+=_nv) {
      count = MIN(10, _nv-_ind[i]);
      memcpy(_pAptr, _Aptr+_ind[i], count * sizeof(double));
    }
  }

  /*for (int i=0; i<_nv; i++) {
    for (int j=0; j<10; j++)
    printf("%.2f ", _pckdA[i*10 + j]);
    printf("\n");
    }*/
    
  delete[] _A; _A = NULL;
}

void MySparseMat::mult(const double* x, double* r) const {
  double *ptr = _pckdA;
  const double *xptr;
  int *iptr = _ind;
  int n,k, ul;
  double res;
  for (n=0; n<_nv; ++n,++r) {
    res=0;
    xptr = x + (*iptr);
    ul = MIN(10, _nv-(*iptr));
    ++iptr;
    for (k=0; k<10; ++k,++ptr,++xptr)
      if (k<ul)
	res += (*ptr) * (*xptr);
    (*r) += res;
  }
}

void MySparseMat::multT(const double* x, double* r) const {
  double *ptr = _pckdAT;
  const double *xptr;
  int *iptr = _ind;
  int n,k, ul;
  double res;
  for (n=0; n<_nv; ++n,++r) {
    res=0;
    xptr = x + (*iptr);
    ul = MIN(10, _nv-(*iptr));
    ++iptr;
    for (k=0; k<10; ++k,++ptr,++xptr)
      if (k<ul)
	res += (*ptr) * (*xptr);
    (*r) += res;
  }
}

void MyCongMat::matVecMult(const double x[], double r[]) const {
  memset(r,0,_numVec*l2*sizeof(double));
  int t,i1=0, j1;
  for (t=0; t<_numVec; ++t,i1+=l2) {
    _A[t].mult(x+i1,r+i1);
  }      
    
  i1=0; j1=l2;
  int i2=l2, j2=0;
  for (t=0; t<_numVec-1; ++t,i1+=l2,j1+=l2,i2+=l2,j2+=l2) { // SPEED: increment ptrs
    _Al[t].mult(x+i1,r+j1);  // SPEED: make a mult2 that does 2 at a time
    _Al[t].multT(x+i2,r+j2);
  }

  
  //printf("result:\n");
  /*
  int rr;
  FILE *fp = fopen("/homes/gws/aseem/r.txt","w");
  for (rr=0; rr<_numVec*l2; rr++) 
    fprintf(fp,"%.5f\n",r[rr]);
  //printf("\n");
  fclose(fp);
  
  fp = fopen("/homes/gws/aseem/x.txt","w");
  for (rr=0; rr<_numVec*l2; rr++) 
    fprintf(fp,"%.5f\n",x[rr]);
  //printf("\n");
  fclose(fp);

  fp = fopen("/homes/gws/aseem/Z.txt","w");
  printExpanded(fp);
  fclose(fp);
  
  exit(0);
  */
  //for (t=0; t<_numVec*l2; t++)
  //printf("%7.2f ",r[t]);
  //printf("\n\n");
}


void MyCongMat::printExpanded(FILE* fp) {
//  int dim = _numVec*l2;
//  double* data = new double[dim*dim];
//  memset(data, 0, dim*dim*sizeof(double));
//  int i, t, offset;
//  for (i=0; i<_numVec; i++) {
//    offset = i*l2;
//    _A[i].drop(data + offset*dim + offset, dim); 
//  }
//
//  int i1=0, j1=l2;
//  int i2=l2, j2=0;
//  for (t=0; t<_numVec-1; ++t,i1+=l2,j1+=l2,i2+=l2,j2+=l2) { 
//    _Al[t].drop(data + dim*j1 + i1, dim);
//    _Al[t].dropT(data + dim*j2 + i2, dim);
//  }
//
//  for (t=0; t<dim; t++) {
//    for (i=0; i<dim; i++) {
//      fprintf(fp,"%.16f ", data[t*dim+i]);
//      //assert(data[t*dim+i] == data[i*dim+t]);
//    }
//    fprintf(fp,"\n");
//  }
//
//  delete[] data;
    return;
}

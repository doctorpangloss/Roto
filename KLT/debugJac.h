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



#ifndef DEBJAC
#define DEBJAC

#include "stdio.h"
#include <vector>

struct DebElem {
  int c;
  double val;
};

class DebugJacobian : public vector<DebElem> {
 public:
  void add(int c, double val) {
    DebElem de;
    de.c = c; de.val = val;
    push_back(de);
  }

  void output(FILE* fp, const int i) {
    for (int j=0; j<size(); ++j)
      fprintf(fp,"%d %d %.10f\n",i+1,(*this)[j].c+1,(*this)[j].val);
  }

  void outputCons(FILE* fp2) {
    fprintf(fp2, "%.10f\n", cons);
  }

  void outputRes(FILE* fp2) {
    fprintf(fp2, "%.10f\n", res);
  }

  double cons;
  double res;

};



class DebugJacobians : public vector<DebugJacobian> {
 public:

  void output() {
    FILE* fp = fopen("/homes/gws/aseem/jac.dat", "w");
    FILE* fp2 = fopen("/homes/gws/aseem/cons.txt", "w");
    FILE* fp3 = fopen("/homes/gws/aseem/res.txt", "w");

    for (int i=0; i<size(); ++i) {
      (*this)[i].output(fp, i);
      (*this)[i].outputCons(fp2);
      (*this)[i].outputRes(fp3);
    }

    fclose(fp);
    fclose(fp2);
    fclose(fp3);
  }
};

#endif

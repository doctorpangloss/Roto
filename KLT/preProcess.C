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



#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "pnmio.h"
#include "klt.h"


static const float pyramid_sigma_fact = 0.9;

int main(int argc, char *argv[]) {
  unsigned char *img;
  int ncols,nrows;
  KLT_TrackingContext tc;
  FILE *fp, *fp1, *fp2, *fp3;
  char *basename;
  char name1[100], name2[100], name3[100];

  assert(argc==2);
  fp = fopen(argv[1], "r");
  img = pgmReadFile(argv[1], NULL, &ncols, &nrows);
  fclose(fp);

  basename = strsep(argv+1,".");
  strcpy(name1,basename);
  strcat(name1,".klt");
  strcpy(name2,basename);
  strcat(name2,".gxklt");
  strcpy(name3,basename);
  strcat(name3,".gyklt");

  fp1 = fopen(name1,"w"); assert(fp1);
  fp2 = fopen(name2,"w"); assert(fp2);
  fp3 = fopen(name3,"w"); assert(fp3);

  tc = KLTCreateTrackingContext();
  
  KLTCreateFile(img, ncols, nrows, tc, fp1,fp2,fp3);
  fclose(fp1); fclose(fp2); fclose(fp3);
  return 0;
}
    

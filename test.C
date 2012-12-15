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



#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <qapplication.h>
#include <qimage.h>
#include <qdatetime.h>
#include "KLT/klt.h"

int main( int argc, char *argv[] ) {

  QApplication a( argc, argv );
  //int aFrame = 10, bFrame = 11; // HERE
  KLT_TrackingContext TC;
  //TC.nPyramidLevels = 4;
  //TC.usePseudo=1;
  //TC.useSmoothness=0;
  //TC.pinLast=0;

  QImage ims[15];
  assert(ims[0].load("framedata/amira-queen/frame001.png")); 
  assert(ims[1].load("framedata/amira-queen/frame002.png"));
  assert(ims[2].load("framedata/amira-queen/frame003.png")); 
  assert(ims[3].load("framedata/amira-queen/frame004.png"));
  assert(ims[4].load("framedata/amira-queen/frame005.png")); 
  assert(ims[5].load("framedata/amira-queen/frame006.png"));
  assert(ims[6].load("framedata/amira-queen/frame007.png")); 
  assert(ims[7].load("framedata/amira-queen/frame008.png"));
  assert(ims[8].load("framedata/amira-queen/frame009.png")); 
  assert(ims[9].load("framedata/amira-queen/frame010.png"));
  assert(ims[10].load("framedata/amira-queen/frame011.png")); 
  assert(ims[11].load("framedata/amira-queen/frame012.png"));
  assert(ims[12].load("framedata/amira-queen/frame013.png")); 
  assert(ims[13].load("framedata/amira-queen/frame014.png"));
  assert(ims[14].load("framedata/amira-queen/frame015.png")); 

  
  const KLT_FullCPyramid* pyrms[15];
  pyrms[0] = new KLT_FullCPyramid(ims[0],&TC);
  pyrms[1] = new KLT_FullCPyramid(ims[1],&TC);
  pyrms[2] = new KLT_FullCPyramid(ims[2],&TC);
  pyrms[3] = new KLT_FullCPyramid(ims[3],&TC);
  pyrms[4] = new KLT_FullCPyramid(ims[4],&TC);
  pyrms[5] = new KLT_FullCPyramid(ims[5],&TC);
  pyrms[6] = new KLT_FullCPyramid(ims[6],&TC);
  pyrms[7] = new KLT_FullCPyramid(ims[7],&TC);
  pyrms[8] = new KLT_FullCPyramid(ims[8],&TC);
  pyrms[9] = new KLT_FullCPyramid(ims[9],&TC);
  pyrms[10] = new KLT_FullCPyramid(ims[10],&TC);
  pyrms[11] = new KLT_FullCPyramid(ims[11],&TC);
  pyrms[12] = new KLT_FullCPyramid(ims[12],&TC);
  pyrms[13] = new KLT_FullCPyramid(ims[13],&TC);
  pyrms[14] = new KLT_FullCPyramid(ims[14],&TC);


  int numP,i, res, elapsed;
  FILE* fp = fopen("dummyfile","r");
  assert(fp);
  fread(&numP,sizeof(int),1,fp);
  Vec2f *Plocs = new Vec2f[numP];
  fread(Plocs,sizeof(Vec2f),numP,fp);
  fclose(fp);
  printf("%d points\n",numP);
  
  fp = fopen("Z","r");
  assert(fp);
  Vec2f *Z = new Vec2f[14*numP];
  res = fread(Z,sizeof(Vec2f),14*numP,fp);
  assert(res==14*numP);
  fclose(fp);
  
  /*printf("Plocs:\n");
  for (i=0; i<numP; i++)
    printf("%f %f, ",Plocs[i].x(), Plocs[i].y());
  printf("\n\n");

  printf("beginning Z:\n");
  for (i=0; i<3*numP; i++)
    printf("%f %f, ",Z[i].x(), Z[i].y());
    printf("\n\n");*/

  QTime time;
  time.start();

  res = TC.longCurveTrack(pyrms, 14, numP, 7, -4, 4, 2, Plocs, Z);

  elapsed = time.elapsed();
  printf("resolution %d\n",res);
  //printf("final Z:\n");
  //for (i=0; i<numP; i++)
  //printf("%f %f, ",Z[i].x(), Z[i].y());
  //printf("\n\n");
  printf("%d milliseconds elapsed\n", elapsed);
 
  return 1;
}
  

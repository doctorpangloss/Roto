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
#include "pnmio.h"
#include "klt.h"


int main(int argc, char *argv[]) {

 KLT_TrackingContext tc;
 int startFrame,endFrame, n;
 char fileRoot[100] = "/big/kitty1/cl_raw720x576-";
 char numberName[100], finalName[100], fnamein[100],fnameout[100];
 FILE *fp1=NULL,*fp2=NULL,*fp3=NULL;
 int nFeatures = 150, indx, ncols, nrows;
 unsigned char *img1, *img2;
 float xloc,yloc,xOut,yOut,tvx,tvy;
 _KLT_Pyramid b,bx,by;
 AseemFeatureList fl;
 
 assert(argc==4);
 tc = KLTCreateTrackingContext();
 startFrame = atoi(argv[1]);
 endFrame = atoi(argv[2]);
 
 img1 = pgmReadFile("/big/kitty1/cl_raw720x576-0000.pgm", NULL, &ncols, &nrows);

 fl = aseemReadFeatures(argv[3]);
 aseemWriteFeaturesToPPM(fl,ncols,nrows,img1,"/big/kitty1/feat0000.ppm");

 img2 = (unsigned char *) malloc(ncols*nrows*sizeof(unsigned char));
 free(img1); img1=NULL;
 
 /* set up initial condition */
 sprintf(numberName, "%s%.4d",fileRoot,n);
 sprintf(finalName, "%s.klt",numberName);
 fp1 = fopen(finalName,"r"); assert(fp1); 
 sprintf(finalName, "%s.gxklt",numberName);
 fp2 = fopen(finalName,"r"); assert(fp2);
 sprintf(finalName, "%s.gyklt",numberName);
 fp3 = fopen(finalName,"r"); assert(fp3);
 KLT_StorePyramids(tc,fp1,fp2,fp3);
 fclose(fp1); fclose(fp2); fclose(fp3);
 
 for (n=startFrame+1; n<=endFrame; n++) {
     sprintf(numberName, "%s%.4d",fileRoot,n);
     sprintf(finalName, "%s.klt",numberName);
     fp1 = fopen(finalName,"r");  assert(fp1);
     b = _KLTReadPyramid(fp1);
     fclose(fp1);
     sprintf(finalName, "%s.gxklt",numberName);
     fp2 = fopen(finalName,"r"); assert(fp2);
     bx = _KLTReadPyramid(fp2);
     fclose(fp1);
     sprintf(finalName, "%s.gyklt",numberName);
     fp3 = fopen(finalName,"r"); assert(fp3);
     by = _KLTReadPyramid(fp3);
     fclose(fp1);
     sprintf(fnamein, "/big/kitty1/cl_raw720x576-%.4d.pgm", n);
     pgmReadFile(fnamein, img2, &ncols, &nrows);
     
     for (indx = 0 ; indx < fl->num; indx++)  {
	 xloc = fl->features[indx].x;
	 xOut = xloc + (tvx=fl->features[indx].vx);
	 yloc = fl->features[indx].y;
	 yOut = yloc + (tvy=fl->features[indx].vy);
	 if (n != startFrame+1) {
	     xOut += .5*fl->features[indx].ax;
	     yOut += .5*fl->features[indx].ay;
	 }
	 KLTPreTrackAFeature(tc,b,bx,by,xloc,yloc,&xOut,&yOut);
	 fl->features[indx].vx = xOut - fl->features[indx].x;
	 fl->features[indx].vy = yOut - fl->features[indx].y;
	 fl->features[indx].ax = fl->features[indx].vx - tvx;
	 fl->features[indx].ay = fl->features[indx].vy - tvy;
	 fl->features[indx].x = xOut;
	 fl->features[indx].y = yOut;
	 printf("v: %f %f, a: %f %f\n",fl->features[indx].vx,fl->features[indx].vy,
		fl->features[indx].ax,fl->features[indx].ay);
     }

     _KLTFreePyramid(tc->pyramid_last);
     _KLTFreePyramid(tc->pyramid_last_gradx);
     _KLTFreePyramid(tc->pyramid_last_grady);
     tc->pyramid_last = b;
     tc->pyramid_last_gradx = bx;
     tc->pyramid_last_grady = by;
     
     printf("did frame %d\n",n);
     sprintf(fnameout, "/big/kitty1/feat%.4d.ppm", n);
     aseemWriteFeaturesToPPM(fl,ncols,nrows,img2,fnameout);
 }
 
 
}
















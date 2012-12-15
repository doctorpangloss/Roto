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



//
// bitmap.h
//
// header file for MS bitmap format
//
//

#ifndef BITMAP_H
#define BITMAP_H

#include <stdio.h>
#include <string.h>

#define BMP_BI_RGB        0L

typedef unsigned short	BMP_WORD; 
typedef unsigned int	BMP_DWORD; 
typedef int				BMP_LONG; 
 
typedef struct { 
	BMP_WORD	bfType; 
	BMP_DWORD	bfSize; 
	BMP_WORD	bfReserved1; 
	BMP_WORD	bfReserved2; 
	BMP_DWORD	bfOffBits; 
} BMP_BITMAPFILEHEADER; 
 
typedef struct { 
	BMP_DWORD	biSize; 
	BMP_LONG	biWidth; 
	BMP_LONG	biHeight; 
	BMP_WORD	biPlanes; 
	BMP_WORD	biBitCount; 
	BMP_DWORD	biCompression; 
	BMP_DWORD	biSizeImage; 
	BMP_LONG	biXPelsPerMeter; 
	BMP_LONG	biYPelsPerMeter; 
	BMP_DWORD	biClrUsed; 
	BMP_DWORD	biClrImportant; 
} BMP_BITMAPINFOHEADER; 

// global I/O routines
extern unsigned char* readBMP( char* fname, int& width, int& height );
extern void writeBMP( char* iname, int width, int height, unsigned char* data ); 

#endif


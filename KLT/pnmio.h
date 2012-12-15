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



/*********************************************************************
 * pnmio.h
 *********************************************************************/

#ifndef _PNMIO_H_
#define _PNMIO_H_

#include <stdio.h>

/**********
 * With pgmReadFile and pgmRead, setting img to NULL causes memory
 * to be allocated
 */

/**********
 * used for reading from/writing to files
 */
unsigned char* pgmReadFile(
  char *fname,
  unsigned char *img,
  int *ncols, 
  int *nrows);
void pgmWriteFile(
  char *fname,
  unsigned char *img,
  int ncols,
  int nrows);
void ppmWriteFileRGB(
  char *fname,
  unsigned char *redimg,
  unsigned char *greenimg,
  unsigned char *blueimg,
  int ncols,
  int nrows);

/**********
 * used for communicating with stdin and stdout
 */
unsigned char* pgmRead(
  FILE *fp,
  unsigned char *img,
  int *ncols, int *nrows);
void pgmWrite(
  FILE *fp,
  unsigned char *img,
  int ncols,
  int nrows);
void ppmWrite(
  FILE *fp,
  unsigned char *redimg,
  unsigned char *greenimg,
  unsigned char *blueimg,
  int ncols,
  int nrows);

#endif

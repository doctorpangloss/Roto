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



/**************************************************************************
 *  interface for texturized strokes. 
 *  
 *  author: Siobhan Quinn 
 *  date  : 04-21-2002
 *  ***********************************************************************/

#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#include <qimage.h>
#include <stdio.h> 
#include <assert.h>
//#include "ImageLib/Image.h"
//#include "ImageLib/FileIO.h"
#include "hertzStrokes/bitmap.h"

#define NUM_TEXTURES 6

struct Pixel ; 

enum Texture_Type { stroke, paper, none }; 
enum fileType { BMP, PNG } ; // add more later 
enum Textures { CHARCOAL, PASTEL, PENCIL, PENCIL_STROKE, CHARCOAL_STROKE, PASTEL_STROKE, PEN  } ; 

class Texture 
{ 

 public: 

  static const int Texture_Widths[NUM_TEXTURES]; 
  static const int Texture_Heights[NUM_TEXTURES];
  static const float Texture_Radii[NUM_TEXTURES]; 
  static char* Texture_FNames[NUM_TEXTURES]; 
  static const Texture_Type Texture_Types[NUM_TEXTURES]; 

  Texture() ;

  // loads the texture 
  Texture( fileType type, GLuint name, char* fileName, int width, int height ); 

  //  Texture( fileType type, char* fileName, int width, int height, int img_width , int img_height ); 
  void TextureBMP( GLuint name , char* fileName, int width , int height ) ; 
  void TexturePNG( GLuint name, char* fileName, int width, int height);
  
 private: 
  
  GLuint _name ;
  int _width; int _height;          // the width/height of the paper texture 

};


#endif 

  




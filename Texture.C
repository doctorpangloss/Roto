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



/********************************************************************
 * Texture.C 
 * loads the paper texture
 * 
 *******************************************************************/ 
#include <stdlib.h>
#include "Texture.h" 

const int Texture::Texture_Widths[NUM_TEXTURES] = { 1024, 1024, 1024, 1024, 2048, 2048 } ; 
const int Texture::Texture_Heights[NUM_TEXTURES] = { 1024, 1024, 1024, 16, 64, 64  } ; 

// in case we don't want to consider the full width of the stroke texture
const float Texture::Texture_Radii[NUM_TEXTURES] = { 1.0, 1.0, 1.0, 0.85, 0.50, 0.50 } ; 

const Texture_Type Texture::Texture_Types[NUM_TEXTURES] = { paper, paper, paper, stroke, stroke, stroke }; 
char* Texture::Texture_FNames[NUM_TEXTURES] = {  "textures/charcoal_cropped.png", 
						 "textures/pastel_cropped.png" , 
						 "textures/pencil_cropped.png" , 
						 "textures/pencil_stroke_1024_16.png",
                                                 "textures/charcoal_stroke_2048_64.png", 
                                                 "textures/pastel_stroke_2048_64.png"};
						 
Texture::Texture() 
{ 
  _width = 0 ; 
  _height = 0 ; 
  _name = 0 ; 
  return ; 
} 

Texture::Texture( fileType ftype, GLuint name,  char* fileName, int width, int height ) 
{

  printf("texture constructor called\n"); 

  _width = width ; 
  _height = height ;

  printf("Texture Constructor: width %d, height %d\n", width, height); 
  switch ( ftype ){ 
  case BMP :
    TextureBMP( name,  fileName, width, height) ; 
    break ; 
  case PNG:
    TexturePNG( name, fileName, width, height);
    break;
  default: 
    printf("Bad image type given\n") ; 
    exit(0) ; 
    break ;
  }

  return ; 
}

 

void Texture::TextureBMP( GLuint name, char* fileName, int width, int height)
{
  unsigned char* dataRGB ;
  unsigned char* data ; 
  _name = name ; 

  printf("file_name %s\n", fileName ) ; 

  // open texture data 
  dataRGB = readBMP( fileName, width, height);  
  if ( dataRGB == NULL )
  {
    printf("bitmap not read\n"); 
    assert(false) ; 
  } 

  // make the texture contain all the given texture information in the alpha
  // channel only. with this we will only map the alpha channel - and not 
  // affect the rgb channels in the textured area. 
  
  data = new unsigned char [ 4 * width * height ] ; 
  assert( data ) ;   

  // calculate the intensity for each pix in the rgb data and copy to the rgba data. 
  for ( int row = 0 ; row < height ; row ++ ) { 
    
    for ( int col = 0 ; col < width ; col ++ ) {

        unsigned char luminance = 0.299 *dataRGB[ 3*( row*width + col) + 0] + 
	                          0.587 *dataRGB[ 3*( row*width + col) + 1] + 
			          0.114 *dataRGB[ 3*( row*width + col) + 2] ;  
	if (luminance < 50  )
	  luminance = 0 ;

	if ( luminance > 205 ) 
	  luminance = 255 ; 

	
	       	data [ 4*( row*width + col ) + 0 ] = 255; 
		data [ 4*( row*width + col ) + 1 ] = 255 ; 
		data [ 4*( row*width + col ) + 2 ] = 255; 
	 	data [ 4*( row*width + col ) + 3 ] = 255-luminance; 

    } 
  }


  glBindTexture( GL_TEXTURE_2D, _name ) ; 
  glPixelStorei( GL_UNPACK_ALIGNMENT , 1 ) ; 



  /*** Have to make sure that blending with colors is happening ***/ 
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE , GL_MODULATE ) ;  
  //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE , GL_DECAL ) ;  
 
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 

  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // length of stroke 
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);  // cross section of stroke

  // when texture area is small, bilinear filter the closest mipmap
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    
  // when texture area is large, bilinear filter the first mipmap
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height, 
		0, GL_RGBA, GL_UNSIGNED_BYTE, 
  		data );

  assert(glGetError() == GL_NO_ERROR);

  delete [] dataRGB ;  
  delete [] data ; 

} 


void Texture::TexturePNG( GLuint name, char* fileName, int width, int height)
{
  
  unsigned char* data ; 
  _name = name ; 

  printf("file_name %s\n", fileName ) ; 

  // open texture data 
  QImage* im = new QImage(fileName);  
  if (!(im->width() == width && im->height() == height)) {
    printf("PNG expected %d %d, got %d %d\n", width, height, 
	   im->width(), im->height());
    exit(0);
  }
  if ( im == NULL )
  {
    printf("PNG not read\n"); 
    assert(false) ; 
  } 

  // make the texture contain all the given texture information in the alpha
  // channel only. with this we will only map the alpha channel - and not 
  // affect the rgb channels in the textured area. 
  
  data = new unsigned char [ 4 * width * height ] ; 
  assert( data ) ;   

  // calculate the intensity for each pix in the rgb data and copy to the rgba data. 
  for ( int row = 0 ; row < height ; row ++ ) { 
    
    for ( int col = 0 ; col < width ; col ++ ) {
        QRgb color = im->pixel(col,row);
        unsigned char luminance = 0.299 * qRed(color) + 
	                          0.587 * qGreen(color) + 
			          0.114 * qBlue(color) ;  
	if (luminance < 50  )
	  luminance = 0 ;

	if ( luminance > 205 ) 
	  luminance = 255 ; 

	
	       	data [ 4*( row*width + col ) + 0 ] = 255; 
		data [ 4*( row*width + col ) + 1 ] = 255 ; 
		data [ 4*( row*width + col ) + 2 ] = 255; 
	 	data [ 4*( row*width + col ) + 3 ] = 255-luminance; 

    } 
  }


  glBindTexture( GL_TEXTURE_2D, _name ) ; 
  glPixelStorei( GL_UNPACK_ALIGNMENT , 1 ) ; 



  /*** Have to make sure that blending with colors is happening ***/ 
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE , GL_MODULATE ) ;  
  //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE , GL_DECAL ) ;  
 
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 

  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // length of stroke 
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);  // cross section of stroke

  // when texture area is small, bilinear filter the closest mipmap
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    
  // when texture area is large, bilinear filter the first mipmap
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height, 
		0, GL_RGBA, GL_UNSIGNED_BYTE, 
  		data );

  assert(glGetError() == GL_NO_ERROR);

  delete im ;  
  delete [] data ; 

} 

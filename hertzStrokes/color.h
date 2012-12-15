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



#ifndef __COLOR_H__
#define __COLOR_H__

//#include "painter.h"
#define NUM_SLICES 20

enum COLORSPACE { RGB_SPACE, XYZ_SPACE, YIQ_SPACE, LAB_SPACE, LUV_SPACE };

#ifndef IMAGE_MODE

//#ifdef _WIN32
#define IMAGE_MODE  GL_RGBA
#define RED_OFFSET 0
#define GREEN_OFFSET 1
#define BLUE_OFFSET 2
#define ALPHA_OFFSET 3

// note: other IMAGE_MODEs are not completely supported
//#else
//#define IMAGE_MODE  GL_ABGR_EXT
//#define RED_OFFSET 3
//#define GREEN_OFFSET 2
//#define BLUE_OFFSET 1
//#define ALPHA_OFFSET 0
//#endif

#endif

inline GLuint makeColor(GLubyte red,GLubyte green, GLubyte blue, GLubyte alpha)
{
  return
    (GLuint(red) << ((3-RED_OFFSET)*8)) |
    (GLuint(green) << ((3-GREEN_OFFSET)*8)) |
    (GLuint(blue) << ((3-BLUE_OFFSET)*8)) |
    (GLuint(alpha) << ((3-ALPHA_OFFSET)*8));
}

inline GLuint makeColor(GLubyte red,GLubyte green, GLubyte blue, float alpha)
{
  return
    (GLuint(red) << ((3-RED_OFFSET)*8)) |
    (GLuint(green) << ((3-GREEN_OFFSET)*8)) |
    (GLuint(blue) << ((3-BLUE_OFFSET)*8)) |
    (GLuint(alpha*0xFF) << ((3-ALPHA_OFFSET)*8));
}

inline GLuint makeColor(GLubyte red,GLubyte green,GLubyte blue)
{
  return makeColor(red,green,blue,GLubyte(0xFF));
}


inline GLubyte getRed(GLuint color) 
{ 
    return (color >> ((3-RED_OFFSET) * 8)) & 0xFF;
}

inline GLubyte getGreen(GLuint color) 
{ 
    return (color >> ((3-GREEN_OFFSET) * 8)) & 0xFF;
}

inline GLubyte getBlue(GLuint color) 
{ 
    return (color >> ((3-BLUE_OFFSET) * 8)) & 0xFF;
}

inline GLubyte getAlpha(GLuint color) 
{ 
    return (color >> ((3-ALPHA_OFFSET) * 8)) & 0xFF;
}




inline float Max(float v1,float v2)
{
  return v1 > v2 ? v1 : v2;
}

inline float Min(float v1,float v2)
{
  return v1 < v2 ? v1 : v2;
}

inline void RGB_To_HSV(GLubyte red, GLubyte green, GLubyte blue,
		       float *h,float *s, float *v)
{
  float r = red/255.0;
  float g = green/255.0;
  float b = blue/255.0;

  float max = Max(r,Max(g,b));
  float min = Min(r,Min(g,b));

  *v = max;
  
  *s = (max != 0) ? ((max-min)/max) : 0;

  if (*s != 0)
    {
      double delta = max-min;
      if (r==max)
	*h=(g-b)/delta;
      else
	if (g==max)
	  *h=2+(b-r)/delta;
	else
	  if (b==max)
	    *h=4+(r-g)/delta;

      *h /= 6;

      if (*h < 0)
	*h += 1;
    }
}
      
inline void HSV_To_RGB(GLubyte * red, GLubyte * green, GLubyte * blue,
		       float h,float s,float v)
{
  if (s == 0)
    {
      *red = *green = *blue = (GLubyte)(v*255);

      return;
    }

  double f,p,q,t,r,g,b;
  int i;

  if (h == 1)
    h = 0;

  h *= 6;

  i = (int)floor(h);

  f = h-i;
  p = v*(1-s);
  q = v*(1-s*f);
  t = v*(1-(s*(1-f)));

  switch(i)
    {
    case 0: r=v; g=t; b=p; break;
    case 1: r=q; g=v; b=p; break;
    case 2: r=p; g=v; b=t; break;
    case 3: r=p; g=q; b=v; break;
    case 4: r=t; g=p; b=v; break;
    case 5: r=v; g=p; b=q; break;
    }

  *red = GLubyte(r*255);
  *green = GLubyte(g*255);
  *blue = GLubyte(b*255);
}


/*
#ifdef _WIN32
inline double drand48()
{
  return double(rand())/RAND_MAX;
}
#endif
*/


inline float luminance(GLubyte red,GLubyte green,GLubyte blue)
{
  //  fprintf(stderr, "r = %d, g= %d, b=%d  l=%f",red,green,blue,
  //	 (30*red + 59*green + 11*blue)/(0xFF*100.0));


  return (30*red + 59*green + 11*blue)/(100.0);
}


inline void RGBtoYIQ(double red, double green, double blue,
		      double &luminance, double &I, double &Q)
{
    // from Foley-Van Dam: NTSC RGB
  luminance = .299*red + .587*green + .114*blue;
  I =         .596*red  - .275*green - .321*blue;
  Q =         .212*red  - .523*green + .311*blue;
}

inline void YIQtoRGB(double luminance, double I, double Q,
		     double & red, double & green, double &blue)
{
  red = luminance + 0.95568806036116*I + 0.61985809445637*Q;
  green = luminance -0.27158179694406*I -0.64687381613840*Q;
  blue = luminance -1.10817732668266*I + 1.70506455991918*Q;
}



#endif




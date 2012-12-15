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



#ifndef __STROKE_HH__
#define __STROKE_HH__

using namespace std;
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include "Point.h"
#include "color.h"
#include "../Texture.h"

typedef long int StrokeNumT;
#define NOSTROKE (StrokeNumT(-1))
#define RAMP_FACTOR 0.6
#define RAMP_FACTOR_TWO 1.2 
#define RAND_OFFSET_WIDTH 0.25
inline StrokeNumT makeStrokeNum(long int num,int) { return num; }
struct _Vertex ; 

struct _CapData { 
  PointR p0; 
  float dx,  dy ; 
  float textureU; 
  bool orientation;
} ;


enum Strip { __left , __middle, __right } ; 

class Stroke
{
private:
  // control polygon
  vector<PointR> _control;

  // limit curve
  vector<PointR> * _limit;  

  // temporary data for subdivision
  vector<PointR> * _temp;

  // points that define the triangle strips to draw.
  // there are 3 adjacent triangle strips in total. 
  // "computeStripPoints" calculates values for these points.  
  vector<Point> left_inner;    // the points that define the left side of the central triangle strip
  vector<Point> left_outer;
  vector<Point> right_inner;   // the points that define the right side of the central triangle strip
  vector<Point> right_outer; 

  vector<Point> inner_cap1 ; 
  vector<Point> outer_cap1 ; 
  vector<Point> inner_cap2 ; 
  vector<Point> outer_cap2 ; 

  _CapData startCap; 
  _CapData endCap; 

  // for stroke textures
  vector<float> textureU ; 
  vector<Point> startCapTexture_inner; 
  vector<Point> startCapTexture_outer; 
  vector<Point> endCapTexture_inner; 
  vector<Point> endCapTexture_outer; 


  // z value
  float _depth;

  // has the limit curve been computed from the current CP?
  bool _computed;  
  
  // have the triangle strip points been computed ? 
  bool _stripComputed ; 

  // number of subdivision levels to use
  int _numLevels; 

  // brush radius
  // if this is changed, then the bounding box will need to be updated
  float _radius;

  // draw round caps at the ends of the stroke?
  bool _cap;

  bool _extendLength;

  // stroke color
  GLuint _color;
  GLuint _adjustedColor;

  // stroke number
  StrokeNumT _num;

  // texture mapping parameters
  bool _useTexture;
  GLuint _textureName ;
  Texture_Type _textureType; 
  int _texture_width, _texture_height; 
  float _texture_radius; 

  float _ustart;
  float _ufreq;

  GLuint _lumAlphaNum;
  GLuint _alphaNum;


  float _tOffset;
  
  static const double  MEAN_FILTER[3] ; 
  

  // private member functions

  //  void discPoint(float x,float y,void (*drawPoint)(int x,int y));

  void drawCap(const _CapData & cap, vector<Point> &  inner_points, vector<Point> &  outer_points, bool orientation ); 

  void subdivideCubicBSpline(vector<PointR> * inputCurve, 
			     vector<PointR> * outputCurve);
  void subdivideFourPoint(vector<PointR> * inputCurve, 
			  vector<PointR> * outputCurve);
  void subdivide(vector<PointR> * inputCurve, 
		 vector<PointR> * outputCurve);
  void computeLimitCurve();
  void scanConvert();

  void drawTriangleStrip( Strip strip) ; 
  
  //  void computeCapPoints ( float cap1_offsets[] , float cap2_offsets, int size ) ; 
  void computeCapPoints( PointR p0, float dx, float dy, float* cap_offsets, bool orientation ); 
  static void blendArray( float arrOut[], float arrIn[], int sizeOut, int sizeIn, int startIndex ) ; 
  static void makeRandomArray( float arrOut[] , int size); 
  static void copyArray( float arrOut[] , float arrIn[], int size ) ;  
  
 

  static GLUquadricObj * qobj;

public:
  Stroke();
  //Stroke(FILE* fp);
  Stroke(float radius,GLuint color,StrokeNumT n);
  Stroke(int cx,int cy,float radius,GLubyte r,GLubyte g,GLubyte b,
	 StrokeNumT n);
  Stroke(int cx,int cy,float radius,GLuint color,StrokeNumT n);
  Stroke(const vector<PointR> & control,float radius,GLuint color,StrokeNumT n);

  Stroke( const Stroke& s);
  Stroke& operator=( const Stroke& s);


  ~Stroke();

  void setUniformRadius(const float t);

  float radius() const { return _radius; }
  float & radius() { return _radius; }
  float depth() const { return _depth; }
  float & depth() { return _depth; }
  bool cap() const { return _cap; }
  bool & cap() { return _cap; }
  bool extendLength() const { return _extendLength; }
  bool & extendLength() { return _extendLength; }
  int & numLevels() { return _numLevels; }
  float toffset() const { return _tOffset; }
  float & toffset() { return _tOffset; }

  vector<PointR> & control() { return _control; }
  const vector<PointR> & control() const { return _control; }  
  StrokeNumT num() const { return _num; }

  GLubyte red() const { return getRed(_color); }
  GLubyte green() const { return getGreen(_color); }
  GLubyte blue() const { return getBlue(_color); }
  GLubyte alpha() const { return getAlpha(_color); }
  GLuint color() const { return _color; }
  GLuint & color() { return _color; }
  GLuint adjustedColor() const { return _adjustedColor; }
  GLuint & adjustedColor() { return _adjustedColor; }

  GLuint & textureName() { return _textureName ; }
  Texture_Type & textureType() { return _textureType ; } 

  float & textureRadius() { return _texture_radius ; } 
  

  int & textureWidth() { return _texture_width ; } 
  int & textureHeight() { return _texture_height ; } 
  

  void clear();
  static void drawLines(vector<PointR> * curve);
  void drawWideLineCurve() ; 
  void drawControlPolygon();
  void drawLineCurve();
  void render();
  void print(FILE * fp = stdout);

  void computeStripPoints( const vector<PointR> *curve) ; 
  

  // Aseem
  //void save(FILE *fp) const;
 private:
  //void load(FILE *fp);
 public:
  //----------------

  void addControlPoint(float x, float y, float r);

  void removeDuplicateControlPoints();

  static void drawDisc(float x,float y, float z, float radius);

  // texture-mapping parameters
  bool useTexture() const { return _useTexture; }
  bool & useTexture() { return _useTexture; }
  float ufreq() const { return _ufreq; }
  float & ufreq() { return _ufreq; }
  float ustart() const { return _ustart; }
  float & ustart()  { return _ustart; }
//  float vfreq() const { return _vfreq; }
//  float & vfreq() { return _vfreq; }
//  float vstart() const { return _vstart; }
//  float & vstart() { return _vstart; }

  GLuint & lumAlphaNum() { return _lumAlphaNum; }
  GLuint lumAlphaNum() const { return _lumAlphaNum; }
  GLuint & alphaNum() { return _alphaNum; }
  GLuint alphaNum() const { return _alphaNum; }

  float arcLength();
};

struct _Vertex { 
  float x ; float y ; float z ;  
  inline _Vertex ( float _x , float _y, float _z ) 
              { y = _y ; x = _x ; z = _z;  } 
}; 



#endif

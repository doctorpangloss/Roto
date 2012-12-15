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



#include <stdio.h>
#include <stdlib.h>

#include "Stroke.h"
//#include "main.h"
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

GLUquadricObj * Stroke::qobj = NULL;

const double Stroke::MEAN_FILTER[3] = { 0.33, 0.34, 0.33 } ;  

Stroke::Stroke()
{
  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false ;

  _numLevels = 3;
   _cap = false;
  //  _cap = true ; 

  _extendLength = false;

  _radius = -1;
  _color = makeColor(0,0,0);
  _adjustedColor = makeColor(0,0,0);

  _texture_width = 1 ; 
  _texture_height= 1 ; 
  _texture_radius = 1.0; 
  

  if (qobj == NULL)
    qobj = gluNewQuadric();
    

  _useTexture = false;

  _lumAlphaNum = 0;
  _alphaNum = 0;
  _depth = 0;
  _tOffset=0;
}

/*
Stroke::Stroke(FILE* fp)
{
  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false; 

  _numLevels = 3;
  //  _cap = false;
  _cap = true ; 
  _extendLength = false;

  _texture_width = 1 ; 
  _texture_height= 1 ; 
  _texture_radius = 1.0; 

  _radius = -1;
  _color = makeColor(0,0,0);
  _adjustedColor = makeColor(0,0,0);

  if (qobj == NULL)
    qobj = gluNewQuadric();
    

  _useTexture = false;
  _lumAlphaNum = 0;
  _Alphanum = 0;
  _depth = 0;
 _tOffset=0;

  load(fp);
}
*/
Stroke::Stroke(float radius,GLuint color,StrokeNumT n)
{
  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false; 

  _numLevels = 3;
  
  _cap = true ;
  //  _cap = false;
  _extendLength = false;

  _radius = radius;
  _color = color;
  _adjustedColor = _color;
 _tOffset=0;

  _texture_width = 1 ; 
  _texture_height= 1 ; 
  _texture_radius = 1.0; 

  _num = n;

  if (qobj == NULL)
    {
      qobj = gluNewQuadric();
    }

  _useTexture = false;
  _lumAlphaNum = 0;
  _alphaNum = 0;
  _depth = 0;
}


Stroke::Stroke(int cx,int cy,float radius,GLubyte r,GLubyte g,GLubyte b,
	       StrokeNumT n)
{
  _control.push_back(PointR(cx,cy,radius));

  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false; 

  _numLevels = 3;
  //  _cap = false;
  _cap = true ; 

  _extendLength = false;
  
  _radius = radius;
  
  _color = makeColor(r,g,b);
  _adjustedColor = makeColor(r,g,b);
  _tOffset=0;
  _num = n;

  _texture_width = 1 ; 
  _texture_height= 1 ; 
  _texture_radius = 1.0; 

  if (qobj == NULL)
    {
      qobj = gluNewQuadric();
    }


  _useTexture = false;
  _lumAlphaNum = 0;
  _alphaNum = 0;
  _depth = 0;
}

Stroke::Stroke(int cx,int cy,float radius,GLuint color,StrokeNumT n)
{
  _control.push_back(PointR(cx,cy,radius));

  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false; 

  _numLevels = 3;
  //  _cap = false;
  _cap = true;
  _extendLength = false;
  
  _radius = radius;
  
  _texture_width = 1 ; 
  _texture_height= 1 ; 
  _texture_radius = 1.0; 

  _color = color;
  _adjustedColor = color;
  _tOffset=0;
  _num = n;
  
  if (qobj == NULL)
    {
      qobj = gluNewQuadric();
    }

  _useTexture = false;
    _lumAlphaNum = 0;
  _alphaNum = 0;
  _depth = 0;
}

Stroke::Stroke(const vector<PointR> & control, float radius,GLuint color,
	       StrokeNumT n) : _control(control)
{
  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false; 

  _numLevels = 3;

  _cap = true ; 
  //  _cap = false;
  _extendLength = false;
  
  _radius = radius;

  _texture_width = 1 ; 
  _texture_height= 1 ; 
  _texture_radius = 1.0; 
  
  _color = color;
  _adjustedColor = color;
  _tOffset=0;
  _num = n;

    if (qobj == NULL)
      {
	qobj = gluNewQuadric();
      }


     _useTexture = false;
    _lumAlphaNum = 0;
  _alphaNum = 0;
  _depth = 0;
}


Stroke::Stroke( const Stroke& s)
{
  _control = s._control;
  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false ; 

  _numLevels = s._numLevels;
  _cap = s._cap;
  _extendLength = s._extendLength;

  _texture_width = s._texture_width ;  
  _texture_height= s._texture_height ; 
  _texture_radius= s._texture_radius ; 

  _radius = s._radius;
  _color = s._color;
  _adjustedColor = s._adjustedColor;
  _tOffset=s._tOffset;

  _num = s._num;

  // SKQ 
  _useTexture = s._useTexture;
  _textureName = s._textureName ; 
  _texture_radius = s._texture_radius; 

  _ufreq = s._ufreq;
 // _vfreq = s._vfreq;
  _ustart = s._ustart;
//  _vstart = s._vstart;
  _depth = s._depth;

  _lumAlphaNum = s._lumAlphaNum;
  _alphaNum = s._alphaNum;

}


Stroke& Stroke::operator=( const Stroke& s)
{
  _control = s._control;
  _limit = _temp = NULL;
  _computed = false;
  _stripComputed = false; 

  _numLevels = s._numLevels;
  _cap = s._cap;
  _extendLength = s._extendLength;

  _texture_width = s._texture_width  ; 
  _texture_height= s._texture_height ; 
  _texture_radius= s._texture_radius ; 

  _radius = s._radius;
  _color = s._color;
  _adjustedColor = s._adjustedColor;
  _tOffset=s._tOffset;
  _num = s._num;


  // SKQ 
  _useTexture = s._useTexture;
  _textureName = s._textureName ; 
  _texture_radius = s._texture_radius; 

  _ufreq = s._ufreq;
 // _vfreq = s._vfreq;
  _ustart = s._ustart;
 // _vstart = s._vstart;
  _depth = s._depth;
    _lumAlphaNum = s._lumAlphaNum;
  _alphaNum = s._alphaNum;


  return *this;
}


void Stroke::clear()
{
  _computed = false;
  _stripComputed = false; 
  _control.erase(_control.begin(),_control.end());
  if (_limit != NULL)
    _limit->erase(_limit->begin(),_limit->end());
}

Stroke::~Stroke()
{
  if (_limit != NULL)
    delete _limit;

  if (_temp != NULL)
    delete _temp;
}

void Stroke::addControlPoint(float x,float y,float r)
{
  if (_control.empty() || _control.back().x != x || _control.back().y != y)
    {
      _control.push_back(PointR(x,y,r));
      _computed = false;
      _stripComputed = false; 
    }
}

void Stroke::drawWideLineCurve( ) 
{


  if (!_computed){ 
    computeLimitCurve();
  }

  float dx,dy,mag;
  PointR p0;
  PointR p1;
  PointR p2;

  if (_limit->empty()){ 
    return;
  }

  p0 = (*_limit)[0];

  if (_limit->size() == 1){

    if (_cap)
      drawDisc(p0.x,p0.y,depth(),p0.r);//radius());
    return;
  }

  
  p1 = (*_limit)[1];


  dx = p1.y - p0.y;
  dy = p0.x - p1.x;
  
  mag = sqrt(dx*dx + dy*dy);
  if (mag != 0){
    dx /= mag;
    dy /= mag; 
  }
  

  glBegin(GL_TRIANGLE_STRIP);
 
  const float z = depth();

  glVertex3f(p0.x + MAX(0,p0.r+_tOffset) * dx, p0.y + MAX(0,p0.r+_tOffset) * dy,z);
  glVertex3f(p0.x - MAX(0,p0.r+_tOffset) * dx, p0.y - MAX(0,p0.r+_tOffset) * dy,z);

  int curve_length  = _limit->size() -1 ; 
  for(int i=1;i< curve_length -1;i++)
    {

      p0 = (*_limit)[i-1];
      p2 = (*_limit)[i+1];

      dx = p2.y - p0.y;
      dy = p0.x - p2.x;


      mag = sqrt(dx*dx + dy*dy);

      if (mag != 0){
	dx /= mag;
	dy /= mag;
      }

      glVertex3f(p0.x + MAX(0,p0.r+_tOffset) * dx, p0.y + MAX(0,p0.r+_tOffset) * dy,z);
      glVertex3f(p0.x - MAX(0,p0.r+_tOffset) * dx, p0.y - MAX(0,p0.r+_tOffset) * dy,z); 

    }

  p0 = (*_limit)[_limit->size()-2];
  p1 = (*_limit)[_limit->size()-1];

  dx = p1.y - p0.y;
  dy = p0.x - p1.x;
    
  mag = sqrt(dx*dx + dy*dy);

  if (mag != 0) {
    dx /= mag;
    dy /= mag;
  }

  glVertex3f(p1.x + MAX(0,p1.r+_tOffset) * dx, p1.y + MAX(0,p1.r+_tOffset) * dy,z);

  glVertex3f(p1.x - MAX(0,p1.r+_tOffset) * dx, p1.y - MAX(0,p1.r+_tOffset) * dy,z);

  glEnd();

}

void Stroke::drawLines(vector<PointR> *curve)
{
  
  glBegin(GL_LINE_STRIP);

  for(unsigned int i=0;i<curve->size();i++)
    {
      PointR & p = (*curve)[i];
      glVertex2f(p.x,p.y);
    }
  glEnd();
}

void Stroke::drawDisc(float x,float y,float z, float radius)
{
  glPushMatrix();
  glTranslatef(x,y,z);
  gluDisk(qobj,0,radius,NUM_SLICES,1);
  glPopMatrix();
}


/******************************
 * 
 * 
 * inputs: 
 *            orientation  - 0 to draw Start Cap, 1 to draw End Cap
 *************************************************************/

void Stroke::computeCapPoints( PointR p0, float dx, float dy, float * cap_offsets, bool orientation ) { 
    

    float dx1, dy1; 
    float radiusI = RAMP_FACTOR*MAX( 0.0, p0.r + _tOffset); 
    float radiusO = RAMP_FACTOR_TWO*MAX( 0.0, p0.r + _tOffset); 
    float capX_inner, capY_inner, capX_outer, capY_outer; 
    float theta; 
    float texturePos; 

    if ( orientation == 0 ) { 
      inner_cap1.clear() ; 
      outer_cap1.clear() ; 
      startCapTexture_inner.clear() ;
      startCapTexture_outer.clear() ;

      texturePos = _ufreq*radiusO; 
      
    }else{
      inner_cap2.clear() ; 
      outer_cap2.clear() ; 
      endCapTexture_inner.clear() ; 
      endCapTexture_outer.clear() ; 
      
      vector<float>::iterator tmp = textureU.end() ; 
      tmp--; 
      texturePos = *tmp;

    } 

    // calculte THETA &  avoid divides by zero 
    if ( dx < 0.000001 && dx > -0.000001 && dy < 0.000001 && dy > -0.000001 ) { 
      
      theta = 0.0 ; 

      float signedDY = dy;
      if ( orientation == 0 ) { 
	signedDY*=-1; 
      } 

      if ( signedDY >= 0 ) { 
	theta+= M_PI; 
      } 

    } else if ( dx < 0.000001 && dx > -0.000001) { 

      // tan undefined
      theta = M_PI / 2.0 ;
      
      // correct the sign      
      float signedDY = dy;
      if ( orientation == 0 ) { 
	signedDY*=-1; 
      } 

      if ( signedDY < 0 ) { 
	theta*=-1; 
      } 
      
    } else if ( orientation == 0 ) { 
      theta = atan2( -dy , -dx ) ; 
    }else{
      theta = atan2( dy, dx ); 
    }

    // debug

    if ( ! finite(theta) ) { 
      printf("dx: %f, dy: %f\n", dx, dy); 
      printf("theta: %f\n", theta) ; 
      assert(finite(theta)); 
    } 

    float texture_theta;




    for ( int slice = 0 ; slice <= NUM_SLICES ; slice ++ ) { 

      dx1 = cos(theta + (float) slice * M_PI / NUM_SLICES);
      dy1 = sin(theta + (float) slice * M_PI / NUM_SLICES);
      
      texture_theta = (float)slice*M_PI/NUM_SLICES ; 

      // no texture cap
      if ( !_useTexture ) { 
	capX_inner = radiusI * dx1 ;
	capY_inner = radiusI * dy1 ; 

	capX_outer = radiusO * dx1 ;
	capY_outer = radiusO * dy1 ; 

      } else if  (_textureType == paper)  { 
	
	// paper texture cap 
	capX_inner =  cap_offsets[slice]*dx1 + radiusI * dx1 ; 
	capY_inner =  cap_offsets[slice]*dy1 + radiusI * dy1 ;

	capX_outer =  cap_offsets[slice]*dx1 + radiusO * dx1 ; 
	capY_outer =  cap_offsets[slice]*dy1 + radiusO * dy1 ;

      }else{
	// stroke texture cap 
	capX_inner = radiusI * dx1 ;
	capY_inner = radiusI * dy1 ; 

	capX_outer = radiusO * dx1 ;
	capY_outer = radiusO * dy1 ; 

	Point textureOuter; 
	Point textureInner; 

	if ( orientation == 0 ) { 

	  textureOuter.x = texturePos -_ufreq*radiusO*sin( texture_theta ) ; 
	  textureOuter.y = 0.5 - 0.5*_texture_radius*cos( texture_theta ) ; 
	  startCapTexture_outer.push_back( textureOuter ) ; 


	  textureInner.x = texturePos - _ufreq*radiusI*sin( texture_theta ) ; 
	  textureInner.y = 0.5 - 0.25*_texture_radius*cos( texture_theta ) ; 
	  startCapTexture_inner.push_back( textureInner ) ; 

	  
        }else{ 

	  textureOuter.x = texturePos -_ufreq*radiusO*sin( texture_theta + M_PI ) ; 
	  textureOuter.y = 0.5 - 0.5*_texture_radius*cos( texture_theta  + M_PI ) ; 
	  endCapTexture_outer.push_back( textureOuter ) ; 


	  textureInner.x = texturePos - _ufreq*radiusI*sin( texture_theta + M_PI) ; 
	  textureInner.y = 0.5 - 0.25*_texture_radius*cos( texture_theta  + M_PI) ; 
	  endCapTexture_inner.push_back( textureInner ) ; 


	}

      } 

      if ( orientation == 0 ) { 
	inner_cap1.push_back( Point( capX_inner, capY_inner ) );  
	outer_cap1.push_back( Point( capX_outer, capY_outer ) );  
      } else { 
	inner_cap2.push_back( Point( capX_inner, capY_inner ) );  
	outer_cap2.push_back( Point( capX_outer, capY_outer ) );  
      }

 
    } 

    if ( orientation == 0 ) {
      startCap.p0 = p0 ; 
      startCap.dx = dx ; 
      startCap.dy = dy ; 
      startCap.orientation = false ; 
    }else{ 
      endCap.p0 = p0; 
      endCap.dx = dx; 
      endCap.dy = dy; 
      endCap.orientation = true; 
    }
} 

void Stroke::setUniformRadius(const float t) {
  for (vector<PointR>::iterator c = _limit->begin(); c!=_limit->end(); ++c)
    c->r = t;
}

void Stroke::computeStripPoints( const vector<PointR> *curve) { 


  left_inner.clear() ; 
  left_outer.clear() ; 
  right_inner.clear() ; 
  right_outer.clear() ; 
  textureU.clear() ; 

  if (curve->empty())
  {
    _stripComputed = true ; 
    return;
    
  }

  // make sure all points are finite 
  //int length = curve->size() ; 
  for ( int i = 0 ; i < curve->size() ; i++ ) {     
    
    PointR p = (*curve)[i];
    
    if ( !finite(p.x) ) { 
      printf("curve point i x value = %f\n", p.x); 
      assert(finite(p.x)); 
    } 
    if ( !finite(p.y) ) { 
      printf("curve point i y value = %f\n", p.x); 
      assert(finite(p.y)); 
    } 
  } 
  

  // get a list of random offsets to use to jitter the line. 
  // MOVE THIS INTO THE IF STATEMENT... 
  float ratio = ((float) _texture_width) / (float) _texture_height ;
  ratio/=_texture_radius; 

  float* left_offsets_buffer = new float [curve->size()]; 
  float* right_offsets_buffer = new float [ curve->size()] ;

  float* left_offsets = new float [curve->size()]; 
  float* right_offsets = new float [ curve->size()] ;

  // random jitters for the caps
  // the first and last jitters of the caps semi circle
  // will be taken from the strip offsets so there is continuity
  // between the two pieces. 

  float* cap1_offsets = new float [ NUM_SLICES +1 ]; 
  float* cap2_offsets = new float [ NUM_SLICES +1 ]; 
  float* cap1_offsets_buffer = new float [ NUM_SLICES -1 ] ; 
  float* cap2_offsets_buffer = new float [ NUM_SLICES -1 ] ; 

  if (_useTexture && (_textureType == paper)){


    makeRandomArray( left_offsets_buffer, curve->size() ); 
    makeRandomArray( right_offsets_buffer, curve->size()); 

    blendArray( left_offsets, left_offsets_buffer, curve->size(),  curve->size(), 0 ) ; 
    blendArray( right_offsets, right_offsets_buffer, curve->size(), curve->size(), 0 ); 


    if ( _cap ) { 
      makeRandomArray( cap1_offsets_buffer, NUM_SLICES-1); 
      makeRandomArray( cap2_offsets_buffer, NUM_SLICES-1); 
      
      blendArray( cap1_offsets, cap1_offsets_buffer,  NUM_SLICES+1, NUM_SLICES-1, 1); 
      blendArray( cap2_offsets, cap2_offsets_buffer,  NUM_SLICES+1, NUM_SLICES-1, 1); 

      // now blend the cap offsets and triangle strip offests so the 
      // triangle strip and the caps will be continuous
      // ( obviously this is not a perfect way to blend them... but not big deal)


      left_offsets[0] = 0.5 * left_offsets[1] + 0.5 * cap1_offsets[0] ; 
      right_offsets[0] = 0.5 * right_offsets[1] + 0.5 * cap1_offsets[NUM_SLICES-2]; 

      left_offsets[curve->size()-1]  = 0.5*left_offsets [curve->size()-2] + 
	                               0.5*cap2_offsets [NUM_SLICES-2] ; 
      right_offsets[curve->size()-1] = 0.5*right_offsets[curve->size()-2] + 
                                       0.5*cap2_offsets [0]; 
     
      cap1_offsets[0] = left_offsets[0] ; 
      cap1_offsets[NUM_SLICES] =  right_offsets[0]; 
     
      cap2_offsets[0] = left_offsets[curve->size()-1]; 
      cap2_offsets[NUM_SLICES] =  right_offsets[curve->size()-1];  

    }
   
  }
  

  delete [] left_offsets_buffer ; 
  delete [] right_offsets_buffer ; 
  delete [] cap1_offsets_buffer ; 
  delete [] cap2_offsets_buffer ; 


  /* save all the points need to calculated the blended sides */ 

  assert(curve != NULL);

  float dx,dy,mag;
  float jitter_left, jitter_right; 
  PointR p0;
  PointR p1;
  PointR p2;

  // always keep track of the last set of dx,dy that !=0 
  // so we can tell the direction of the stroke even if there is
  // a pause in stroke movement
  float dirDX = 0.0; 
  float dirDY = -1.0;  

  if ( curve->size() == 0 ) { return ; } 

  p0 = (*curve)[0];

  _ufreq = 1.0 / ( 2.0*(float)p0.r*ratio ) ; 
  float cur_textureU = _ufreq * p0.r ; 
  textureU.push_back(cur_textureU) ; 
  
  if (curve->size() == 1)
  {
    if (_cap)
      drawDisc(p0.x,p0.y,depth(),p0.r);//radius());
    
    delete [] left_offsets ; 
    delete [] right_offsets ;
    delete [] cap1_offsets ; 
    delete [] cap2_offsets ; 

    _stripComputed = true; 
    return;
  }

  // curve->size is at least 2
  p1 = (*curve)[1];

  dx = p1.y - p0.y;
  dy = p0.x - p1.x;



  mag = sqrt(dx*dx + dy*dy);
 
  if (mag != 0){
    dx /= mag;
    dy /= mag; 
  }
  // last changed direction
  if ( (dx !=0.0 || dy != 0.0) && finite(dx) && finite(dy) ) { 
    
    dirDX = dx ; 
    dirDY = dy ; 
  } 


  // compute the start cap
  computeCapPoints( p0, dirDX, dirDY, cap1_offsets, 0 ); 
  startCap.textureU = cur_textureU ;  


  // compute the triangle strip
  float X2, Y2, X1, Y1, tmpX, tmpY ; 

  if ( _useTexture && _textureType == paper) {
      jitter_left = left_offsets[0] ; 
      jitter_right = right_offsets[0] ;
  }

  if ( _useTexture && _textureType == stroke ){
    float dist = sqrt((p1.x-p0.x)*(p1.x-p0.x)+(p1.y-p0.y)*(p1.y-p0.y));
    cur_textureU+=_ufreq * dist ;
    
    textureU.push_back(cur_textureU); 
  }

  if ( _useTexture && _textureType == paper ) { 
    X2 =  p0.x + jitter_right*dx + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
    Y2 =  p0.y + jitter_right*dy + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ; 
    
    tmpX = p0.x + jitter_right*dx + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ; 
    tmpY = p0.y + jitter_right*dy + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ; 

    right_outer.push_back(Point( tmpX , tmpY )); 

    
  }else{ 
    X2 =  p0.x + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
    Y2 =  p0.y + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ;
    
    tmpX = p0.x + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ; 
    tmpY = p0.y + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ; 

    right_outer.push_back(Point( tmpX , tmpY ) );
    
  }

  right_inner.push_back(Point( X2, Y2 )) ; 


  if ( _useTexture && _textureType == paper ) { 
    
    X1 =  p0.x - jitter_left*dx - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
    Y1 =  p0.y - jitter_left*dy - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ; 

    tmpX = p0.x - jitter_left*dx - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ; 
    tmpY = p0.y - jitter_left*dy - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ; 

    left_outer.push_back(Point( tmpX, tmpY )); 

		
  } else { 

    X1 =  p0.x - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
    Y1 =  p0.y - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ; 

    tmpX =  p0.x - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ; 
    tmpY =  p0.y - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ; 

    left_outer.push_back(Point( tmpX, tmpY ));		    
  } 

  left_inner.push_back(Point( X1, Y1 )); 
  
  // computing strip points
  for( unsigned int i=1;i<curve->size()-1;i++){
      p0 = (*curve)[i-1];
      p2 = (*curve)[i+1];

      float jitter_left, jitter_right ; 

      if ( _useTexture && _textureType == paper ) {
	  jitter_left  = left_offsets[i-1];
	  jitter_right = right_offsets[i-1];
      }
     
 
      dx = p2.y - p0.y;
      dy = p0.x - p2.x;

      mag = sqrt(dx*dx + dy*dy);

      if (mag != 0) {
	dx /= mag;
	dy /= mag;
      }
 
      if ( dx !=0.0 || dy != 0.0 ) { 
	dirDX = dx ; 
	dirDY = dy ; 
      } 
     

      // textureU
      if (_useTexture && _textureType == stroke ){
	p1 = (*curve)[i]; 
	float dist = sqrt((p1.x-p0.x)*(p1.x-p0.x)+(p1.y-p0.y)*(p1.y-p0.y));
	cur_textureU+=_ufreq * dist ; 

	textureU.push_back(cur_textureU); 
	
      } 


      if ( _useTexture && (_textureType == paper) ) { 
	
	  X2 =  p0.x + jitter_right*dx + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
	  Y2 =  p0.y + jitter_right*dy + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ; 
	  right_outer.push_back(Point( p0.x + jitter_right*dx + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ,  
				       p0.y + jitter_right*dy + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ) ); 
      } else { 

	X2 =  p0.x + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
	Y2 =  p0.y + RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ; 
	right_outer.push_back(Point( p0.x + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ,  
				     p0.y + RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ) ); 
      } 
      right_inner.push_back(Point( X2, Y2 ));

      if ( _useTexture && _textureType == paper ) { 

	X1 =  p0.x - jitter_left*dx - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
	Y1 =  p0.y - jitter_left*dy - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ; 
	left_outer.push_back(Point( p0.x - jitter_left*dx - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ,  
				    p0.y - jitter_left*dy - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ) ); 
      }else{ 

	X1 =  p0.x - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dx ; 
	Y1 =  p0.y - RAMP_FACTOR*MAX(0,p0.r+_tOffset) * dy ; 
	left_outer.push_back(Point( p0.x - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dx ,  
				    p0.y - RAMP_FACTOR_TWO*MAX(0,p0.r+_tOffset) * dy ) ); 
      }
      left_inner.push_back(Point( X1, Y1 )); // save this point
  }
  

  p0 = (*curve)[curve->size()-2];
  p1 = (*curve)[curve->size()-1];
  
  // textureU
  if (_useTexture && _textureType == stroke ){

    float dist = sqrt((p1.x-p0.x)*(p1.x-p0.x)+(p1.y-p0.y)*(p1.y-p0.y));
    cur_textureU+=_ufreq * dist ; 
    textureU.push_back(cur_textureU); 
  }


  if ( _useTexture && _textureType == paper ) { 
    jitter_right = right_offsets[curve->size()-1]; 
    jitter_left = left_offsets[curve->size()-1]; 
  }

  dx = p1.y - p0.y;
  dy = p0.x - p1.x;

  mag = sqrt(dx*dx + dy*dy);

  if (mag != 0){
      dx /= mag;
      dy /= mag;
  }

  if ( (dx !=0.0 || dy != 0.0) && finite(dx) && finite(dy) ) { 
    dirDX = dx ; 
    dirDY = dy ; 
  }
  

  if ( _useTexture && _textureType == paper ) {   
    X2 = p1.x + jitter_right*dx + RAMP_FACTOR*MAX(0,p1.r+_tOffset) * dx ; 
    Y2 = p1.y + jitter_right*dy + RAMP_FACTOR*MAX(0,p1.r+_tOffset) * dy ; 
    right_outer.push_back(Point( p1.x + jitter_right*dx + RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dx ,  
				 p1.y + jitter_right*dy + RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dy ) ); 
  }else{ 
    X2 = p1.x + RAMP_FACTOR*MAX(0,p1.r+_tOffset) * dx ; 
    Y2 = p1.y + RAMP_FACTOR*MAX(0,p1.r+_tOffset) * dy ; 
    right_outer.push_back(Point( p1.x + RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dx ,  
				 p1.y + RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dy ) ); 
  }

  right_inner.push_back( Point( X2, Y2));


  if ( _useTexture && _textureType == paper ) { 
    X1 = p1.x - jitter_left*dx - RAMP_FACTOR*MAX(0, p1.r+_tOffset) * dx ;  
    Y1 = p1.y - jitter_left*dy - RAMP_FACTOR*MAX(0, p1.r+_tOffset) * dy ; 
    left_outer.push_back(Point( p1.x - jitter_left*dx - RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dx ,  
				p1.y - jitter_left*dy - RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dy ) ); 
  }else{ 
    X1 = p1.x - RAMP_FACTOR*MAX(0, p1.r+_tOffset) * dx ;  
    Y1 = p1.y - RAMP_FACTOR*MAX(0, p1.r+_tOffset) * dy ; 
    left_outer.push_back(Point( p1.x - RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dx ,  
				p1.y - RAMP_FACTOR_TWO*MAX(0,p1.r+_tOffset) * dy ) );
  } 

  left_inner.push_back( Point(X1, Y1)); 

  // compute end cap points
  computeCapPoints( p1, dirDX, dirDY, cap2_offsets, 1 ); 
  endCap.textureU = cur_textureU ;  

  _stripComputed = true; 

}

void Stroke::drawTriangleStrip( Strip strip ) 
{

  vector<Point>::iterator left_itr; 
  vector<Point>::iterator right_itr; 
  vector<Point>::iterator end; 

  GLboolean blended = glIsEnabled(GL_BLEND);

  float left_alpha ; 
  float right_alpha ; 

  const float z = depth();
  float curve_length ; 


  float textureX1;      // the left edge of the texture; 
  float textureX2;      // the right edge of the texture ; 
  vector<float>::iterator textureU_itr = textureU.begin() ; // for stroke texture


  float X1, Y1, X2, Y2; 

  if ( left_outer.size() <= 0 || left_inner.size() <= 0  || right_inner.size() <= 0 || right_outer.size() <= 0 ) { 
    printf("returning size==0\n"); 
   
   return  ; 
  } 

  float cur_color[4]; 
  glGetFloatv( GL_CURRENT_COLOR, cur_color) ; 
 
  switch( strip ) { 

  case __left : 

    left_itr = left_outer.begin() ; 
    right_itr = left_inner.begin() ; 
    end = left_outer.end() ; 

    curve_length = left_outer.size() ; 
    assert( left_outer.size() == left_inner.size() ); 

    left_alpha = 1.0 ;  // G!
    right_alpha = 1.0 ; 
        
    textureX2 = 0.5 - _texture_radius*0.5;
    textureX1 = 0.5 - _texture_radius*0.25;  

    break ; 

  case __middle :
    
    left_itr = left_inner.begin() ; 
    right_itr = right_inner.begin() ; 
    end = left_inner.end() ; 

    curve_length = left_inner.size() ; 
    assert( left_inner.size() == right_inner.size() ); 

    left_alpha = 1.0 ; 
    right_alpha = 1.0 ;

    textureX2 = 0.5 - _texture_radius*0.25 ; 
    textureX1 = 0.5 + _texture_radius*0.25 ; 

    break ; 

  case __right : 

    left_itr = right_inner.begin() ; 
    right_itr = right_outer.begin() ; 
    end = right_inner.end() ; 

    curve_length = right_outer.size() ; 
    assert( right_outer.size() == right_inner.size() ); 

    left_alpha = 1.0 ; 
    right_alpha = 1.0 ; // G! 0


    textureX2 = 0.5 + _texture_radius*0.25 ; 
    textureX1 = 0.5 + _texture_radius*0.5  ;

    break; 

  }

  

  
  //glEnable( GL_BLEND ) ;  //G!
  
  glBegin(GL_TRIANGLE_STRIP);

    X2 = left_itr->x ; 
    Y2 = left_itr->y ; 
    
    left_itr++ ; 
    glColor4f( cur_color[0], cur_color[1], cur_color[2], left_alpha ) ; 
 
    if (_useTexture){

      if ( _textureType == paper ) { 
	glTexCoord2f( X2 / _texture_width , Y2 / _texture_height ) ;
      }else{
	if ( _textureType == stroke ) { 
	  glTexCoord2f( *textureU_itr , textureX2 ); 
	}
      } 
    }
    
 
    glVertex3f( X2, Y2, z); 
  

    X1 = right_itr->x ; 
    Y1 = right_itr->y ; 
    right_itr++ ; 
    
    glColor4f( cur_color[0], cur_color[1], cur_color[2], right_alpha ); 
    
    if (_useTexture){
      if ( _textureType == paper ) { 
	glTexCoord2f( X1 / _texture_width ,Y1/ _texture_height ) ;
      }else{ 
	if ( _textureType == stroke ) {  
	  glTexCoord2f( *textureU_itr , textureX1 ); 
	}
      }
    }    
    
    glVertex3f( X1, Y1, z); 
    
    if ( _textureType == stroke && _useTexture  ) {  
      textureU_itr++;  
    }

    while (  left_itr != end ) { 
      
      X2 = left_itr->x ; 
      Y2 = left_itr->y ; 
      
      left_itr++ ; 
      glColor4f( cur_color[0], cur_color[1], cur_color[2], left_alpha ); 
      if (_useTexture){
	if ( _textureType == paper ) { 
	  glTexCoord2f(X2 / _texture_width , Y2 /_texture_height ) ;
	}else{ 
	  if ( _textureType == stroke ) {  
	    glTexCoord2f( *textureU_itr , textureX2 ); 
	    // glTexCoord2f( textureX2, *textureU_itr ); 
	  }
	}
      }
      
      glVertex3f( X2, Y2, z); 
      
      X1 = right_itr->x ; 
      Y1 = right_itr->y ; 
      

      right_itr++ ; 
      glColor4f( cur_color[0], cur_color[1], cur_color[2], right_alpha ); 
      
      if (_useTexture){
	if ( _textureType == paper ) { 
	  glTexCoord2f( X1 / _texture_width , Y1 / _texture_height ) ;
	}else{ 
	  if ( _textureType == stroke ) {  

	    glTexCoord2f( *textureU_itr , textureX1 ); 
	  }
	}
      }  
      
      glVertex3f( X1, Y1, z); 
      if ( _useTexture &&  _textureType == stroke ) { 
	textureU_itr++ ;
      } 
      
    }
    
  glEnd();

  glColor4f( cur_color[0], cur_color[1], cur_color[2], cur_color[3]); 

  if ( !blended ) { 
    glDisable(GL_BLEND); 
  }

  assert(!glGetError()); 

} 

/* draws the stroke */ 
void Stroke::scanConvert( )
{

  GLboolean text2D_enabled = glIsEnabled(GL_TEXTURE_2D) ; 

  if (_useTexture){
    
    glBindTexture(GL_TEXTURE_2D, _textureName ); 
    glEnable( GL_TEXTURE_2D ); 

  }

  //  _cap= false; 
  _cap = true; 


  if (_cap){
    drawCap(startCap, inner_cap1, outer_cap1, 0 ); 
  }

   drawTriangleStrip( __middle ) ;

   drawTriangleStrip( __right ) ;

   drawTriangleStrip( __left ) ; 
   
  if (_cap){
    drawCap(endCap, inner_cap2, outer_cap2, 1 ); 
  }


  if ( !text2D_enabled ) { 
    glDisable(GL_TEXTURE_2D); 
  }

}

void Stroke::drawControlPolygon()
{
  drawLines(&_control);
}

void Stroke::drawLineCurve()
{
  if (!_computed){ 
    computeLimitCurve();
  }

  drawLines(_limit);
}

void Stroke::render()
{
  // SKQ
  float cur_color[4]; 

  glGetFloatv( GL_CURRENT_COLOR, cur_color);  

  if (!_computed){ 
    computeLimitCurve();
  }

  if ( ! _stripComputed ) { 

    computeStripPoints( _limit ) ; 
  }

  scanConvert( );
}

float Stroke::arcLength()
{
  if (!_computed)
	  computeLimitCurve();
  
  float length = 0;
 
  int limitLength  = _limit->size() -1 ; 

  for (int i=0 ; i < limitLength ; i++ ) 
  {
     
    PointR & p0 = (*_limit)[i];
    PointR & p1 = (*_limit)[i+1];
    
    float dist = sqrt((p1.x-p0.x)*(p1.x-p0.x)+(p1.y-p0.y)*(p1.y-p0.y));
    length += dist;
    
  }

  if (_cap && (limitLength >=2))
    {
      length += (*_limit)[0].r + (*_limit)[_limit->size()-1].r;

      //	  length += 2*_radius;
    }else{ 
      length = 0.0 ; 
    }
  
  return length;
}

void Stroke::subdivideCubicBSpline(vector<PointR> * inputCurve, 
				   vector<PointR> * outputCurve)
{
  outputCurve->erase(outputCurve->begin(),outputCurve->end());

  if (inputCurve->size() < 1)
    return;

  PointR pi0;
  PointR pi1;
  PointR pi2;

  pi0 = (*inputCurve)[0];

  outputCurve->push_back(pi0);

  if (inputCurve->size() == 1)
    return;

  if (inputCurve->size() == 2)
    {
      pi1 = (*inputCurve)[1];

      outputCurve->push_back(pi1);

      return;
    }

  pi1 = (*inputCurve)[1];

  outputCurve->push_back((pi0 + pi1)/2);

  for(unsigned int i=1;i<inputCurve->size()-1;i++)
    {
      pi0 = (*inputCurve)[i-1];
      pi1 = (*inputCurve)[i];
      pi2 = (*inputCurve)[i+1];

      outputCurve->push_back((pi0 + pi1*6 + pi2)/8);

      outputCurve->push_back( (pi1 + pi2)/2);
    }
	
  outputCurve->push_back(pi2);

			 
}

/*
void Stroke::subdivideFourPoint(vector<PointR> * inputCurve, 
				vector<PointR> * outputCurve)
{
  outputCurve->erase(outputCurve->begin(),outputCurve->end());

  if (inputCurve->size() < 1)
    return;

  PointR pi0;
  PointR pi1;
  PointR pi2;
  PointR pi3;

  if (inputCurve->size() == 1)
    {
      pi0 = (*inputCurve)[0];
      outputCurve->push_back(pi0);kPoint(pi0.x,pi0.y));

      return;
    }

  if (inputCurve->size() == 2)
    {
      pi0 = (*inputCurve)[0];
      pi1 = (*inputCurve)[1];
	
      outputCurve->push_back(Point(pi0.x,pi0.y));
      outputCurve->push_back(Point((pi0.x+pi1.x)/2,(pi0.y+pi1.y)/2));
      outputCurve->push_back(Point(pi1.x,pi1.y));
	
      return;
    }

  pi0 = (*inputCurve)[0];
  pi1 = (*inputCurve)[1];

  Point piminus1(2*pi0.x - pi1.x,2*pi0.y - pi1.y);

  pi0 = (*inputCurve)[inputCurve->size()-1];
  pi1 = (*inputCurve)[inputCurve->size()-2];

  Point piplus1(2*pi0.x - pi1.x,2*pi0.y - pi1.y);
    
  for(int i=0;i<inputCurve->size()-1;i++)
    {
      pi0 = (i==0 ? piminus1 : (*inputCurve)[i-1]);
      pi1 = (*inputCurve)[i];
      pi2 = (*inputCurve)[i+1];
      pi3 = (i==inputCurve->size()-2? piplus1:(*inputCurve)[i+2]);

      outputCurve->push_back(Point( pi1.x, pi1.y));

      outputCurve->push_back(Point( (-pi0.x + 9*pi1.x + 9*pi2.x - pi3.x)/16,
				    (-pi0.y + 9*pi1.y + 9*pi2.y - pi3.y)/16));
    }

  pi0 = inputCurve->back();

  outputCurve->push_back(Point(pi0.x,pi0.y));
}
*/

void Stroke::subdivide(vector<PointR> * inputCurve, vector<PointR> * outputCurve)
{
  subdivideCubicBSpline(inputCurve,outputCurve);

  /*
    switch(curveType)
    {
    case CUBIC_BSPLINE:
    subdivideCubicBSpline(inputCurve,outputCurve);
    break;

    case FOUR_POINT:
    subdivideFourPoint(inputCurve,outputCurve);
    break;

    default:
    fprintf(stderr, "Illegal subdivision scheme selected\n");
    exit(-1);
    }
  */
}

void Stroke::computeLimitCurve()
{

  if (_limit == NULL)
    _limit = new vector<PointR>();

  if (_temp == NULL)
    _temp = new vector<PointR>();

  subdivide(&_control,_limit);

  if (_extendLength && _control.size() > 1)
  {
	PointR & p0 = _control[0];
	PointR & p1 = _control[1];

	PointR v = p0 - p1;
	v.normalize();
	v = v*_radius;
	_limit->insert(_limit->begin(),p0 + v);


	PointR &pn_1 = _control[_control.size()-2];
	PointR &pn =   _control[_control.size()-1];

	v = pn - pn_1;
	v.normalize();
	v = v*_radius;
	_limit->push_back(pn + v);
  }

  for(int i=0;i<_numLevels/2;i++)
    {
      subdivide(_limit,_temp);
      subdivide(_temp,_limit);
    }

   _computed = true;
}
  
  

void Stroke::print(FILE * fp)
{
  fprintf(fp,"Stroke %d: Color = %X. Radius %f. ",
	 (long int)_num,color(),radius());

  fprintf(fp,"Points = ");
  for(vector<PointR>::iterator it = _control.begin(); it!=_control.end(); ++it)
    fprintf(fp,"[%f,%f] ",(*it).x,(*it).y);

  fprintf(fp,"\n\n");
}

/*
void Stroke::save(FILE *fp) const {
  int dummy = _control.size();
  fwrite(&dummy, sizeof(int), 1, fp); // numPoints
  for(vector<PointR>::const_iterator it = _control.begin(); it!=_control.end(); ++it)
    fwrite(it, sizeof(PointR), 1, fp);
}

void Stroke::load(FILE *fp) {
  int dummy, res;
  _numLevels=0;
  res = fread(&dummy, sizeof(int), 1, fp); // numPoints
  assert(res==1);
  _control.reserve(dummy);
  PointR ptr;
  for (int i=0; i<dummy; i++) {
    res = fread(&ptr, sizeof(PointR), 1, fp);
    assert(res==1);
    _control.push_back(ptr);
  }
}
*/


void Stroke::removeDuplicateControlPoints()
{
  //int startsize = _control.size();

  vector<PointR>::iterator it = _control.begin();
  
  while (it != _control.end())
    {
      vector<PointR>::iterator next = it;   ++ next;
      
      if (next == _control.end())
	return;
      
      if ((*it) == (*next))
	_control.erase(it);
      
      it = next;
    }
}

/* makes an array of random numbers between 0 and 1 */
void Stroke::makeRandomArray( float arrOut[] , int size ) {

  float offset, rand_sign ; 
  for ( int i = 0 ; i < size ; i++) { 
    
    offset = 1.0* rand() / RAND_MAX ; 
    rand_sign = ( 1.5*rand()/ RAND_MAX ); 
   
    if ( rand_sign > 0.5 ){
      offset = -1.0*offset ; 
    }
    arrOut[i] = offset; 
  }
}

void Stroke::copyArray( float arrOut[] , float arrIn[], int size){ 
  
  for ( int i = 0 ; i < size ; i++) { 
    arrOut[i] = arrIn[i]; 
  } 
} 

void Stroke::blendArray( float arrOut[], float arrIn[],  int sizeOut, int sizeIn,  int startIndex  ) { 
  
  int mid = 1 ; 

  assert( startIndex >= 0 && startIndex < sizeOut  ) ; 
  assert( sizeOut >= startIndex + sizeIn ) ; 
  
  for ( int j = 0 ; j < startIndex ; j++ ) {  
    arrOut[j] = 0.0; 
  } 

  // don't blend the first point
  if ( sizeIn > 0 ){ 
    arrOut[startIndex] = arrIn[0] ; 
  }

  int outI = startIndex + 1 ; 
  for ( int i = 1 ; i < sizeIn-1 ; i++ ){ 
    
    arrOut[outI] = 0.0;  
    for ( int neigh = -1; neigh <= 1; neigh ++) {
      arrOut[outI] += MEAN_FILTER[mid+neigh]*arrIn[i+neigh];
    }
    outI ++ ; 
  } 
  
  // don't blend last point
  arrOut[outI] = arrIn[sizeIn-1] ;
  outI ++ ; 

  while( outI < sizeOut ) { 
    arrOut[outI] = 0.0;
    outI++; 
  } 

}  



void Stroke::drawCap(const _CapData & cap, vector<Point> &  inner_points, vector<Point> &  outer_points, bool orientation ) { 

  PointR p0 = cap.p0; 

  //float textureR ; 
  float cur_color[4] ; 
  glGetFloatv( GL_CURRENT_COLOR, cur_color);

  glColor4f( cur_color[0], cur_color[1], cur_color[2], 1.0 ) ;
  GLboolean blended ; 
  glGetBooleanv( GL_BLEND, &blended ) ; 

  float texturePos = -1.0 ; 
  vector<float>::iterator textureU_itr;  

  if ( _textureType == stroke && textureU.size() > 0 ) { 
    if ( orientation == 0 ) { 
      // drawing the start cap 
      texturePos = _ufreq*MAX(0.0, (cap.p0).r); 
    }else{ 
      // drawing the end cap 
      textureU_itr = textureU.end() ;
      textureU_itr--; 
      texturePos = *textureU_itr; 
    } 
  }

  glPushMatrix() ; 
  glTranslatef(p0.x, p0.y, _depth); 
  


  vector<Point>::iterator inner = inner_points.begin();
  vector<Point>::iterator outer ; 
  vector<Point>::iterator inner_end = inner_points.end() ; 
  vector<Point>::iterator outer_end = outer_points.end() ; 
  vector<Point>::iterator cap_inner , cap_outer; 
  
  if ( orientation == 0 ) { 
    cap_inner = startCapTexture_inner.begin() ; 
    cap_outer = startCapTexture_outer.begin() ; 
  } else { 
    cap_inner = endCapTexture_inner.begin() ; 
    cap_outer = endCapTexture_outer.begin() ; 
  } 

  glBegin( GL_TRIANGLE_FAN );

    if (_useTexture){
      
      if ( _textureType == paper ) { 
	glTexCoord2f(p0.x / _texture_width , p0.y / _texture_height); 
	
      }else if ( _textureType == stroke ) { 
		
	// stroke texture
	glTexCoord2f( texturePos , 0.5 ); 

      }
    }
    
    
    glVertex3i(0,0,0);
    

    while ( inner!= inner_end ){ 
      
      if (_useTexture){
	
	if ( _textureType == paper ) { 
	  glTexCoord2f( (p0.x + inner->x) / _texture_width , (p0.y + inner->y) / _texture_height ) ;
	}else{
	  if ( _textureType == stroke ) { 
	    
	    // stroke texture
	    //glTexCoord2f(textureU + _ufreq * p0.r * (-sin(i*M_PI/NUM_SLICES)+1)/2, (-cos(i*M_PI/NUM_SLICES)+1)/2 )
	    //	    glTexCoord2f(cap_inner->x, cap_inner->y); 
	    glTexCoord2f(cap_inner->x, cap_inner->y); 
	    cap_inner++; 
	  }
	} 
      }
      
      glVertex3f( inner->x, inner->y, 0); 
      inner++ ;

    } 
  glEnd() ; 

  
  if ( orientation == 0 ) { 
    cap_inner = startCapTexture_inner.begin() ; 
  } else { 
    cap_inner = endCapTexture_inner.begin() ; 
  } 

  inner = inner_points.begin() ;
  outer = outer_points.begin() ; 

  glBegin( GL_TRIANGLE_STRIP ); 

    while ( outer!= outer_end && inner != inner_end ){ 
 

      // inner point
      glColor4f( cur_color[0], cur_color[1], cur_color[2], 1.0 ) ; 
      if ( _useTexture  ) { 
	
	if ( _textureType == paper ) { 
	  glTexCoord2f( ( p0.x + inner->x)/_texture_width, (p0.y + inner->y)/_texture_height ) ; 
	} else if ( _textureType == stroke ) { 

	  glTexCoord2f( cap_inner->x, cap_inner->y ) ; 
	  cap_inner++; 
	} 
      }
      glVertex3f( inner->x, inner->y, 0) ; 
      
      glColor4f( cur_color[0], cur_color[1], cur_color[2], 0.0 ) ; 
      // outer point
      if ( _useTexture ) { 
	if ( _textureType == paper ) { 
	  glTexCoord2f( ( p0.x + outer->x)/_texture_width, (p0.y + outer->y)/_texture_height ) ; 
	} else if ( _textureType == stroke ) { 

	  glTexCoord2f( cap_outer->x, cap_outer->y ) ; 
	  cap_outer++; 
	}
      }
      
      glVertex3f( outer->x, outer->y, 0); 
      
      
      outer++ ; 
      inner++ ; 
    } 
  glEnd() ; 


  /**********

  glColor4f(1.0, 0.0, 0.0, 1.0) ; 
  glPointSize(4.0);
  
  // debug 
  inner = inner_points.begin();
  outer = outer_points.begin();
  glBegin(GL_POINTS) ; 
  while ( inner != inner_end ) { 
    glVertex3f( inner->x, inner->y, 0);
    glVertex3f( outer->x, outer->y, 0); 
      
    inner++ ; 
    outer++ ; 
   
  } 
  
  glEnd() ; 

  ***********/



  glPopMatrix() ; 
  glColor4f( cur_color[0], cur_color[1], cur_color[2], cur_color[3] ) ; 

} 

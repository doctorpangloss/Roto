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



#ifndef FITCURVES_H
#define FITCURVES_H


#include "GraphicsGems.h"					
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef Point2 *BezierCurve;

/* Forward declarations */
void		FitCurve(Point2 *d, int nPts, double error);
double aseemFitCurve(Point2* d,int nPts,Point2* res);
static	void		FitCubic();
static	double		*Reparameterize();
static	double		NewtonRaphsonRootFind();
static	Point2		BezierII();
static	double 		B0(), B1(), B2(), B3();
static	Vector2		ComputeLeftTangent();
static	Vector2		ComputeRightTangent();
static	Vector2		ComputeCenterTangent();
static	double		ComputeMaxError();
static	double		*ChordLengthParameterize();
static	BezierCurve	GenerateBezier();
static	Vector2		V2AddII();
static	Vector2		V2ScaleIII();
static	Vector2		V2SubII();




#endif

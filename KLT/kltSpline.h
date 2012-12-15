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



// to be included in klt.h


public:

void setupSplineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, MultiSplineData* mt, bool redo);

private:

void splineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, bool redo);

void safeSplineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, const int level);

bool innerSplineTrack(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, const int level);

double createSplineMatrices(const KLT_FullCPyramid** pyrms, const KLT_FullPyramid** pyrmsE, ZVec* Z, const int level, CSplineKeeper* keep/*,
			    CSplineKeeper *imgKeep, CSplineKeeper *edgeKeep*/); 

void initSplineEdgeMins(const KLT_FullPyramid** pyrmsE, const int level);

void putGradOnJ(Ubcv& J1, Ubcv& J2, Ubcv& J3, const double grad[24], const int vars[4]);
void putVarsOnJ(Ubcv& J1, Ubcv& J2, const TrackSample& ds, const int vars[4], const double coef);
void putVarsOnJ(Ubcv& J, const TrackSample& ds, const int vars[4], const double c1, const double c2);

// Given a TrackSample with loc & normal, calculate location j off
void splineLoc(const TrackSample& ds, int j, float invsubs, Vec2f* loc);

void noSubSplineLoc(const TrackSample& ds, int j, float subs, Vec2f* loc);

// Calculate 2x8 derivative K matrix for sample, j off tangent
void calculateSplineK(const TrackSample& s, double k, int j, double invsubs, double K[16]);

void doSplineShapeInterp();

void doOneStep(CSplineKeeper* keep, double* x);

bool _useD1;
bool _useD2;
bool _useD0;

QTime time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11;
int tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9, tt10, tt11;
FILE* fpo;

public:
MultiSplineData* _mts;




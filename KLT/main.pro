TEMPLATE = lib

LIBS   += 
INCLUDEPATH += /homes/gws/aseem/boost
CONFIG += qt thread staticlib X11
QMAKE_CXXFLAGS_RELEASE += -DNDEBUG -ffast-math

HEADERS = base.h kernels.h klt_util.h pyramid.h error.h klt.h error.h myMontage.h bezSpline.h contCorr.h\
linearSolver.h myIMatrix.h mysparsemat.h keeper.h multiTrackData.h multiKeeper.h debugJac.h \
HB_Sweep.h HB_OneCurve.h cubicEval.h  contCorr.h  genKeeper.h\
multiSplineData.h splineKeeper.h kltSpline.h contEdgeMin.h buildingSplineKeeper.h \
diagMatrix.h multiDiagMatrix.h obsCache.h samples.h

SOURCES = klt.C error.c kernels.C klt_util.C pyramid.C linearSolver.C \
mysparsemat.C keeper.C multiKeeper.C HB_Sweep.C HB_OneCurve.C bezSpline.C \
multiSplineData.C splineKeeper.C kltSpline.C contCorr.C buildingSplineKeeper.C \
diagMatrix.C multiDiagMatrix.C obsCache.C

TARGET = KLT

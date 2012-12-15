TEMPLATE	= app
CONFIG		+= qt thread opengl x11 
INCLUDEPATH 	+= myCurveFit /homes/gws/aseem/boost
LIBS		+= KLT/libKLT.a  
QMAKE_CXXFLAGS_RELEASE += -DNDEBUG

HEADERS		= mymainwindow.h centralWidget.h imgSequence.h \
rotoscope.h dynArray.h jl_vectors.h llist.h geigerCorresponder.h\
imageGraph.h liveWire.h heap.h rotoCurves.h corresponder.h \
rotoCorresponder.h draw.h drawPath.h drawCurves.h rotoPath.h \
abstractPath.h interModule.h myCurveFit/FitCurves.h \
bBoxf2D.h KLT/klt.h rangeDialog.h drawCorresponder.h gcnode.h \
pixmaps.h joint.h Texture.h hertzStrokes/Stroke.h hertzStrokes/bitmap.h trackGraph.h \
KLT/contCorr.h KLT/bezSpline.h KLT/talkFitCurve.h drawContCorr.h transformer.h \
rotoRegion.h

SOURCES		= main.C mymainwindow.C centralWidget.C rotoscope.C \
imageGraph.C liveWire.C rotoCurves.C imgSequence.C corresponder.C \
draw.C drawPath.C drawCurves.C rotoPath.C hertzStrokes/Stroke.cpp \
geigerCorresponder.C abstractPath.C interModule.C \
myCurveFit/FitCurves.c myCurveFit/GGVecLib.c rangeDialog.C \
drawCorresponder.C rotoCorresponder.C Texture.C hertzStrokes/bitmap.cpp \
trackGraph.C drawContCorr.C transformer.C rotoRegion.C

TARGET		= main

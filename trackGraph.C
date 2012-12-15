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



#include "trackGraph.h"

bool TrackGraph::pathAlreadyThere(const RotoPath* rp) {
  PathV::const_iterator ok = find(_paths.begin(), _paths.end(), rp);
  return (ok != _paths.end());  
}

void TrackGraph::addPath(RotoPath* rp) {
  _paths.push_back(rp);
  processJointRecords(rp);
}

void TrackGraph::processConstraint(const RotoPath* constrained, const RotoPath* constrainee, 
		       const int edside, const int eeside) {
  //assert(_paths.back() == constrained); // should have JUST added constrained path

  PathV::const_iterator  pi = find(_paths.begin(), _paths.end(),constrained);
  JointRecord& jr = _joints[pi - _paths.begin()];
  jr.fixer[edside]=constrainee;
  jr.fixerSide[edside] = eeside;
  //jr.fixedLoc[edside] = constrainee->getElement( eeside * (constrainee->getNumElements()-1));
}


// CORRECTNESS: should remove joints if not in first keyframe
void TrackGraph::processJointRecords(const RotoPath* rp) {
  Joints::const_iterator c;
  JointRecord newJR;

  for (c=rp->bJoints()->begin(); c!=rp->bJoints()->end(); ++c) {
    PathV::const_iterator pc = find(_paths.begin(), _paths.end(), c->_rp);
    if (pc != _paths.end()) {
      traceVarLoc(pc, &newJR, c->_side, 0);      
      break;
    }
    else {
      newJR.which[0] = -1; newJR.side[0] = -1;
    }
  }

  for (c=rp->eJoints()->begin(); c!=rp->eJoints()->end(); ++c) {
    PathV::const_iterator pc = find(_paths.begin(), _paths.end(), c->_rp);
    if (pc != _paths.end()) {
      traceVarLoc(pc, &newJR, c->_side, 1);      
      break;
    }
    else {
      newJR.which[1] = -1; newJR.side[1] = -1;
    }
  }

  _joints.push_back(newJR);
}


void TrackGraph::traceVarLoc(PathV::const_iterator pc, JointRecord* jr, int pcside, const int side) {
  int index = pc - _paths.begin(), nindex, npcside;
  if (_joints[index].which[pcside]!=-1) {
    assert(_joints[index].side[pcside]!=-1);
    nindex = _joints[index].which[pcside];
    npcside = _joints[index].side[pcside];
  }
  else {
    nindex = index; npcside = pcside;
  }
  assert(_joints[nindex].which[npcside]==-1);
  assert(_joints[nindex].side[npcside]==-1);

  jr->which[side] = nindex;
  jr->side[side] = npcside;  
}


void TrackGraph::printStructure() {
  printf("%d curves\n",_paths.size());

  for (int i=0; i<_joints.size(); i++) {
    printf("%d: %d %d, %d %d\n",i,_joints[i].which[0], _joints[i].side[0],
	   _joints[i].which[1], _joints[i].side[1]);
  }

}


void TrackGraph::getKey0Paths(PathV* putHere, const int numFrames) {
  assert(putHere->size()==0);
  assert(_paths.size()>0);
  *putHere = _paths;
  uint i,j;
  for (i=0; i<numFrames; ++i)  // loop over time
    for (j=0; j < putHere->size(); ++j) { // loop over curves

      assert(((*putHere)[j])->prevC()); 
      (*putHere)[j] = ((*putHere)[j])->prevC();

      // bring constraining curves back to key0, too
      JointRecord& jr = _joints[j];
      if (jr.fixer[0] != NULL)
	jr.fixer[0] = jr.fixer[0]->prevC();
      if (jr.fixer[1] != NULL)
	jr.fixer[1] = jr.fixer[1]->prevC();
    }      
}

/*
void TrackGraph::startMulti(MultiTrackData* mt) {
  if (_key0paths.size()==0)
    getKey0Paths(&_key0paths, mt->_numFrames);
  assert(_key0paths.size() == _joints.size());
  mt->_nCurves = _key0paths.size();
  mt->_trackWidths = new Vec2i[mt->_nCurves];
  mt-> _crs = new CurveRecord[mt->_nCurves];
  

  int tot=0, i, j; // tot is number of curve points
  for (i=0; i<_key0paths.size(); ++i) {
    tot += _key0paths[i]->getNumElements();
    mt->_crs[i]._nPoints = _key0paths[i]->getNumElements();
    mt->_trackWidths[i].Set(_key0paths[i]->lowHeight(), _key0paths[i]->highHeight());
    mt->_crs[i]._useEdges = _key0paths[i]->trackEdges();
  }
  
  
  int cv=0, ccv;   // variable counters
  int totCurveVar=0;
  int which;
  for (i=0; i<_key0paths.size(); ++i) {
    which = _joints[i].which[0];
    if (which ==-1) {
      mt->_crs[i]._vars[0][0]=i;
      mt->_crs[i]._vars[0][1]=0;
      mt->_crs[i]._vars[0][2]=cv;
      ccv=1;
    }
    else {
      ccv=0;
      assert(which < i);
      if (_joints[i].side[0]==0) {
	for (j=0; j<3; j++) {
	  mt->_crs[i]._vars[0][j] = mt->_crs[which]._vars[0][j];
	  assert(mt->_crs[i]._vars[0][j] != -1);
	}
      }

      else {

	for (j=0; j<3; j++) {
	  mt->_crs[i]._vars[0][j] = mt->_crs[which]._vars[3][j];
	  assert(mt->_crs[i]._vars[0][j] != -1);
	}
      }
    }

    mt->_crs[i]._vars[2][0] = mt->_crs[i]._vars[1][0] = i;
    mt->_crs[i]._vars[1][1]=ccv;
    mt->_crs[i]._vars[1][2]=cv + ccv;
    ccv += _key0paths[i]->getNumElements()-3;
    mt->_crs[i]._vars[2][1]=ccv;
    mt->_crs[i]._vars[2][2]=cv + ccv++;

    which = _joints[i].which[1];
    if (which==-1) {
      mt->_crs[i]._vars[3][0]=i;
      mt->_crs[i]._vars[3][1]=ccv;
      mt->_crs[i]._vars[3][2]=cv + ccv++;
    }
    else {
      assert(which < i);
      if (_joints[i].side[1]==0) {
	for (j=0; j<3; j++) {
	  mt->_crs[i]._vars[3][j] = mt->_crs[which]._vars[0][j];
	  assert(mt->_crs[i]._vars[3][j] != -1);
	}
      }
      else {
	for (j=0; j<3; j++) {
	  mt->_crs[i]._vars[3][j] = mt->_crs[which]._vars[3][j];
	  assert(mt->_crs[i]._vars[3][j] != -1);
	}
      }
    }

    mt->_crs[i]._nVars = ccv;
    if (i==0) mt->_crs[0]._packedStart = 0;
    else mt->_crs[i]._packedStart = 
	   mt->_crs[i-1]._packedStart + mt->_crs[i-1]._nVars*(mt->_numFrames+1);
    cv += ccv*(mt->_numFrames-1);
    totCurveVar += ccv;
  }

  
  //assert(ccp==tot);
  assert(cv>0);
  printf("%d curve control points (one frame), %d variables (all frames)\n", tot, cv );
  mt->_nVars = totCurveVar;
}

// copy data from rotoPaths to Z in mt
void TrackGraph::finishStartingMulti(MultiTrackData* mt) {

  mt->_Z = new Vec2f[(mt->_numFrames+1) * mt->_nVars];

  int i,j, offset=0, t, numEm1;
  PathV currPaths = _key0paths;
  RotoPath *prev, *curr;

  for (i=0; i<mt->_nCurves; ++i) {  // loop over curves
    curr = _key0paths[i];
    JointRecord& jr = _joints[i];
    numEm1 = curr->getNumElements()-1;
    offset=0;
    for (t=0; t<=mt->_numFrames; ++t, offset+=mt->_crs[i]._nVars) { // loop over time
      assert(curr);

      // 0th point
      if (t>0 && t<mt->_numFrames && (jr.fixer[0] != NULL || curr->endFixed(0))) { // constraining curve, or nudged
	FixedCurvePoint fcp;
	fcp._frame = t-1;  fcp._i = 0;
	mt->_crs[i]._fixedLocs.push_back(fcp);
      }
      mt->_Z[mt->packed_index(i,0,t)] = 
	(t==0) ? curr->getElement(0) : prev->getNextLoc(0);
      
      // middle points
      for (j=1; j <numEm1; ++j) {
	mt->_Z[mt->packed_index(i,j) + offset] = 
	  (t==0) ? curr->getElement(j) : prev->getNextLoc(j);
      } // end loop inner curve points
      
      // last point
      if (t>0 && t<mt->_numFrames && (jr.fixer[1] != NULL || curr->endFixed(1))) { // constraining curve, or nudged
	FixedCurvePoint fcp;
	fcp._frame = t-1;  fcp._i = numEm1;
	mt->_crs[i]._fixedLocs.push_back(fcp);
      }
      mt->_Z[mt->packed_index(i,numEm1,t)] = 
	(t==0) ? curr->getElement(numEm1) : prev->getNextLoc(numEm1);
      

      prev = curr;
      curr = curr->nextC(); 
      
      if (jr.fixer[0] != NULL)
	jr.fixer[0] = jr.fixer[0]->nextC();
      if (jr.fixer[1] != NULL)
	jr.fixer[1] = jr.fixer[1]->nextC();
      
    } // end loop over time
  } // end loop over curves
  //printf("u\n");

}
*/
/*
// copy locations from Z in mt to actual rotoPaths
void TrackGraph::copyLocs(MultiTrackData* mt) {
  int i,j, offset=0, t, numEm1;
  PathV currPaths = _key0paths;
  RotoPath *curr;

  for (i=0; i<mt->_nCurves; ++i) {
    curr = _key0paths[i]->nextC();
    numEm1 = _key0paths[i]->getNumElements()-1;
    offset = mt->_crs[i]._nVars;
    
    for (t=1; t<mt->_numFrames; ++t, offset+=mt->_crs[i]._nVars) {
      assert(curr);

      curr->setElement(0, mt->_Z[mt->packed_index(i,0,t)]);

      for (j=1;  j<numEm1; ++j) {
	curr->setElement(j, mt->_Z[mt->packed_index(i,j) + offset]);
      }

      curr->setElement(numEm1, mt->_Z[mt->packed_index(i,numEm1,t)]);
      curr = curr->nextC();
    }
    //offset -= (mt->_numFrames-2) * mt->_crs[i]._nVars;
  }


}
*/

void TrackGraph::clear() {
  _paths.clear();
  _key0paths.clear();
  _joints.clear();
}

void TrackGraph::buildBackReconcileJoints(int bFrame, int aFrame) {
  RotoPath::buildBackReconcileJoints(_paths,bFrame,aFrame);
}


void TrackGraph::buildMulti(MultiSplineData* mts) {
  if (_key0paths.size()==0)
    getKey0Paths(&_key0paths, mts->_numFrames);
  assert(_key0paths.size() == _joints.size());
  mts->_nCurves = _key0paths.size();
  mts->_trackWidths = new Vec2i[mts->_nCurves]; 
  //mts->_edgemins = new ContEdgeMin[mts->_nCurves];
  mts->_edgeMins.reserve(mts->_nCurves);

  mts->_conts.reserve(mts->_numFrames*mts->_nCurves);
  mts->_splines.reserve((mts->_numFrames+1)*mts->_nCurves);
  mts->_numSegs.reserve((mts->_numFrames+1)*mts->_nCurves);

  int vc=0, var=0, i, j, numControls, which, n; 
  PathV currPaths = _key0paths;
  RotoPath *curr;
  for (i=0; i<mts->_nCurves; ++i) {
    curr = _key0paths[i];
    JointRecord& jr = _joints[i];

    mts->_trackWidths[i].Set(_key0paths[i]->lowHeight(), _key0paths[i]->highHeight());
    if (_key0paths[i]->trackEdges()) 
      mts->_useEdges.push_back(true);
    else 
      mts->_useEdges.push_back(false);
    mts->_edgeMins.push_back(vector<double>());
    //mts->_edgemins[i].setUseEdges(_key0paths[i]->trackEdges());

    for (j=0; j<=mts->_numFrames; ++j) {
      assert(curr);
      BezSpline* spline = new BezSpline(); 
      numControls = curr->getNumControls(); 
      assert(numControls>0);
      spline->setMapLength(numControls);
      assert(!spline->hasZMap());
      //spline->_ZMap = new int[numControls];
      spline->createZMap();
      if (j>0 && j<mts->_numFrames)
	spline->createVarMap();
	//spline->_varMap = new int[numControls];
      mts->_splines.push_back(spline);
      mts->_numSegs.push_back(curr->getNumSegs());
      if (j<mts->_numFrames) {
	assert(curr->getNextCont());
	mts->_conts.push_back(curr->getNextCont());  
	//printf("curve %d time %d %x identity: %d\n",i,j,curr->getNextCont(), curr->getNextCont()->identity());
      }
      printf("curve %d time %d curve %x, ends %d %d, fixers %x %x\n",i,j,curr,curr->endFixed(0), curr->endFixed(1),
	     jr.fixer[0], jr.fixer[1]);

      // First control
      if (j>0 && j<mts->_numFrames && (jr.fixer[0] != NULL || curr->endFixed(0))) { // constraining curve, or nudged
	FixedControl fcp;
	fcp._t = j;  fcp._n = 0; fcp._c = i;
	mts->_fixedControls.push_back(fcp);
      }

      which = _joints[i].which[0];
      if (which==-1) {// resident first point
	mts->_Z.push_back(*(curr->getCtrl(0)));
	spline->ZMap(0) = vc++;
	if (j>0 && j<mts->_numFrames)
	  spline->varMap(0) = var++;
      }
      else {
	spline->ZMap(0) = mts->end_tc(j,which,_joints[i].side[0]);
	if (j>0 && j<mts->_numFrames)
	  spline->varMap(0) = mts->end_tc_var(j,which,_joints[i].side[0]);
      }

      // Middle controls
      for (n=1; n<=numControls-2; ++n) {
	if (j>0 && j<mts->_numFrames && curr->internalFixed(n)) {
	  FixedControl fcp;
	  fcp._t = j;  fcp._n = n; fcp._c = i;
	  mts->_fixedControls.push_back(fcp);
	}
	mts->_Z.push_back(*(curr->getCtrl(n)));
	spline->ZMap(n) = vc++;
	if (j>0 && j<mts->_numFrames)
	  spline->varMap(n) = var++;
      }
      
      // Last Control
      if (j>0 && j<mts->_numFrames && (jr.fixer[1] != NULL || curr->endFixed(1))) { // constraining curve, or nudged
	FixedControl fcp;
	fcp._t = j;  fcp._n = numControls-1; fcp._c = i;
	mts->_fixedControls.push_back(fcp);
      }

      which = _joints[i].which[1];
      if (which==-1) {// resident last point
	mts->_Z.push_back(*(curr->getCtrl(numControls-1)));
	spline->ZMap(numControls-1) = vc++;
	if (j>0 && j<mts->_numFrames)
	  spline->varMap(numControls-1) = var++;
      }
      else {
	spline->ZMap(numControls-1) = mts->end_tc(j,which,_joints[i].side[1]);
	if (j>0 && j<mts->_numFrames)
	  spline->varMap(numControls-1) = mts->end_tc_var(j,which,_joints[i].side[1]);
      }

      curr = curr->nextC(); 


      if (jr.fixer[0] != NULL)
	jr.fixer[0] = jr.fixer[0]->nextC();
      if (jr.fixer[1] != NULL)
	jr.fixer[1] = jr.fixer[1]->nextC();
    }  // loop over time
  }  // loop over curves

  assert(vc==(int)mts->_Z.size());


  //assert(vc % (mts->_numFrames+1) == 0);
  //mts->_nVars = vc / (mts->_numFrames+1);

  mts->finishInit();

}


void TrackGraph::copyLocs(MultiSplineData* mt) { 
  int c,t, numEm;
  PathV currPaths = _key0paths;
  RotoPath *curr;

  for (c=0; c<mt->_nCurves; ++c) {
    curr = _key0paths[c]->nextC();
    //printf("curr has %d controls\n",curr->getNumControls());
    numEm = _key0paths[c]->getNumControls();

    for (t=1; t<mt->_numFrames; ++t) {  
      curr->takeBez(mt->getSpline(c,t));
      //for (n=0; n<numEm; ++n)
      //curr->setControl(mt->_Z[mt->tcn(t,c,n)], n);

      //curr->handleNewBezCtrls();
      curr = curr->nextC();

    } // t
    
  } // c

}

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



#include "liveWire.h"

LiveWire::LiveWire(Vec2i seedPoint, ImageGraph *graph) : 
  _graph(graph) {

  _watch = false;
  _curr = NULL;
  _w = _graph->width();
  _h = _graph->height();
  seedPoint.Clamp(0,_w-1, _h-1);
  assert(seedPoint.x()>=0 && seedPoint.x() < _w);
  assert(seedPoint.y()>=0 && seedPoint.y() < _h);

  // set up minimimum GraphNodePtr for Heap
  GraphNode minDummy;
  minDummy._totalCost = -1;
  GraphNodePtr minDummyPtr(&minDummy);

  // now follows Zhang pseudocode exactly
  pq = new BinaryHeap<GraphNodePtr>(minDummyPtr);
  _graph->init();  
  _graph->getNode(seedPoint.x(),seedPoint.y())->_prevNode = NULL;
  _graph->getNode(seedPoint.x(),seedPoint.y())->_totalCost = 0;
  GraphNodePtr newPtr(_graph->getNode(seedPoint.x(),seedPoint.y()));
  pq->Insert(newPtr);
  
  // while (continueIniting()) {} // to get rid of qt, put this back in
  
  continueIniting();
  _timer = new QTimer(this);
  connect(_timer, SIGNAL(timeout()), this,SLOT(continueIniting()));
  _timer->start(0);

}

bool LiveWire::continueIniting() {
  if (pq->IsEmpty()) {
    _timer->stop();
    return false;
  }
  
  GraphNode* prevVal;
  if (_watch)
    prevVal = _graph->getNode(_watchLoc.x(), _watchLoc.y())->_prevNode;
  

  for (int count=0; count<10000 && !pq->IsEmpty(); ++count) {

  GraphNodePtr q;
  int i,j, di, dj, lnum, offset;
  pq->DeleteMin(q);
  q._ptr->_state = EXPANDED;
  i = q._ptr->_loc.x();
  j = q._ptr->_loc.y();
  offset=(j-1)*_w+i-1;
  
  // iterate over 8-neighbors
  for (dj=-1; dj<=1; ++dj, offset+=_w-3) {
    for (di=-1; di<=1; ++di, ++offset) {
    if (j+dj<0 || j+dj>=_h)
      continue;
      if ((dj==0 && di==0) ||
	  (i+di<0 || i+di>=_w))
	continue;
      
      // establish link number from q to r
      if (dj==-1)
	lnum = -di+2;
      else if (dj==0) {
	if (di==-1) lnum=4;
	else lnum=0;
      }
      else
	lnum = di+6;
      
      GraphNode* r = _graph->getNode(offset);
      float newCost = q._ptr->_totalCost + q._ptr->_linkCost[lnum];
      if (r->_state != EXPANDED) {
	if (r->_state == INITIAL) {
	  r->_prevNode = q._ptr;
	  r->_totalCost = newCost;
	  r->_state = ACTIVE;
	  GraphNodePtr rPtr(r);
	  pq->Insert(rPtr);
	}
	else {
	  assert(r->_state == ACTIVE);
	  if ( newCost < r->_totalCost) {
	    r->_prevNode = q._ptr;
	    r->_totalCost = newCost;
	  }
	}
      } // not expanded
    } // loop di
  } // loop dj
 
  }

  if (_watch && prevVal != _graph->getNode(_watchLoc.x(), _watchLoc.y())->_prevNode)
    emit interestChanged();
  return (!pq->IsEmpty());
 
}

LiveWire::~LiveWire() {
  _timer->stop();
  //delete _timer;
  delete pq;
}

void LiveWire::startTracePath(Vec2i freePoint) {  
  freePoint.Clamp(0,_w-1, _h-1);
  _curr = _graph->getNode(freePoint.x(),freePoint.y());
}

const Vec2i* LiveWire::prevStep() {   
  _curr = _curr->_prevNode;
  if (_curr == NULL) return NULL;
  else
    return &(_curr->_loc);
}

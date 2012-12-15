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



#ifndef IMAGEGRAPH_H
#define IMAGEGRAPH_H

#include <qimage.h>
#include "jl_vectors.h"

enum NodeState {INITIAL, ACTIVE, EXPANDED};

class GraphNode {

 public:

  GraphNode() {
    _prevNode = NULL;
    _state = INITIAL;
  }

  float _linkCost[8];
  NodeState _state;
  float _totalCost;
  GraphNode* _prevNode;
  Vec2i _loc;
};

class GraphNodePtr {
 public:
  
  GraphNodePtr(GraphNode* i) { _ptr = i;}
  GraphNodePtr() {}
  
  GraphNodePtr operator=(const GraphNodePtr& other) {
    _ptr = other._ptr; 
  return *this;
  }
  
  int operator<(const GraphNodePtr& other) {
    return (_ptr->_totalCost < other._ptr->_totalCost);
  }

  GraphNode* _ptr;
};

class ImageGraph {


 public:
  
  ImageGraph(const QImage im);
  ImageGraph(FILE* fp);
  ~ImageGraph() { delete[] _graph; }

  void init();

  int width() { return _w;}
  int height() { return _h;}

  GraphNode* getNode(const int i, const int j) {
    assert(i>=0 && i<_w && j>=0 && j<_h);
    return _graph+j*_w+i;
  }

  GraphNode* getNode(const int i) {
    assert(i>=0 && i<_maxI); // SPEED
    return _graph+i;
  }

  void save(const char* fileName) const;

 private:
  
  Vec3i paddedImageQuery(const QImage im, int i, int j);
  float horizLink(const QImage im, const int x1, const int y1, 
		   const int x2, const int y2,
		   const int x3, const int y3, 
		   const int x4, const int y4);
  float diagLink(const QImage im, const int x1, const int y1, 
		   const int x2, const int y2);
 
  int _w, _h, _maxI;
  GraphNode* _graph;

};

#endif

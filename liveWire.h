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



#ifndef LIVEWIRE_H
#define LIVEWIRE_H

#include "jl_vectors.h"
#include "imageGraph.h"
#include "heap.h"
#include <qobject.h>
#include <qtimer.h>

class LiveWire : QObject {

  Q_OBJECT
  
 public:
  
  LiveWire(Vec2i seedPoint, ImageGraph *graph);
  ~LiveWire();


  // Iterator
  void startTracePath(Vec2i freePoint);
  const Vec2i* prevStep();
  
  void registerInterest(Vec2i loc, QObject* obj) {
    _watch=true; _watchLoc = loc;
    connect(this, SIGNAL(interestChanged()), obj, SLOT(wireValChanged()));
  }

  public slots:
    
    bool continueIniting();

  signals:
    
    void interestChanged();

 private:
  ImageGraph* _graph;
  int _w, _h;
  QTimer* _timer;
  BinaryHeap<GraphNodePtr> *pq;
  bool _watch;
  Vec2i _watchLoc;

  // iterator
  GraphNode* _curr;

};



#endif

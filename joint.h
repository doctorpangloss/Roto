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



#ifndef JOINT_H
#define JOINT_H

#ifdef __APPLE__
using namespace std;
#endif
#include <vector>

// Each rotopath gets a joint object pointing to the other


class RotoPath;

class Joint {
 public:
  
  Joint() {_rp=NULL; }
  Joint(RotoPath* r, const int s) : _rp(r), _side(s) {}

  // this exists only to provide a strict weak ordering for <set>, 
  // so that it can sort the joints
  int operator<(const Joint& j) const { 
    return (_rp < j._rp); 
  }

  Joint(const Joint& j) { _rp=j._rp; _side=j._side; }
  bool operator==(const Joint& j) const { return (_rp==j._rp && _side==j._side); }

  bool rpEqual(const RotoPath* o) const { return o==_rp; }

  RotoPath* _rp;
  int _side; // 0 for beginning, 1 for end

};


typedef std::vector<Joint> Joints;

void makeJointsConsistent(Joints& jnts);

typedef std::set<Joint> JointSet;
JointSet JointSet_union(JointSet& js1, JointSet& js2);

bool jointsOk(const Joints& j);

#endif

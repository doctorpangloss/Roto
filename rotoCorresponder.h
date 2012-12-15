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



#ifndef ROTOCORRESPONDER_H
#define ROTOCORRESPONDER_H
//#include "abstractPath.h"
#include "geigerCorresponder.h"

class RotoCorresponder : public GeigerCorresponder {
 public:

  RotoCorresponder(AbstractPath* X, AbstractPath* Y) 
    : GeigerCorresponder(X, Y) {}

 protected:

  double linkCost(const int i, const int j, const int ip, const int jp) const;
  
  double internalCost(const int i, const int j) const;


};

#endif

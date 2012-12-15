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



/*********************************************************************
 * error.h
 *********************************************************************/

#ifndef _KLTERROR_H_
#define _KLTERROR_H_

#ifdef __cplusplus
extern "C" {
#endif 

#include <stdio.h>
#include <stdarg.h>

void KLTError(char *fmt, ...);
void KLTWarning(char *fmt, ...);

#ifdef __cplusplus
}
#endif 

#endif


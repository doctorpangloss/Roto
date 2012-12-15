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
 * base.h
 *********************************************************************/

#ifndef _BASE_H_
#define _BASE_H_

/*
#ifndef uchar
#define uchar unsigned char
#endif

#ifndef schar
#define schar signed char
#endif

#ifndef uint
#define uint unsigned int
#endif

#ifndef ushort
#define ushort unsigned short
#endif

#ifndef ulong
#define ulong unsigned long
#endif
*/

#define max(a,b)	((a) > (b) ? (a) : (b))
#define min(a,b)	((a) < (b) ? (a) : (b))
#define max3(a,b,c)	((a) > (b) ? max((a),(c)) : max((b),(c)))
#define min3(a,b,c)	((a) < (b) ? min((a),(c)) : min((b),(c)))

#endif


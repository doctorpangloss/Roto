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



#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

int main(int argc, char* argv[]) {

  int block_size = atoi(argv[1]), i;
  int _cs = atoi(argv[2]);
  int hold = atoi(argv[3]) - 1;
  hold = MIN(hold, block_size/2-1);   
  block_size = MAX(0, block_size - 2*hold); 
  int _block_size;
  int npb = (block_size-1) / MAX(1,(block_size-1) / _cs);
  _block_size = (block_size-1) / npb + 1;  
  assert(_block_size > 1);

  for (i=0; i<hold; i++)
    printf("%d: %d %.5f\n",i,i,1.);

  int left = (block_size-1) % npb;       // these pts must be distributed into last buckets
  printf("%d: %d %.5f\n",hold,hold,1.);
  int bucket=0, ct=1;                    // which bucket we're dumping into, # points in current bucket
  double k;
  for (i=1; i<block_size; i++) {
    if (bucket >= (_block_size-left-1)) {  // into last 'left' buckets
      if (ct >= npb+1) {
	++bucket;
	ct=0;
	k=1.;
      }
      else
	k = 1. - double(ct) / double(npb+1);
    }
    else {
      if (ct >= npb) {
	++bucket;
	ct=0;
	k=1.;
      }
      else
	k = 1. - double(ct) / double(npb);
    }
    
    printf("%d: %d %.5f\n",i+hold,bucket+hold,k);
    assert (bucket < _block_size);
	
    ++ct;
  }

  _block_size += 2*hold-1;  printf("later block size: %d\n", _block_size);
  for (i=0; i<hold; i++)
    printf("%d: %d %.5f\n",i+block_size+hold,i+bucket+hold+1,1.);

  assert(k==1.);
  assert(bucket = _block_size-1);

  return 1;
}

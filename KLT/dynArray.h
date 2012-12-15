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



#ifndef DYNARRAY_H
#define DYNARRAY_H

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

template<class T, int initSize> class DynArray {

  // overview: dynamically resizeable, growable array

public:

  // constructor
  DynArray() {
    count = 0;
    space = initSize;
    assert(space>0);
    data = new T[initSize]; }

  inline DynArray(const DynArray& other);
  inline DynArray<T, initSize>& operator=(const DynArray<T, initSize>& other);
  

  // destructor
  ~DynArray() { delete[] data; }

  int add(const T& p) {
    // Adds T p to array, returns indice T p is located at
    if (count == space) return hardAdd(p);
    data[count] = p;
    return count++;
 }

  T getElement(const int i) const {
    assert (i < count && i>-1);
    return data[i];
  }

  const T& getElementConstRef(const int i) const {
    assert (i < count && i>-1);
    return data[i];
  }

  int addDefault() {
    if (count == space) {
      int newSpace = (int)ceil(1.25 * double(space));  // get new size
      assert(newSpace>0);
      T* newEls = new T[newSpace];
      memcpy(newEls,data,space*sizeof(T)); // copy over from old array
      delete[] data;                          // delete old array
      data = newEls; 
      space = newSpace;
    }
    return count++;
  }

  const T* getPointerToElement(const int i) const {
    assert (i < count && i>-1);
    return (data+i);
  }

  T* getPointerToElement(const int i) {
    if (i < count && i>-1)
      return (data+i);
    else return NULL;
  }
  
  void setElement(const int i, const T& p) {
    assert (i < count && i>-1);
    data[i] = p;
  }

  int getNumElements() const {return count; }

  void resetElements() {count = 0; }

  void deleteLastElement() {
    if (count > 1)
      --count;
  }

  int findElement(const T& in) {
    for (int i=0; i<count; i++)
      if (data[i]==in) 
	return i;
    return -1;
  }

  // dangerous(changes integer pointers) and slow
  void removeElement(int i) {
    assert (i < count && i>-1);
    for (int k=i; k<count-1; k++)
      data[k] = data[k+1];
    --count;
  }

  //following is also dangerous, since data is not initialized
  void addNElements(int n) {
    if (count+n < space)
      count = count+n;
    else {
      int newSpace = count+n+1; //(int)ceil(1.25* double(count+n));
      T* newEls = new T[newSpace];
      memcpy(newEls,data,space*sizeof(T)); // copy over from old array
      delete[] data;                          // delete old array
      data = newEls;
      space = newSpace;
      assert(space>0);
      count += n;
    }
  }

  void ensureCapacity(int n) {
    if (space >= n) return;
    T* newEls = new T[n];
    memcpy(newEls,data,space*sizeof(T)); // copy over from old array
    delete[] data;                          // delete old array
    data = newEls;
    space = n;
  }
  
  // following use with caution
  T* getData() {return data; }
  void setCount(const int i) { count = i; }
  
  const T* getData() const {return data; }
  

  //DynArray<T,initSize>  *next, *prev;

protected:

  int hardAdd(const T& p) {
    int newSpace = (int)ceil(1.25 * double(space));  // get new size
    assert(newSpace > 0);
    T* newEls = new T[newSpace];
    memcpy(newEls,data,space*sizeof(T)); // copy over from old array
    delete[] data;                          // delete old array
    data = newEls; 
    data[count] = p;                        // put in new data
    space = newSpace;
    return count++;
  }
  
  T* data;
  int count;     // current number of datum
  int space;     // current size of array

};

template<class T, int initSize> inline DynArray<T, initSize>::DynArray(const DynArray& other) {
  //delete data;
  space = other.count;
  if (space == 0) 
    space = initSize;
  count = other.count;
  data = new T[space];
  memcpy(data,other.data, count*sizeof(T)); // don't care about (space-count) data
}

template<class T, int initSize> inline DynArray<T, initSize>& DynArray<T, initSize>::operator=(const DynArray<T, initSize>& other) {
  delete[] data;
  space = other.count;
  if (space == 0) 
    space = initSize;
  count = other.count;
  data = new T[space];
  memcpy(data,other.data, count*sizeof(T)); // don't care about (space-count) data
  return *this;
}


#endif

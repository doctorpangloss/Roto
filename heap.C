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






// Routine to allocate the heap array

template <class Etype>
void
BinaryHeap<Etype>::GetArray( int NewMaxSize )
{
  Array = new Etype [ NewMaxSize + 1 ];
}

// Constructor for BinaryHeap

template <class Etype>
BinaryHeap<Etype>::BinaryHeap( const Etype & MinVal ) :
  MaxSize( DefaultSize ), CurrentSize( 0 ), OrderOK( 1 )
{
  GetArray( MaxSize );
  Array[ 0 ] = MinVal;
}

// If heap is full, double heap array

template <class Etype>
void
BinaryHeap<Etype>::CheckSize( )
{
  if( CurrentSize == MaxSize )
    {
      Etype *Old = Array;
      GetArray( MaxSize * 2 );
      for( int i = 0; i <= MaxSize; i++ )
	Array[ i ] = Old[ i ];
      delete [ ] Old;
      MaxSize *= 2;
    }
}

// Add X into the heap without maintaining order

template <class Etype>
void
BinaryHeap<Etype>::Toss( Etype & X )
{
  CheckSize( );
  Array[ ++CurrentSize ] = X;
 
  //X.setHeapPtr(CurrentSize);
  if( X < Array[ CurrentSize / 2 ] )
    OrderOK = 0;
}

// Insert X into heap and if heap order is being maintained,
// percolate X up as needed

template <class Etype>
void
BinaryHeap<Etype>::Insert( Etype & X )
{
  if( OrderOK == 0 )
    {
      Toss( X );
      return;
    }
  
  CheckSize( );
  
  // Percolate up
  int Hole = ++CurrentSize;
  for( ; X < Array[ Hole / 2 ] && Hole>1; Hole /= 2 ) {
    Array[ Hole ] = Array[ Hole / 2 ];
    //Array[Hole].setHeapPtr(Hole);
  }
  Array[ Hole ] = X;
  //Array[Hole].setHeapPtr(Hole);
}

template <class Etype>
void BinaryHeap<Etype>::update(int i) {
  int Hole = i;
  Etype X = Array[i];
  for( ; X < Array[ Hole / 2 ] && Hole>1; Hole /= 2 ) {
    assert(Hole != 0);
    Array[ Hole ] = Array[ Hole / 2 ];
    //Array[Hole].setHeapPtr(Hole);
  }
  Array[ Hole ] = X;
  //Array[Hole].setHeapPtr(Hole);
  //checkHeapProperty(); // DEBUG
}

template <class Etype>
bool BinaryHeap<Etype>::checkHeapProperty() {
  for (int i=CurrentSize; i>1; i--)
    if (Array[i/2] > Array[i])
      return 0;
  return 1;
}

// Return minimum item in the heap
// Call FixHeap first if necessary

template <class Etype>
const Etype &
BinaryHeap<Etype>::FindMin( )
{
  if (IsEmpty( )) {
    printf("Binary heap is empty\n" );
    exit(0);
  }
  
  if( OrderOK == 0 )
    FixHeap( );
  return Array[ 1 ];
}

// Delete the minimum item and place it in X

template <class Etype>
void
BinaryHeap<Etype>::DeleteMin( Etype & X )
{
  X = FindMin( );
  //X.setHeapPtr(0);
  Array[ 1 ] = Array[ CurrentSize-- ];
  PercolateDown( 1 );
}

// Delete the minimum item; throw it away
// NOTE: It would be better to write an additional
// private member to consolidate the common work between
// the two forms of DeleteMin.

template <class Etype>
void
BinaryHeap<Etype>::DeleteMin( )
{
  Etype X;
  DeleteMin( X );
}

// Private member to percolate down in the heap

template <class Etype>
void
BinaryHeap<Etype>::PercolateDown( int Hole )
{
  int Child;
  Etype Tmp = Array[ Hole ];
  
  for( ; Hole * 2 <= CurrentSize; Hole = Child )
    {
      Child = Hole * 2;
      if( Child != CurrentSize &&
	  Array[ Child + 1 ] < Array[ Child ] )
	Child++;
      if( Array[ Child ] < Tmp ) {
	Array[ Hole ] = Array[ Child ];
	//Array[Hole].setHeapPtr(Hole);
      }
      else
	break;
    }
  Array[ Hole ] = Tmp;
  //Array[Hole].setHeapPtr(Hole);

}

// Linear time FixHeap member

template <class Etype>
void
BinaryHeap<Etype>::FixHeap( )
{
  for( int i = CurrentSize / 2; i > 0; i-- )
    PercolateDown( i );
  OrderOK = 1;
  //checkHeapProperty(); // DEBUG
}



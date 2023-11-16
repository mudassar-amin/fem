/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of TTL.
 *
 * TTL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * TTL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with TTL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using TTL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the TTL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#ifndef _API_H
#define _API_H

//===========================================================================
//
//  Application Programming Interface to TTL (API)
//
//===========================================================================


/** \page api Application Programming Interface to TTL (API)

There are two interface channels between TTL and the application data structure:
\ref darttypesec and \ref traitstypesec, which are template arguments to the TTL
functions.
Function templates in TTL use either \c DartType or both \c DartType and \c TraitsType
depending on the complexity of the algorithms.
\c DartType is used for topological queries, and \c TraitsType is
used for topological modifiers and geometric calculations.
The documentation in TTL tells which functionality
is needed in \c DartType and \c TraitsType for each function template.

TTL consists of two namespaces with generic functions: \ref ttl contains
the main interface and \ref ttl_util contains utility functions that can also be used
by the application programmer. The figure below shows the (simple) architecture
of an application using TTL.


\image html api.gif


In the following we explain how \c DartType and \c TraitsType
must be implemented as an interface to the actual data structure.
For more details see \ref publications "literature".

\section darttypesec DartType

The generic functions in TTL navigates in the application data structure through a
"topological element" call a \e dart, which must be implemented as a struct or a class by the
application programmer. A dart can be considered as a unique triple \e d = (\e V, \e E, \e T),
where \e V is a node of the edge \e E, and \e V and \e E is a node and an edge of the
triangle \e T; see figure a) below where a dart is indicated as an arrow. 

    
\image html dart.gif


TTL expects that three functions, which we call \e alpha-iterators, are present in the
dart class: \c alpha0(), \c alpha1() and \c alpha2(). These functions reposition the
dart in the triangulation as shown in figure b) above and return a reference to the modified
dart itself. Thus, \c alpha0(), \c alpha1() and \c alpha2() change the node, the edge and the
triangle of the triple \e d = (\e V, \e E, \e T) respectively.

\note
- The \e dart is not part of the data structure. It is a "dynamic"
  element that "moves" around in the actual application data structure when
  the alpha-iterators are applied to it. Thus, this mechanism does
  not require any extra memory other than that occupied by the (few) darts
  that are involved in a function template of TTL that is called.

In addition to the alpha-iterators, TTL expects that standard class member functions such
as constructors, assignment operators and the like are also implemented
in the dart class.

The following syntax is required:

\code
class MyDart {
  ...
public:
  // Constructors and destructors
  ...
  MyDart(const MyDart& dart) {...} // copy constructor
  MyDart& operator= (const MyDart& dart) {...} // assignment operator

  // comparing dart objects
  bool operator==(const MyDart& dart) const {...}
  bool operator!=(const MyDart& dart) const {...}

  // alpha-iterators
  MyDart& alpha0() {...; return *this;}
  MyDart& alpha1() {...; return *this;}
  MyDart& alpha2() {...; return *this;}
  };
\endcode

Thus, the alpha-iterators change the content of the dart and
return a reference to the dart itself.

\note
- Note that \c alpha2() must return the dart itself without changing its content
  (when checked with the \c == operator) if the edge associated with the dart
  is at the boundary of the triangulation. The source code of
  \ref ttl::isBoundaryEdge illustrates this:

\code
template <class DartType>
  bool isBoundaryEdge(const DartType& dart) {
    
  DartType dart_iter = dart;
  if (dart_iter.alpha2() == dart)
    return true;
  else
    return false;
}
\endcode

Consult the definition and source code of class \ref hed::Dart, which implements the dart class
for the half-edge data structure, for a thorough example.

\section traitstypesec TraitsType

\c TraitsType can be a static class or a struct that contains functionality
required by the functions in TTL. The purpose with \c TraitsType is twofold:
- Let the application programmer provide topological modifiers
  on the actual data structure that are not resolved by TTL,
  e.g., remove a triangle from the triangulation,
- Let the application programmer provide basic geometric calculations and thus
  control the level of accuracy for such calculations.

For example, assume that a function in TTL requires a scalar product between
vectors in the plane, represented as darts, with the following syntax:

\code 
real_type scalarProduct2d(const DartType& d1, const DartType& d2)
\endcode

Then this function must be implemented in the traits class together
with the definition of \c real_type:

\code
struct MyTraitsType {
  typedef double real_type;
  ...
  static real_type scalarProduct2d(const MyDart& v1, const MyDart& v2) {
    return ::my_scalarProduct2d(v1,v2);
  }
  ...
  static bool swapTestDelaunay(const MyDart& dart) {
    if (dart.getEdge().isConstrained())
      return false;
    else
      return ttl::swapTestDelaunay<MyTraitsType>(dart);
  }
};
\endcode
  
The second function, \c swapTestDelaunay, implements the Delaunay swapping
test. The function desides if the edge associated with the given dart
should be swapped according to the \e circumcircle criterion
(or equivalently according to the \e MaxMin angle criterion).
In the example above, it is first examined if the edge is constrained, in which case
false is returned. Then the test is directed back to TTL again since the
circumcircle test is present there; see ttl::swapTestDelaunay.
Of course, the user could also make his/her own implementation of the
circumcircle test, e.g., by using robust schemes with exact arithmetic
as Jonathan Richard Shewchuk does (http://www-2.cs.cmu.edu/~quake/robust.html).

This example also illustrates the syntax
\c functionName<MyTraitsType>(...arguments...) when calling
function templates using a traits class.

Actually, scalar product between vectors, as well as cross product and other
point and vector algebra, are also present in the the namespace \ref ttl_util.
Assume that \c MyDart::x() and \c MyDart::y() deliver the coordinates in the
plane of the node associated with a dart. Then the scalar product
in \c MyTraitsType can be implementes as:

\code
static real_type scalarProduct2d(const MyDart& v1, const MyDart& v2) {
  MyDart v10 = v1; v10.alpha0();
  Dart v20 = v2; v20.alpha0();
  return ttl_util::scalarProduct2d(v10.x()-v1.x(), v10.y()-v1.y(), 
                                   v20.x()-v2.x(), v20.y()-v2.y());
}
\endcode

Consult the definition and source code of struct \ref hed::TTLtraits,
which implements the traits class for the half-edge data structure.
This example shows which members that are required in the traits class
for running incremental Delaunay triangulation and removing nodes in
a triangulation.

*/

# endif // _API_H

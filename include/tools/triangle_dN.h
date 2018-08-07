/*
 
 2018, Oded Stein, Eitan Grinspun and Keenan Crane
 
 This file is part of the code for "Developability of Triangle Meshes".
 
 The code for "Developability of Triangle Meshes" is free software: you can
 redistribute it and/or modify it under the terms of the GNU General Public
 License as published by the Free Software Foundation, either version 2 of the
 License, or (at your option) any later version.
 
 The code for "Developability of Triangle Meshes" is distributed in the hope
 that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with the code for "Developability of Triangle Meshes". If not,
 see <https://www.gnu.org/licenses/>.
 
 */


#ifndef DEVELOPABLEFLOW_TRIANGLE_DN
#define DEVELOPABLEFLOW_TRIANGLE_DN

#include <Eigen/Core>

#include <igl/igl_inline.h>

//Returns the derivatives of the normals with respect to every vertex

template <typename derivedV, typename derivedScalar, typename derivedMat>
IGL_INLINE void triangle_dN(
                                           const Eigen::PlainObjectBase<derivedV>& eij, //edge eij
                                           const Eigen::PlainObjectBase<derivedV>& ejk, //edge ejk
                                           const Eigen::PlainObjectBase<derivedV>& eki, //edge eki
                                           const Eigen::PlainObjectBase<derivedV>& N, //normal
                                           const derivedScalar& dA, //double area
                                           Eigen::PlainObjectBase<derivedMat>& dNdi, //derivative dNdi
                                           Eigen::PlainObjectBase<derivedMat>& dNdj, //derivative dNdj
                                           Eigen::PlainObjectBase<derivedMat>& dNdk); //derivative dNdk



#ifndef IGL_STATIC_LIBRARY
#  include "triangle_dN.cpp"
#endif

#endif

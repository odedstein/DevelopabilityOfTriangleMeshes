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


#ifndef DEVELOPABLEFLOW_UNITVEC_DERIVATIVE
#define DEVELOPABLEFLOW_UNITVEC_DERIVATIVE

#include <Eigen/Core>

#include <igl/igl_inline.h>

//Returns the derivatives of the normalized edge with respect to every vertex

template <typename derivedV, typename derivedMat>
IGL_INLINE void unitvec_derivative(
                                const Eigen::PlainObjectBase<derivedV>& e, //vector eij (not normalized)
                                Eigen::PlainObjectBase<derivedMat>& dedTip, //derivative dedTip
                                Eigen::PlainObjectBase<derivedMat>& dedTail); //derivative dedTail



#ifndef IGL_STATIC_LIBRARY
#  include "unitvec_derivative.cpp"
#endif

#endif

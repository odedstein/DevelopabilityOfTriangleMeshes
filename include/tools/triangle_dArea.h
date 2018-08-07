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


#ifndef DEVELOPABLEFLOW_TRIANGLE_DAREA
#define DEVELOPABLEFLOW_TRIANGLE_DAREA

#include <Eigen/Core>

#include <igl/igl_inline.h>

//Returns the derivatives of the DOUBLE AREAS with respect to vertex positions (important! double areas!)

template <typename derivedV, typename derivedMat>
IGL_INLINE void triangle_dArea(
                               const Eigen::PlainObjectBase<derivedV>& vi, //position of vert i
                               const Eigen::PlainObjectBase<derivedV>& vj, //position of vert j
                               const Eigen::PlainObjectBase<derivedV>& vk, //position of vert k
                               const Eigen::PlainObjectBase<derivedV>& N, //face normal
                               Eigen::PlainObjectBase<derivedMat>& dAdi, //derivative dAdi
                               Eigen::PlainObjectBase<derivedMat>& dAdj, //derivative dAdj
                               Eigen::PlainObjectBase<derivedMat>& dAdk); //derivative dAdk



#ifndef IGL_STATIC_LIBRARY
#  include "triangle_dArea.cpp"
#endif

#endif

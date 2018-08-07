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


#ifndef DEVELOPABLEFLOW_TRIANGLE_COSDIHEDRAL_ANGLE_H
#define DEVELOPABLEFLOW_TRIANGLE_COSDIHEDRAL_ANGLE_H

#include <igl/igl_inline.h>

#include <Eigen/Core>

template <typename derivedIC, typename derivednormals, typename derivedTT, typename deriveddihedral>
IGL_INLINE void triangle_cosdihedral_angle(
                                           const Eigen::PlainObjectBase<derivedIC>& IC, //IC from unique_simplices
                                           const Eigen::PlainObjectBase<derivedTT>& TT, //TT from triangle_triangle_adjacency
                                           const Eigen::PlainObjectBase<derivednormals>& normals, //per-face normals
                                           Eigen::PlainObjectBase<deriveddihedral>& dihedral); //cos of dihedral angle at each edge


#ifndef IGL_STATIC_LIBRARY
#  include "triangle_cosdihedral_angle.cpp"
#endif

#endif

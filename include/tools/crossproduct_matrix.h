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


#ifndef CROSSPRODUCT_MATRIX_H
#define CROSSPRODUCT_MATRIX_H

#include <Eigen/Core>

#include <igl/igl_inline.h>

template <typename derivedMatrix, typename derivedVector>
IGL_INLINE void crossproduct_matrix(const Eigen::PlainObjectBase<derivedVector>& v, //vector
                                    Eigen::PlainObjectBase<derivedMatrix>& J); //matrix return value


#ifndef IGL_STATIC_LIBRARY
#  include "crossproduct_matrix.cpp"
#endif

#endif

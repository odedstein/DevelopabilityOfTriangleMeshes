#ifndef IGL_KOPP_H
#define IGL_KOPP_H

#include <Eigen/Core>

#include <igl/igl_inline.h>

// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------

template<typename derivedScalar>
IGL_INLINE void eigendecomp(const Eigen::PlainObjectBase<derivedScalar>& A,
                            Eigen::PlainObjectBase<derivedScalar>& D,
                            Eigen::PlainObjectBase<derivedScalar>& V);

template<typename Scalar>
IGL_INLINE int dsyevc3(Scalar A[3][3], Scalar w[3]);

template<typename Scalar>
IGL_INLINE int dsyevh3(Scalar A[3][3], Scalar Q[3][3], Scalar w[3]);

template<typename Scalar>
IGL_INLINE int dsyevq3(Scalar A[3][3], Scalar Q[3][3], Scalar w[3]);

template<typename Scalar>
IGL_INLINE void dsytrd3(Scalar A[3][3], Scalar Q[3][3], Scalar d[3], Scalar e[2]);

#ifndef IGL_STATIC_LIBRARY
#  include "kopp.cpp"
#endif

#endif

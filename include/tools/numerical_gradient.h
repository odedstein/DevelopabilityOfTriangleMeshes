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


#ifndef NUMERICAL_GRADIENT_H
#define NUMERICAL_GRADIENT_H

#include <functional>

#include <Eigen/Core>

#include <igl/igl_inline.h>

//Computes the numerical gradient of a function

template <typename energyScalar, typename derivedPosition, typename stepsizeScalar, typename derivedGradient>
IGL_INLINE void numerical_gradient(const std::function<energyScalar(const Eigen::PlainObjectBase<derivedPosition>&)>& energyFct, //vector
                                   const Eigen::PlainObjectBase<derivedPosition>& currentPosition, //the current position at which to evaluate the gradient
                                   const stepsizeScalar& stepsize, //numerical gradient step size
                                   Eigen::PlainObjectBase<derivedGradient>& numericalGradient); //matrix return value


#ifndef IGL_STATIC_LIBRARY
#  include "numerical_gradient.cpp"
#endif

#endif

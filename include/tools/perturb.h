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


#ifndef DEVELOPABLEFLOW_PERTURB_H
#define DEVELOPABLEFLOW_PERTURB_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>
#include <string>


template <typename derivedV, typename indexType, typename derivedScalar>
IGL_INLINE void perturb(
                        Eigen::PlainObjectBase<derivedV>& V, //Vertices
                        const std::vector<std::vector<indexType> >& VF, //Vertex face neighbors
                        const std::vector<bool>& isB, //Boundary vertices
                        const derivedScalar& strength = 0.1, //strength of perturbation
                        const int& randomSeed = 0); //random seed for perturbation



#ifndef IGL_STATIC_LIBRARY
#  include "perturb.cpp"
#endif

#endif

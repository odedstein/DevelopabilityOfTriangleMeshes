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


#ifndef DEVELOPABLEFLOW_MAX_HINGE_ENERGY_H
#define DEVELOPABLEFLOW_MAX_HINGE_ENERGY_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>

template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void max_hinge_energy_and_grad(
                                      const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                      const Eigen::PlainObjectBase<derivedF>& F, //TT from triangle_triangle_adjacency
                                      const std::vector<std::vector<indexType> >& VF, //VF from vertex-triangle adjacency
                                      const std::vector<std::vector<indexType> >& VFi, //VFi from vertex-triangle adjacency
                                      const std::vector<bool>& isB, //isB from is_border_vertex
                                      Eigen::PlainObjectBase<derivedEnergy>& energy, //energy return val
                                      Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad); //energy grad return val


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy>
IGL_INLINE void max_hinge_energy(
                             const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                             const Eigen::PlainObjectBase<derivedF>& F, //TT from triangle_triangle_adjacency
                             const std::vector<std::vector<indexType> >& VF, //VF from vertex-triangle adjacency
                             const std::vector<std::vector<indexType> >& VFi, //VFi from vertex-triangle adjacency
                             const std::vector<bool>& isB, //isB from is_border_vertex
                             Eigen::PlainObjectBase<derivedEnergy>& energy); //energy return val

#ifndef IGL_STATIC_LIBRARY
#  include "max_hinge_energy.cpp"
#endif

#endif

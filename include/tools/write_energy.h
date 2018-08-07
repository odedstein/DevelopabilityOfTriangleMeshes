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


#ifndef DEVELOPABLEFLOW_WRITE_ENERGY_H
#define DEVELOPABLEFLOW_WRITE_ENERGY_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>
#include <string>


template <typename derivedV, typename derivedF, typename indexType, typename derivedP, typename t_res>
IGL_INLINE void write_energy(
                            const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                            const Eigen::PlainObjectBase<derivedF>& F, //Faces
                            const std::vector<std::vector<indexType> >& VF, //VF from vertex-triangle adjacency
                            const std::vector<std::vector<indexType> >& VFi, //VFi from vertex-triangle adjacency
                            const std::vector<bool>& isB, //isB from is_border_vertex
                            const Eigen::PlainObjectBase<derivedP>& p, //direction along which energy will be sampled
                            const std::string& filename, //file to write the energy values to
                            const t_res& resolution = 1e-7, //resolution of the energy writing
                            EnergyType energyType = ENERGY_TYPE_HINGE);//Which energy to use for the step: 0 normal, 1 midWidth (according to order in this file)



#ifndef IGL_STATIC_LIBRARY
#  include "write_energy.cpp"
#endif

#endif

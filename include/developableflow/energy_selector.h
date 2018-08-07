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


#ifndef DEVELOPABLEFLOW_ENERGY_SELECTOR_H
#define DEVELOPABLEFLOW_ENERGY_SELECTOR_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>

//Selects the correct energy and applies regularizers


enum EnergyType {
    ENERGY_TYPE_HINGE = 0,
    ENERGY_TYPE_MINWIDTH = 1,
    ENERGY_TYPE_OLDHINGE = 2,
    ENERGY_TYPE_OLDMINWIDTH = 3,
    ENERGY_TYPE_PAIRWISENORMALS = 4,
    ENERGY_TYPE_MAXPAIRWISENORMALS = 5,
    ENERGY_TYPE_NUMS = 6
};


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void energy_selector(const EnergyType energyType, //Which energy type to use
                                const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                const Eigen::PlainObjectBase<derivedF>& F, //TT from triangle_triangle_adjacency
                                const std::vector<std::vector<indexType> >& VF, //VF from vertex-triangle adjacency
                                const std::vector<std::vector<indexType> >& VFi, //VFi from vertex-triangle adjacency
                                const std::vector<bool>& isB, //isB from is_border_vertex
                                Eigen::PlainObjectBase<derivedEnergy>& energy, //energy return val
                                Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad); //energy grad return val



#ifndef IGL_STATIC_LIBRARY
#  include "energy_selector.cpp"
#endif

#endif

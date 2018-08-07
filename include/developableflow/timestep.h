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


#ifndef DEVELOPABLEFLOW_TIMESTEP_H
#define DEVELOPABLEFLOW_TIMESTEP_H

#include <igl/igl_inline.h>
#include "energy_selector.h"

#include <Eigen/Core>
#include <vector>

//Perform one timestep in the direction p, starting with the guess t. At the end, t will contain the actual length of the timestep, and p will contain a new direction. V will contain a new position, energy and energygrad will contain energy and gradient at the new place.
//mode is the type of line search that should be done, step is the type of stepping that should be done
//Returns 0 if timestep completes, returns -1 if linesearch exceeded maxium tries

enum Linesearch {
    LINESEARCH_NONE = 0,
    LINESEARCH_OVERTON = 1,
    LINESEARCH_BACKTRACK = 2,
    LINESEARCH_NUMS = 3
};

enum StepType {
    STEP_TYPE_GRADDESC = 0,
    STEP_TYPE_LBFGS = 1,
    STEP_TYPE_NUMS = 2
};


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad, typename derivedT, typename derivedP>
IGL_INLINE int hinge_timestep(
                        Eigen::PlainObjectBase<derivedV>& V, //Vertices
                        const Eigen::PlainObjectBase<derivedF>& F, //TT from triangle_triangle_adjacency
                        const std::vector<std::vector<indexType> >& VF, //VF from vertex-triangle adjacency
                        const std::vector<std::vector<indexType> >& VFi, //VFi from vertex-triangle adjacency
                        const std::vector<bool>& isB, //isB from is_border_vertex
                        derivedT& t, //initial time guess, contains actual time step at the end
                        Eigen::PlainObjectBase<derivedP>& p, //search direction, will be updated with new search direction in the end
                        Eigen::PlainObjectBase<derivedEnergy>& energy, //energy return val
                        Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad, //energy grad return val
                        Linesearch mode = LINESEARCH_NONE, //the type of line search to use
                        StepType type = STEP_TYPE_GRADDESC); //Which step method to use.

template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad, typename derivedT, typename derivedP>
IGL_INLINE int minwidth_timestep(
                                 Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                 const Eigen::PlainObjectBase<derivedF>& F, //TT from triangle_triangle_adjacency
                                 const std::vector<std::vector<indexType> >& VF, //VF from vertex-triangle adjacency
                                 const std::vector<std::vector<indexType> >& VFi, //VFi from vertex-triangle adjacency
                                 const std::vector<bool>& isB, //isB from is_border_vertex
                                 derivedT& t, //initial time guess, contains actual time step at the end
                                 Eigen::PlainObjectBase<derivedP>& p, //search direction, will be updated with new search direction in the end
                                 Eigen::PlainObjectBase<derivedEnergy>& energy, //energy return val
                                 Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad, //energy grad return val
                                 Linesearch mode = LINESEARCH_NONE, //the type of line search to use
                                 StepType type = STEP_TYPE_GRADDESC); //Which step method to use.


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad, typename derivedT, typename derivedP>
IGL_INLINE int timestep(
                             Eigen::PlainObjectBase<derivedV>& V, //Vertices
                             const Eigen::PlainObjectBase<derivedF>& F, //TT from triangle_triangle_adjacency
                             const std::vector<std::vector<indexType> >& VF, //VF from vertex-triangle adjacency
                             const std::vector<std::vector<indexType> >& VFi, //VFi from vertex-triangle adjacency
                             const std::vector<bool>& isB, //isB from is_border_vertex
                             derivedT& t, //initial time guess, contains actual time step at the end
                             Eigen::PlainObjectBase<derivedP>& p, //search direction, will be updated with new search direction in the end
                             Eigen::PlainObjectBase<derivedEnergy>& energy, //energy return val
                             Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad, //energy grad return val
                             Linesearch mode = LINESEARCH_NONE, //the type of line search to use
                             StepType type = STEP_TYPE_GRADDESC, //Which step method to use.
                        EnergyType energyType = ENERGY_TYPE_HINGE); //Which energy to use for the step: 0 normal, 1 midWidth (according to order in this file)



#ifndef IGL_STATIC_LIBRARY
#  include "timestep.cpp"
#endif

#endif

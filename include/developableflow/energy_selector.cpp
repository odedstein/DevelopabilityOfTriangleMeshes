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


#include "energy_selector.h"

#include "hinge_energy.h"
#include "max_hinge_energy.h"
#include "old_hinge_energy.h"
#include "old_max_hinge_energy.h"
#include "hingepairs_energy.h"
#include "maxhingepairs_energy.h"

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>


//Infinity
#define INFTY std::numeric_limits<double>::infinity()


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void energy_selector(const EnergyType energyType,
                                const Eigen::PlainObjectBase<derivedV>& V,
                                const Eigen::PlainObjectBase<derivedF>& F,
                                const std::vector<std::vector<indexType> >& VF,
                                const std::vector<std::vector<indexType> >& VFi,
                                const std::vector<bool>& isB,
                                Eigen::PlainObjectBase<derivedEnergy>& energy,
                                Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad)
{
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergyGrad::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef Eigen::SparseMatrix<t_V_s> t_Vs;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    
    switch(energyType) {
        case ENERGY_TYPE_HINGE:
            hinge_energy_and_grad(V, F, VF, VFi, isB, energy, energyGrad);
            break;
        case ENERGY_TYPE_MINWIDTH:
            max_hinge_energy_and_grad(V, F, VF, VFi, isB, energy, energyGrad);
            break;
        case ENERGY_TYPE_OLDHINGE:
            old_hinge_energy_and_grad(V, F, VF, VFi, isB, energy, energyGrad);
            break;
        case ENERGY_TYPE_OLDMINWIDTH:
            old_max_hinge_energy_and_grad(V, F, VF, VFi, isB, energy, energyGrad);
            break;
        case ENERGY_TYPE_PAIRWISENORMALS:
            hingepairs_energy_and_grad(V, F, VF, VFi, isB, energy, energyGrad);
            break;
        case ENERGY_TYPE_MAXPAIRWISENORMALS:
            maxhingepairs_energy_and_grad(V, F, VF, VFi, isB, energy, energyGrad);
            break;
        default:
            std::cout << "Such an energy type does not exist." << std::endl;
    }
    
    
    assert(energy==energy && "There are nans in the energy");
    assert(energyGrad==energyGrad && "There are nans in the energyGrad");
}

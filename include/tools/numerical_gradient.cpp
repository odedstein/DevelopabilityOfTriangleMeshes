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


#include "numerical_gradient.h"

template <typename energyScalar, typename derivedPosition, typename stepsizeScalar, typename derivedGradient>
IGL_INLINE void numerical_gradient(const std::function<energyScalar(const Eigen::PlainObjectBase<derivedPosition>&)>& energyFct,
                                   const Eigen::PlainObjectBase<derivedPosition>& currentPosition,
                                   const stepsizeScalar& stepsize,
                                   Eigen::PlainObjectBase<derivedGradient>& numericalGradient)
{
    
    const energyScalar currentEnergy = energyFct(currentPosition);
    
    Eigen::PlainObjectBase<derivedPosition> updatedPositionFwd = currentPosition;
    Eigen::PlainObjectBase<derivedPosition> updatedPositionBwd = currentPosition;
    numericalGradient.resize(currentPosition.rows(), currentPosition.cols());
    for(int ind=0; ind<currentPosition.size(); ++ind) {
        if(ind!=0) {
            updatedPositionFwd(ind-1) = currentPosition(ind-1);
            updatedPositionBwd(ind-1) = currentPosition(ind-1);
        }
        updatedPositionFwd(ind) += stepsize;
        const energyScalar newEnergyFwd = energyFct(updatedPositionFwd);
        updatedPositionBwd(ind) -= stepsize;
        const energyScalar newEnergyBwd = energyFct(updatedPositionBwd);
        
        //Forward difference
        //numericalGradient(ind) = (newEnergyFwd - currentEnergy)/stepsize;
        //Central difference
        numericalGradient(ind) = (newEnergyFwd - newEnergyBwd)/(2.*stepsize);
    }
    
}

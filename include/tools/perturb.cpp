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


#include "perturb.h"

#include <stdlib.h>

#include <Eigen/Core>


template <typename derivedV, typename indexType, typename derivedScalar>
IGL_INLINE void perturb(
                        Eigen::PlainObjectBase<derivedV>& V,
                        const std::vector<std::vector<indexType> >& VF,
                        const std::vector<bool>& isB,
                        const derivedScalar& strength,
                        const int& randomSeed)
{
    std::srand(randomSeed);
    
    for(int i=0; i<V.rows(); ++i) {
#ifdef IGNORE_VALENCE_3
        if(VF[i].size()<4)
            continue;
#endif
#ifdef FIX_BOUNDARIES
        if(isB[i])
            continue;
#endif
        
        for(int j=0; j<3; ++j)
            V(i,j) += strength*((double) std::rand()-0.5*RAND_MAX) / ((double) RAND_MAX);
    }
}

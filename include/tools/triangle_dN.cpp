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


#include "triangle_dN.h"

template <typename derivedV, typename derivedScalar, typename derivedMat>
IGL_INLINE void triangle_dN(
                            const Eigen::PlainObjectBase<derivedV>& eij,
                            const Eigen::PlainObjectBase<derivedV>& ejk,
                            const Eigen::PlainObjectBase<derivedV>& eki,
                            const Eigen::PlainObjectBase<derivedV>& N,
                            const derivedScalar& dA,
                            Eigen::PlainObjectBase<derivedMat>& dNdi,
                            Eigen::PlainObjectBase<derivedMat>& dNdj,
                            Eigen::PlainObjectBase<derivedMat>& dNdk)
{
    dNdi = ejk.cross(N)*N.transpose()/dA;
    dNdj = eki.cross(N)*N.transpose()/dA;
    dNdk = eij.cross(N)*N.transpose()/dA;
    
    assert(dNdi==dNdi);
    assert(dNdj==dNdj);
    assert(dNdk==dNdk);
}

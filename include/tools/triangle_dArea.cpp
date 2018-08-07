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


#include "triangle_dArea.h"


template <typename derivedV, typename derivedMat>
IGL_INLINE void triangle_dArea(
                               const Eigen::PlainObjectBase<derivedV>& vi,
                               const Eigen::PlainObjectBase<derivedV>& vj,
                               const Eigen::PlainObjectBase<derivedV>& vk,
                               const Eigen::PlainObjectBase<derivedV>& N,
                               Eigen::PlainObjectBase<derivedMat>& dAdi,
                               Eigen::PlainObjectBase<derivedMat>& dAdj,
                               Eigen::PlainObjectBase<derivedMat>& dAdk)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V;
    typedef typename derivedMat::Scalar t_Mat_s;
    typedef Eigen::Matrix<t_Mat_s, 3, 3> t_Mat;
    
    dAdi = N.cross(vk-vj);
    dAdj = N.cross(vi-vk);
    dAdk = N.cross(vj-vi);
}

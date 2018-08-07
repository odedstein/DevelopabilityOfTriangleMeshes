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


#include "crossproduct_matrix.h"


template <typename derivedMatrix, typename derivedVector>
IGL_INLINE void crossproduct_matrix(const Eigen::PlainObjectBase<derivedVector>& v,
                                    Eigen::PlainObjectBase<derivedMatrix>& J)
{
    J.resize(3,3);
    J << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
}

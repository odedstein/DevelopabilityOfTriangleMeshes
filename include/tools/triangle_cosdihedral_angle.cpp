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


#include "triangle_cosdihedral_angle.h"

template <typename derivedIC, typename derivednormals, typename derivedTT, typename deriveddihedral>
IGL_INLINE void triangle_cosdihedral_angle(
                                        const Eigen::PlainObjectBase<derivedIC>& IC,
                                        const Eigen::PlainObjectBase<derivedTT>& TT,
                                        const Eigen::PlainObjectBase<derivednormals>& normals,
                                        Eigen::PlainObjectBase<deriveddihedral>& dihedral)
{
    typedef Eigen::Matrix<typename derivedIC::Scalar, Eigen::Dynamic, 3> t_IC;
    typedef Eigen::Matrix<typename derivedTT::Scalar, Eigen::Dynamic, 3> t_TT;
    typedef Eigen::Matrix<typename derivednormals::Scalar, Eigen::Dynamic, 3> t_normals;
    typedef Eigen::Matrix<typename deriveddihedral::Scalar, Eigen::Dynamic, 1> t_dihedral;
    
    typedef Eigen::Matrix<typename derivednormals::Scalar, 3, 1> t_vec3;
    
    int numEdges = IC.maxCoeff()+1;
    int numFaces = TT.rows();
    
    dihedral = t_dihedral(numEdges);
    dihedral.setZero();
    
    for(int face=0; face<numFaces; ++face) {
        const t_vec3& n1 = normals.row(face);
        
        for(int j=0; j<3; ++j) {
            //j here is in normal edge convention. in triangle_triangle_adjacency convention, this is ttj
            const int ttj = (j+1)%3;
            const typename derivedTT::Scalar adjface = TT(face,ttj);
            //If this is a bdry edge, skip it
            if(adjface<0)
                continue;
            //Global edge index from IC
            const typename derivedIC::Scalar& edgeIndex = IC(j*numFaces+face);
            //dihedral angle reference
            typename deriveddihedral::Scalar& cosdihedralAngle = dihedral(edgeIndex);
            if(cosdihedralAngle==0) {
                //This edge has not been set yet. Set it.
                const t_vec3& n2 = normals.row(adjface);
                
                cosdihedralAngle = n1.dot(n2);
            }
        }
    }
    
}

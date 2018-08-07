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


#include "curvature_energy.h"

#include <tools/kopp.h>
#include <tools/triangle_dTheta.h>
#include <tools/crossproduct_matrix.h>

#include <igl/internal_angles.h>
#include <igl/per_face_normals.h>
#include <igl/squared_edge_lengths.h>
#include <igl/doublearea.h>
#include <igl/cotmatrix_entries.h>


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void curvature_energy_and_grad(
                                      const Eigen::PlainObjectBase<derivedV>& V,
                                      const Eigen::PlainObjectBase<derivedF>& F,
                                      const std::vector<std::vector<indexType> >& VF,
                                      const std::vector<std::vector<indexType> >& VFi,
                                      const std::vector<bool>& isB,
                                      Eigen::PlainObjectBase<derivedEnergy>& energy,
                                      Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3t;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergy::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    
    //Precomputation
    t_V edgesSq;
    igl::squared_edge_lengths(V, F, edgesSq);
    t_V edges = edgesSq.array().sqrt();
    t_V angles;
    igl::internal_angles_using_squared_edge_lengths(edgesSq, angles);
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
    t_Vv doubleAreas;
    igl::doublearea(edges, doubleAreas);
    
    energy = t_energy(V.rows());
    energyGrad = t_energyGrad::Zero(V.rows(), 3);
    for(int vert=0; vert<V.rows(); ++vert) {
        if(isB[vert]) {
            energy(vert) = 0;
            continue;
        }
        
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        t_V_s anglesum = 0.;
        t_V_s A = 0.;
        for(int f=0; f<adjacentFaces.size(); ++f) {
            const int& face = adjacentFaces[f];
            const t_F_i& i = vert;
            const t_F_i& j = F(face, (adjacentFacesi[f]+1)%3);
            const t_F_i& k = F(face, (adjacentFacesi[f]+2)%3);
            const t_V_s& theta = angles(face, adjacentFacesi[f]);
            
            if(doubleAreas(face)<1e-4)
                continue;
            
            //Add angle sum and barycentric area
            anglesum += theta;
            A += doubleAreas(face)/6.;
        }
        
        const t_V_s pi = 2*asin(1.);
        t_V_s K = 2.*pi - anglesum;
        
        energy(vert) = K*K/A;
        
        for(int f=0; f<adjacentFaces.size(); ++f) {
            const int& face = adjacentFaces[f];
            const t_F_i& i = vert;
            const t_F_i& j = F(face, (adjacentFacesi[f]+1)%3);
            const t_F_i& k = F(face, (adjacentFacesi[f]+2)%3);
            const t_V_s& theta = angles(face, adjacentFacesi[f]);
            const t_V3& normal = faceNormals.row(face);
            
            if(doubleAreas(face)<1e-4)
                continue;
            
            //Compute edges
            t_V3 eki = V.row(i) - V.row(k);
            t_V3 eij = V.row(j) - V.row(i);
            
            //Angle derivatives
            t_V3t dThetadi, dThetadj, dThetadk;
            triangle_dTheta(eij, eki, normal, dThetadi, dThetadj, dThetadk);
            
            //Edge derivatives
            t_V33 Jeij, Jeik;
            crossproduct_matrix(eij, Jeij);
            crossproduct_matrix((-eki).eval(), Jeik);
            const t_V33 dPdi = -Jeij+Jeik, dPdj = -Jeik, dPdk = Jeij;
            
            //Add together
            energyGrad.row(i) += -2*K/A*dThetadi - normal.transpose()*K*K/6/A/A*dPdi;
            energyGrad.row(j) += -2*K/A*dThetadj - normal.transpose()*K*K/6/A/A*dPdj;
            energyGrad.row(k) += -2*K/A*dThetadk - normal.transpose()*K*K/6/A/A*dPdk;
        }
    }

}


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy>
IGL_INLINE void curvature_energy(
                             const Eigen::PlainObjectBase<derivedV>& V,
                             const Eigen::PlainObjectBase<derivedF>& F,
                             const std::vector<std::vector<indexType> >& VF,
                             const std::vector<std::vector<indexType> >& VFi,
                             const std::vector<bool>& isB,
                             Eigen::PlainObjectBase<derivedEnergy>& energy)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    
    //Precomputation
    t_V edgesSq;
    igl::squared_edge_lengths(V, F, edgesSq);
    t_V edges = edgesSq.array().sqrt();
    t_V angles;
    igl::internal_angles_using_squared_edge_lengths(edgesSq, angles);
    t_Vv doubleAreas;
    igl::doublearea(edges, doubleAreas);
    
    energy = t_energy(V.rows());
    for(int vert=0; vert<V.rows(); ++vert) {
        if(isB[vert]) {
            energy(vert) = 0;
            continue;
        }
        
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        t_V_s anglesum = 0.;
        t_V_s A = 0.;
        for(int f=0; f<adjacentFaces.size(); ++f) {
            const int& face = adjacentFaces[f];
            const t_F_i& i = vert;
            const t_F_i& j = F(face, (adjacentFacesi[f]+1)%3);
            const t_F_i& k = F(face, (adjacentFacesi[f]+2)%3);
            const t_V_s& theta = angles(face, adjacentFacesi[f]);
            
            //Add angle sum and barycentric area
            anglesum += theta;
            A += doubleAreas(face)/6.;
        }
        
        const t_V_s pi = 2*asin(1.);
        t_V_s K = 2.*pi - anglesum;
        
        energy(vert) = K*K/A;
        
    }
    
}

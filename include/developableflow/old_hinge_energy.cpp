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


#include "old_hinge_energy.h"

#include <tools/kopp.h>
#include <tools/triangle_dN.h>
#include <tools/triangle_dTheta.h>

#include <igl/internal_angles.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/squared_edge_lengths.h>
#include <igl/doublearea.h>

#include <vector>
#include <list>



template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void old_hinge_energy_and_grad(
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
    t_V angles(V.rows(), 3);
    igl::internal_angles(V, F, angles);
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
    t_Vv dAreas;
    igl::doublearea(V, F, dAreas);
    t_V vertexNormalsRaw = t_V::Zero(V.rows(), 3);
    t_V vertexNormals = t_V(V.rows(), 3);
    for(int v=0; v<V.rows(); ++v) {
        for(const int& f : VF[v])
            vertexNormalsRaw.row(v) += dAreas(f)*faceNormals.row(f);
        vertexNormals.row(v) = vertexNormalsRaw.row(v).normalized();
    }
    
    
    //Actually compute energy
    energy = t_energy(V.rows());
    energyGrad = t_energyGrad::Zero(V.rows(), 3);
    for(int vert=0; vert<V.rows(); ++vert) {
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || VF[vert].size()<4) {
#else
            if(isB[vert]) {
#endif
                energy(vert) = 0;
                continue;
            }
            
            t_V33 mat = t_V33::Zero();
            
            const t_V3& Nv = vertexNormals.row(vert);
            const std::vector<indexType>& adjacentFaces = VF[vert];
            const std::vector<indexType>& adjacentFacesi = VFi[vert];
            for(int f=0; f<adjacentFaces.size(); ++f) {
                const int& face = adjacentFaces[f];
                const int& j = adjacentFacesi[f];
                const t_V_s& theta = angles(face, j);
                const t_V3& normal = faceNormals.row(face);
                
                mat += theta*normal*normal.transpose();
            }
            
            //Do eigendecomposition
            t_V33 diag, vecs;
            eigendecomp(mat, diag, vecs);
            t_V3 x;
            //        if(fabs(vecs.col(0).dot(Nv)) > fabs(vecs.col(1).dot(Nv)) && fabs(vecs.col(0).dot(Nv)) > fabs(vecs.col(2).dot(Nv))) {
            //0 is the third dimension. Ignore it
            //            energy(vert) = diag(1,1);
            //            x = vecs.col(1);
            //        } else {
            energy(vert) = diag(0,0);
            x = vecs.col(0);
            //        }
            
            //Compute gradient
            for(int f=0; f<adjacentFaces.size(); ++f) {
                const int& face = adjacentFaces[f];
                const t_V_s& theta = angles(face, adjacentFacesi[f]);
                const t_V3& normal = faceNormals.row(face);
                const t_V_s xdotN = normal.dot(x);
                
                //Get the vertices nomenclature right
                const int& i = vert;
                const int& j = F(face, (adjacentFacesi[f]+1)%3);
                const int& k = F(face, (adjacentFacesi[f]+2)%3);
                
                //Compute edges
                t_V3 ejk = V.row(k) - V.row(j);
                t_V3 eki = V.row(i) - V.row(k);
                t_V3 eij = V.row(j) - V.row(i);
                
                //Double area
                t_V_s dA = eij.cross(eki).norm();
                
                //Normal derivatives
                t_V33 dNdi, dNdj, dNdk;
                triangle_dN(eij, ejk, eki, normal, dA, dNdi, dNdj, dNdk);
                
                //Angle derivatives
                t_V3t dThetadi, dThetadj, dThetadk;
                triangle_dTheta(eij, eki, normal, dThetadi, dThetadj, dThetadk);
                
                //Putting it all together
                energyGrad.row(i) += xdotN*xdotN*dThetadi + 2*theta*xdotN*x.transpose()*dNdi;
                energyGrad.row(j) += xdotN*xdotN*dThetadj + 2*theta*xdotN*x.transpose()*dNdj;
                energyGrad.row(k) += xdotN*xdotN*dThetadk + 2*theta*xdotN*x.transpose()*dNdk;
            }
        }
    }
    
    
    template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedMinCurvatureDirs>
    IGL_INLINE void old_hinge_energy(
                                     const Eigen::PlainObjectBase<derivedV>& V,
                                     const Eigen::PlainObjectBase<derivedF>& F,
                                     const std::vector<std::vector<indexType> >& VF,
                                     const std::vector<std::vector<indexType> >& VFi,
                                     const std::vector<bool>& isB,
                                     Eigen::PlainObjectBase<derivedEnergy>& energy,
                                     Eigen::PlainObjectBase<derivedMinCurvatureDirs>& minCurvatureDirs)
    {
        typedef typename derivedV::Scalar t_V_s;
        typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
        typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
        typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
        typedef typename derivedF::Scalar t_F_i;
        typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
        typedef typename derivedEnergy::Scalar t_energy_s;
        typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
        typedef typename derivedMinCurvatureDirs::Scalar t_mincurvaturedirs_s;
        typedef Eigen::Matrix<t_mincurvaturedirs_s, Eigen::Dynamic, 3> t_MinCurvatureDirs;
        
        //Precomputation
        t_V angles(V.rows(), 3);
        igl::internal_angles(V, F, angles);
        t_V faceNormals(F.rows(), 3);
        igl::per_face_normals(V, F, faceNormals);
        
        //Actually compute energy
        energy = t_energy(V.rows());
        minCurvatureDirs = t_MinCurvatureDirs(V.rows(), 3);
        for(int vert=0; vert<V.rows(); ++vert) {
#ifdef IGNORE_VALENCE_3
            if(isB[vert] || VF[vert].size()<4) {
#else
                if(isB[vert]) {
#endif
                    energy(vert) = 0;
                    continue;
                }
                
                t_V33 mat = t_V33::Zero();
                
                const std::vector<indexType>& adjacentFaces = VF[vert];
                const std::vector<indexType>& adjacentFacesi = VFi[vert];
                for(int f=0; f<adjacentFaces.size(); ++f) {
                    const int& face = adjacentFaces[f];
                    const t_V_s& theta = angles(face, adjacentFacesi[f]);
                    const t_V3& normal = faceNormals.row(face);
                    
                    mat += theta*normal*normal.transpose();
                }
                
                //Do eigendecomposition
                t_V33 D, V;
                eigendecomp(mat, D, V);
                energy(vert) = D(0,0);
                minCurvatureDirs.row(vert) = V.col(0);
            }
            
        }
        
        template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy>
        IGL_INLINE void old_hinge_energy(
                                     const Eigen::PlainObjectBase<derivedV>& V,
                                     const Eigen::PlainObjectBase<derivedF>& F,
                                     const std::vector<std::vector<indexType> >& VF,
                                     const std::vector<std::vector<indexType> >& VFi,
                                     const std::vector<bool>& isB,
                                     Eigen::PlainObjectBase<derivedEnergy>& energy)
        {
            Eigen::PlainObjectBase<derivedV> dummy;
            old_hinge_energy(V, F, VF, VFi, isB, energy, dummy);
        }

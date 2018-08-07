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


#include "hinge_energy.h"

#include <tools/kopp.h>
#include <tools/triangle_dN.h>
#include <tools/triangle_dTheta.h>
#include <tools/crossproduct_matrix.h>

#include <igl/internal_angles.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/squared_edge_lengths.h>
#include <igl/doublearea.h>
#include <igl/parallel_for.h>

#include <vector>
#include <list>
#include <thread>
#include <mutex>

#include <Eigen/Sparse>


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void hinge_energy_and_grad(
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
    typedef Eigen::SparseMatrix<t_V_s> t_Vs;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3t;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergy::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    
    const t_V33 Id = t_V33::Identity();
    const auto macos = [] (const t_V_s& x)->t_V_s {
        return x > 1. ? acos(1.) : (x < -1. ? acos(-1.) : acos(x));
    };
    
    
    //Precomputation
    t_V angles(V.rows(), 3);
    igl::internal_angles(V, F, angles);
    
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
    
    t_Vv doubleAreas(F.rows());
    igl::doublearea(V, F, doubleAreas);
    t_Vv vertDoubleAreaSum = t_Vv::Zero(V.rows());
    for(int v=0; v<V.rows(); ++v)
        for(const int& f : VF[v])
            vertDoubleAreaSum(v) += doubleAreas(f);
    
    t_V vertexNormalsRaw = t_V::Zero(V.rows(), 3);
    t_V vertexNormals = t_V(V.rows(), 3);
    for(int v=0; v<V.rows(); ++v) {
        for(const int& f : VF[v])
            vertexNormalsRaw.row(v) += doubleAreas(f)*faceNormals.row(f);
        vertexNormals.row(v) = vertexNormalsRaw.row(v).normalized();
    }
    
    //Actually compute energy
    energy = t_energy(V.rows());
    energyGrad = t_energyGrad::Zero(V.rows(), 3);
    
    //Parallel stuff
    std::vector<std::vector<std::pair<t_F_i, t_V3t> > > energyGradsToAccumulate;
    const auto add_energyGrad = [&energyGrad, &energyGradsToAccumulate] (const int& thread, const t_F_i& i, const t_V3t& gi, const t_F_i& j, const t_V3t& gj, const t_F_i& k, const t_V3t& gk) {
        energyGradsToAccumulate[thread].emplace_back(i, gi);
        energyGradsToAccumulate[thread].emplace_back(j, gj);
        energyGradsToAccumulate[thread].emplace_back(k, gk);
    };
    
    const auto prep_loop = [&energyGradsToAccumulate, &V] (const int& threadNum) {
        energyGradsToAccumulate.resize(threadNum);
    };
    
    const auto handle_vertex = [&] (const t_F_i& vert, const int& thread) {
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || VF[vert].size()<4) {
#else
        if(isB[vert]) {
#endif
            energy(vert) = 0;
            return;
        }
        
        const t_V3& Nv = vertexNormals.row(vert);
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        t_V33 mat = t_V33::Zero();
        for(int f=0; f<adjacentFaces.size(); ++f) {
            const t_V_s& thetaf = angles(adjacentFaces[f], adjacentFacesi[f]);
            const t_V3& Nf = faceNormals.row(adjacentFaces[f]);
            
            //if(Nv.cross(Nf).norm() > 1e-6) {
            if(1.+Nv.cross(Nf).norm()!=1.) {
                const t_V3 Nfw = Nv.cross(Nf).cross(Nv).normalized()*macos(Nv.dot(Nf));
                mat += thetaf*Nfw*Nfw.transpose();
                assert(mat==mat);
            }
        }
        
        //Do eigendecomposition
        t_V33 diag, vecs;
        eigendecomp(mat, diag, vecs);
        t_V3 x;
        if(fabs(vecs.col(0).dot(Nv)) > 0.1) {
            //0 is the third dimension. Ignore it
            energy(vert) = diag(1,1);
            x = vecs.col(1);
        } else {
            energy(vert) = diag(0,0);
            x = vecs.col(0);
        }
        assert(energy(vert)==energy(vert) && x==x);
        //if(1.+energy(vert)==1.)
        //    continue;
        
        
        //Compute gradient
        for(int f=0; f<adjacentFaces.size(); ++f) {
            const t_V_s& thetaf = angles(adjacentFaces[f], adjacentFacesi[f]);
            const t_V3& Nf = faceNormals.row(adjacentFaces[f]);
            const t_V_s sinphif = Nv.cross(Nf).norm();
            const t_V_s phif = macos(Nv.dot(Nf));
            const t_V_s tanphif = sinphif/Nv.dot(Nf);
            const t_V3 nuf = Nv.cross(Nf).normalized();
            const t_V3 muvf = nuf.cross(Nv);
            const t_V3 muff = nuf.cross(Nf);
            const t_V3 Nfw = muvf*phif;
            
            //Get the vertices nomenclature right
            const int& i = vert;
            const int& jf = F(adjacentFaces[f], (adjacentFacesi[f]+1)%3);
            const int& kf = F(adjacentFaces[f], (adjacentFacesi[f]+2)%3);
            
            //Compute edges
            const t_V3 ejkf = V.row(kf) - V.row(jf);
            const t_V3 ekif = V.row(i) - V.row(kf);
            const t_V3 eijf = V.row(jf) - V.row(i);
            
            //Double area
            const t_V_s& dAf = doubleAreas(adjacentFaces[f]);
            
            //Angle derivatives
            t_V3t dThetafdi, dThetafdj, dThetafdk;
            triangle_dTheta(eijf, ekif, Nf, dThetafdi, dThetafdj, dThetafdk);
            
            //Normal derivatives
            t_V33 dNfdi, dNfdj, dNfdk;
            triangle_dN(eijf, ejkf, ekif, Nf, dAf, dNfdi, dNfdj, dNfdk);
            
            if(1.+sinphif!=1.) {
                //Geometric grad in Nfw
                const t_V_s xdotNfw = Nfw.dot(x);
                const t_V3t factorf = 2*xdotNfw*thetaf*x.transpose()*(muvf*muff.transpose() + phif/sinphif*nuf*nuf.transpose());
                const t_V3t factorv = -2*xdotNfw*thetaf*x.transpose()*((muvf+phif*Nv)*muvf.transpose() + phif/tanphif*nuf*nuf.transpose())/vertexNormalsRaw.row(vert).norm(); //The id-Nv*Nv.transpose() is actually irrelevant, since it just maps to the vectors muv, nu anyways
                assert(factorf==factorf && factorv==factorv);
                
                add_energyGrad(thread, i, xdotNfw*xdotNfw*dThetafdi + factorf*dNfdi, jf, xdotNfw*xdotNfw*dThetafdj + factorf*dNfdj, kf, xdotNfw*xdotNfw*dThetafdk + factorf*dNfdk);
                
                for(int g=0; g<adjacentFaces.size(); ++g) {
                    const int& face = adjacentFaces[g];
                    const int& i = vert;
                    const int& j = F(face, (adjacentFacesi[g]+1)%3);
                    const int& k = F(face, (adjacentFacesi[g]+2)%3);
                    
                    t_V33 Jeij, Jeik;
                    crossproduct_matrix((V.row(j)-V.row(i)).eval(), Jeij);
                    crossproduct_matrix((V.row(k)-V.row(i)).eval(), Jeik);
                    const t_V33 dPdi = -Jeij+Jeik, dPdj = -Jeik, dPdk = Jeij;
                    
                    add_energyGrad(thread, i, factorv*dPdi, j, factorv*dPdj, k, factorv*dPdk);
                }
            }
            
        }
    };
    
    const auto accumulate = [&energyGrad, &energyGradsToAccumulate] (const int& thread) {
        for(const std::pair<t_F_i, t_V3t>& i : energyGradsToAccumulate[thread])
            energyGrad.row(i.first) += i.second;
    };
    
#ifndef PARALLEL_COMPUTATION
    //SERIAL VERSION
    prep_loop(1);
    for(int vert=0; vert<V.rows(); ++vert)
        handle_vertex(vert, 0);
    accumulate(0);
#else
    //PARALLEL VERSION
    igl::parallel_for(V.rows(), prep_loop, handle_vertex, accumulate);
#endif
    
}


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedMinCurvatureDirs>
IGL_INLINE void hinge_energy(
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
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedMinCurvatureDirs::Scalar t_mincurvaturedirs_s;
    typedef Eigen::Matrix<t_mincurvaturedirs_s, Eigen::Dynamic, 3> t_MinCurvatureDirs;
    
    const auto macos = [] (const t_V_s& x)->t_V_s {
        return x > 1. ? acos(1.) : (x < -1. ? acos(-1.) : acos(x));
    };
    
    //Precomputation
    t_V angles(V.rows(), 3);
    igl::internal_angles(V, F, angles);
    
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
    
    t_Vv doubleAreas(F.rows());
    igl::doublearea(V, F, doubleAreas);
    t_Vv vertDoubleAreaSum = t_Vv::Zero(V.rows());
    for(int v=0; v<V.rows(); ++v)
        for(const int& f : VF[v])
            vertDoubleAreaSum(v) += doubleAreas(f);
    
    t_V vertexNormalsRaw = t_V::Zero(V.rows(), 3);
    t_V vertexNormals = t_V(V.rows(), 3);
    for(int v=0; v<V.rows(); ++v) {
        for(const int& f : VF[v])
            vertexNormalsRaw.row(v) += doubleAreas(f)*faceNormals.row(f);
        vertexNormals.row(v) = vertexNormalsRaw.row(v).normalized();
    }
    
    
    //Actually compute energy
    energy = t_energy(V.rows());
    minCurvatureDirs = t_MinCurvatureDirs(V.rows(), 3);
    const auto handle_vertex = [&] (const t_F_i& vert) {
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || VF[vert].size()<4) {
#else
            if(isB[vert]) {
#endif
                energy(vert) = 0;
                minCurvatureDirs.row(vert).setZero();
                return;
            }
            
            const t_V3& Nv = vertexNormals.row(vert);
            const std::vector<indexType>& adjacentFaces = VF[vert];
            const std::vector<indexType>& adjacentFacesi = VFi[vert];
            
            t_V33 mat = t_V33::Zero();
            for(int f=0; f<adjacentFaces.size(); ++f) {
                const t_V_s& thetaf = angles(adjacentFaces[f], adjacentFacesi[f]);
                const t_V3& Nf = faceNormals.row(adjacentFaces[f]);
                
                //if(Nv.cross(Nf).norm() > 1e-6) {
                if(1.+Nv.cross(Nf).norm()!=1.) {
                    const t_V3 Nfw = Nv.cross(Nf).cross(Nv).normalized()*macos(Nv.dot(Nf));
                    mat += thetaf*Nfw*Nfw.transpose();
                    assert(mat==mat);
                }
            }
            
            //Do eigendecomposition
            t_V33 diag, vecs;
            eigendecomp(mat, diag, vecs);
            if(fabs(vecs.col(0).dot(Nv)) > 0.1) {
                //0 is the third dimension. Ignore it
                energy(vert) = diag(1,1);
                minCurvatureDirs.row(vert) = vecs.col(1);
            } else {
                energy(vert) = diag(0,0);
                minCurvatureDirs.row(vert) = vecs.col(0);
            }
            assert(energy(vert)==energy(vert) && minCurvatureDirs.row(vert)==minCurvatureDirs.row(vert));
        };
        
#ifndef PARALLEL_COMPUTATION
        //SERIAL VERSION
        for(int vert=0; vert<V.rows(); ++vert)
            handle_vertex(vert);
#else
        //PARALLEL VERSION
        igl::parallel_for(V.rows(), handle_vertex);
#endif
        
    }
    
    
    template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy>
    IGL_INLINE void hinge_energy(
                                 const Eigen::PlainObjectBase<derivedV>& V,
                                 const Eigen::PlainObjectBase<derivedF>& F,
                                 const std::vector<std::vector<indexType> >& VF,
                                 const std::vector<std::vector<indexType> >& VFi,
                                 const std::vector<bool>& isB,
                                 Eigen::PlainObjectBase<derivedEnergy>& energy)
    {
        Eigen::PlainObjectBase<derivedV> dummy;
        hinge_energy(V, F, VF, VFi, isB, energy, dummy);
    }

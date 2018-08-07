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


#include "max_hinge_energy.h"

#include <tools/kopp.h>

#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/per_face_normals.h>
#include <igl/squared_edge_lengths.h>

#include <set>

#define INFTY std::numeric_limits<double>::infinity()

#undef HINGE_SAMPLE_MIDARCS

template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void max_hinge_energy_and_grad(
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
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 2> t_E;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 1> t_Fv;
    typedef Eigen::Matrix<t_F_i, 2, 1> t_F2;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> t_Ebool;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergy::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    
    
    const t_V33 Id = t_V33::Identity();
    const auto macos = [] (const t_V_s& x)->t_V_s {
        return x > 1. ? acos(1.) : (x < -1. ? acos(-1.) : acos(x));
    };
    
    //Precomputation
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
    
    //Energy
    energy = t_energy(V.rows());
    t_E normalIndices(V.rows(), 2);
    t_E normalSigns(V.rows(), 2);
    t_E maxIndices(V.rows(), 2);
    t_V us(V.rows(), 3);
    const auto handle_vertex_energy = [&] (const t_F_i& vert) {
        //Boundary vertices contribute no energy
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
        
        energy(vert) = INFTY;
        
        //Loop through all pairs
        const auto handle_pair = [&](const int& i, const int& j, const t_V_s& signi, const t_V_s& signj) {
            const t_V3& Nfi = faceNormals.row(adjacentFaces[i]), &Nfj = faceNormals.row(adjacentFaces[j]);
            //if(1.+Nv.cross(Nfi).norm()==1. || 1.+Nv.cross(Nfj).norm()==1.)
            if(Nv.cross(Nfi).norm() < 1e-6 || Nv.cross(Nfj).norm() < 1e-6)
                return;
            const t_V3 Nfwi = signi*Nv.cross(Nfi).cross(Nv).normalized()*macos(Nv.dot(Nfi));
            const t_V3 Nfwj = signj*Nv.cross(Nfj).cross(Nv).normalized()*macos(Nv.dot(Nfj));
            const t_V3 base = Nfwi-Nfwj;
            //if(1.+base.norm()==1.)
            if(base.norm() < 1e-6)
                return;
            const t_V3 u = Nv.cross(base.normalized());
            
            t_F2 localMaxface = t_F2(-1,-1);
            t_V_s localEnergy = -1;
            for(int k1=0; k1<adjacentFaces.size(); ++k1) {
                const t_V3& Nfk1 = faceNormals.row(adjacentFaces[k1]);
                const t_V3 Nfwk1 = 1.+Nv.cross(Nfk1).norm()==1. ? t_V3::Zero() : (Nv.cross(Nfk1).cross(Nv).normalized()*macos(Nv.dot(Nfk1))).eval();
                for(int k2=k1; k2<adjacentFaces.size(); ++k2) {
                    const t_V3& Nfk2 = faceNormals.row(adjacentFaces[k2]);
                    const t_V3 Nfwk2 = 1.+Nv.cross(Nfk2).norm()==1. ? t_V3::Zero() : (Nv.cross(Nfk2).cross(Nv).normalized()*macos(Nv.dot(Nfk2))).eval();
                    const t_V_s udN = u.dot(Nfwk1-Nfwk2);
                    if(localEnergy < udN*udN) {
                        localEnergy = udN*udN;
                        localMaxface << k1, k2;
                    }
                }
            }
            
            if(localEnergy >= 0 && localEnergy < energy(vert)) {
                energy(vert) = localEnergy;
                normalIndices.row(vert) << i, j;
                normalSigns.row(vert) << signi, signj;
                maxIndices.row(vert) = localMaxface;
                us.row(vert) = u;
            }
        };
        
        for(int i=0; i<adjacentFaces.size(); ++i) {
            for(int j=i+1; j<adjacentFaces.size(); ++j) {
                handle_pair(i, j, 1, 1);
                handle_pair(i, j, 1, -1);
            }
            handle_pair(i, i, 1, -1);
        }
        
        if(energy(vert) == INFTY)
            energy(vert) = 0;
    };
    
#ifndef PARALLEL_COMPUTATION
    //SERIAL VERSION
    for(int vert=0; vert<V.rows(); ++vert)
        handle_vertex_energy(vert);
#else
    //PARALLEL VERSION
    igl::parallel_for(V.rows(), handle_vertex_energy);
#endif
    
    
    //Gradient
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
    
    const auto handle_vertex_grad = [&] (const t_F_i& vert, const int& thread) {
        //Boundary vertices contribute no energy
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || 1.+energy(vert)==1. || VF[vert].size()<4) {
#else
        if(isB[vert] || 1.+energy(vert)==1.) {
#endif
            energy(vert) = 0;
            return;
        }
        
        const t_V3& Nv = vertexNormals.row(vert);
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        const int& facei = adjacentFaces[normalIndices(vert,0)], &facej = adjacentFaces[normalIndices(vert,1)], &facek1 = adjacentFaces[maxIndices(vert,0)], &facek2 = adjacentFaces[maxIndices(vert,1)];
        const t_V3& Nfi = faceNormals.row(facei), &Nfj = faceNormals.row(facej), &Nfk1 = faceNormals.row(facek1), &Nfk2 = faceNormals.row(facek2);
        const t_V_s sinphii = Nv.cross(Nfi).norm(), sinphij = Nv.cross(Nfj).norm(), sinphik1 = Nv.cross(Nfk1).norm(), sinphik2 = Nv.cross(Nfk2).norm();
        const t_V_s phii = macos(Nv.dot(Nfi)), phij = macos(Nv.dot(Nfj)), phik1 = macos(Nv.dot(Nfk1)), phik2 = macos(Nv.dot(Nfk2));
        const t_V_s tanphii = sinphii/Nv.dot(Nfi), tanphij = sinphij/Nv.dot(Nfj), tanphik1 = sinphik1/Nv.dot(Nfk1), tanphik2 = sinphik2/Nv.dot(Nfk2);
        const t_V3 nui = Nv.cross(Nfi).normalized(), nuj = Nv.cross(Nfj).normalized(), nuk1 = Nv.cross(Nfk1).normalized(), nuk2 = Nv.cross(Nfk2).normalized();
        const t_V3 muvi = nui.cross(Nv), muvj = nuj.cross(Nv), muvk1 = nuk1.cross(Nv), muvk2 = nuk2.cross(Nv);
        const t_V3 mufi = nui.cross(Nfi), mufj = nuj.cross(Nfj), mufk1 = nuk1.cross(Nfk1), mufk2 = nuk2.cross(Nfk2);
        const t_V3 Nfwi = 1.+sinphii==1. ? t_V3::Zero() : (normalSigns(vert,0)*muvi*phii).eval(), Nfwj = 1.+sinphij==1. ? t_V3::Zero() : (normalSigns(vert,1)*muvj*phij).eval(), Nfwk1 = 1.+sinphik1==1. ? t_V3::Zero() : (muvk1*phik1).eval(), Nfwk2 = 1.+sinphik2==1. ? t_V3::Zero() : (muvk2*phik2).eval();
        const t_V3& u = us.row(vert);
        
        //Get the vertices nomenclature right
        const int& ii = vert, &ij = vert, &ik1 = vert, &ik2 = vert;
        const int& ji = F(facei, (adjacentFacesi[normalIndices(vert,0)]+1)%3), &jj = F(facej, (adjacentFacesi[normalIndices(vert,1)]+1)%3), &jk1 = F(facek1, (adjacentFacesi[maxIndices(vert,0)]+1)%3), &jk2 = F(facek2, (adjacentFacesi[maxIndices(vert,1)]+1)%3);
        const int& ki = F(facei, (adjacentFacesi[normalIndices(vert,0)]+2)%3), &kj = F(facej, (adjacentFacesi[normalIndices(vert,1)]+2)%3), &kk1 = F(facek1, (adjacentFacesi[maxIndices(vert,0)]+2)%3), &kk2 = F(facek2, (adjacentFacesi[maxIndices(vert,1)]+2)%3);
        
        //Compute edges
        t_V3 ejki = V.row(ki) - V.row(ji), ejkj = V.row(kj) - V.row(jj), ejkk1 = V.row(kk1) - V.row(jk1), ejkk2 = V.row(kk2) - V.row(jk2);
        t_V3 ekii = V.row(ii) - V.row(ki), ekij = V.row(ij) - V.row(kj), ekik1 = V.row(ik1) - V.row(kk1), ekik2 = V.row(ik2) - V.row(kk2);
        t_V3 eiji = V.row(ji) - V.row(ii), eijj = V.row(jj) - V.row(ij), eijk1 = V.row(jk1) - V.row(ik1), eijk2 = V.row(jk2) - V.row(ik2);
        
        //Double area
        const t_V_s& dAi = doubleAreas(facei), &dAj = doubleAreas(facej), &dAk1 = doubleAreas(facek1), &dAk2 = doubleAreas(facek2);
        
        //Normal derivatives
        t_V33 dNdii, dNdji, dNdki, dNdij, dNdjj, dNdkj, dNdik1, dNdjk1, dNdkk1, dNdik2, dNdjk2, dNdkk2;
        triangle_dN(eiji, ejki, ekii, Nfi, dAi, dNdii, dNdji, dNdki);
        triangle_dN(eijj, ejkj, ekij, Nfj, dAj, dNdij, dNdjj, dNdkj);
        triangle_dN(eijk1, ejkk1, ekik1, Nfk1, dAk1, dNdik1, dNdjk1, dNdkk1);
        triangle_dN(eijk2, ejkk2, ekik2, Nfk2, dAk2, dNdik2, dNdjk2, dNdkk2);
        
        //Altitude gradients
        t_V33 Jnv, Jnimnj;
        crossproduct_matrix(Nv, Jnv);
        crossproduct_matrix((Nfwi-Nfwj).eval(), Jnimnj);
        const t_V33 dudi = (Id-u*u.transpose())/(Nfwi-Nfwj).norm()*normalSigns(vert,0)*Jnv;
        const t_V33 dudj = -(Id-u*u.transpose())/(Nfwi-Nfwj).norm()*normalSigns(vert,1)*Jnv;
        const t_V33 dudv = -(Id-u*u.transpose())/(Nfwi-Nfwj).norm()*Jnimnj;
        assert(dudi==dudi && dudj==dudj && dudv==dudv);
        
        //Geometric grad
        const t_V3 Nfwk = Nfwk1-Nfwk2;
        const t_V_s udotNfwk = u.dot(Nfwk);
        const t_V3t factorfi = 2*udotNfwk*Nfwk.transpose()*dudi*(muvi*mufi.transpose() + phii/sinphii*nui*nui.transpose());
        const t_V3t factorfj = 2*udotNfwk*Nfwk.transpose()*dudj*(muvj*mufj.transpose() + phij/sinphij*nuj*nuj.transpose());
        const t_V3t factorfk1 = 2*udotNfwk*u.transpose()*(muvk1*mufk1.transpose() + phik1/sinphik1*nuk1*nuk1.transpose());
        const t_V3t factorfk2 = -2*udotNfwk*u.transpose()*(muvk2*mufk2.transpose() + phik2/sinphik2*nuk2*nuk2.transpose());
        const t_V3t factorvi = -2*udotNfwk*Nfwk.transpose()*dudi*((muvi+phii*Nv)*muvi.transpose() + phii/tanphii*nui*nui.transpose())/vertexNormalsRaw.row(vert).norm();
        const t_V3t factorvj = -2*udotNfwk*Nfwk.transpose()*dudj*((muvj+phij*Nv)*muvj.transpose() + phij/tanphij*nuj*nuj.transpose())/vertexNormalsRaw.row(vert).norm();
        const t_V3t factorvk1 = -2*udotNfwk*u.transpose()*((muvk1+phik1*Nv)*muvk1.transpose() + phik1/tanphik1*nuk1*nuk1.transpose())/vertexNormalsRaw.row(vert).norm();
        const t_V3t factorvk2 = 2*udotNfwk*u.transpose()*((muvk2+phik2*Nv)*muvk2.transpose() + phik2/tanphik2*nuk2*nuk2.transpose())/vertexNormalsRaw.row(vert).norm();
        const t_V3t factorv = 2*udotNfwk*Nfwk.transpose()*dudv/vertexNormalsRaw.row(vert).norm();
        assert(factorfi==factorfi && factorfj==factorfj && factorfk1==factorfk1 && factorfk2==factorfk2 && factorvi==factorvi && factorvj==factorvj && factorvj==factorvj && factorvk1==factorvk1 && factorvk2==factorvk2 && factorv==factorv);
        
        add_energyGrad(thread, ii, factorfi*dNdii, ji, factorfi*dNdji, ki, factorfi*dNdki);
        add_energyGrad(thread, ij, factorfj*dNdij, jj, factorfj*dNdjj, kj, factorfj*dNdkj);
        add_energyGrad(thread, ik1, factorfk1*dNdik1, jk1, factorfk1*dNdjk1, kk1, factorfk1*dNdkk1);
        add_energyGrad(thread, ik2, factorfk2*dNdik2, jk2, factorfk2*dNdjk2, kk2, factorfk2*dNdkk2);
        
        for(int g=0; g<adjacentFaces.size(); ++g) {
            const int& face = adjacentFaces[g];
            const int& i = vert;
            const int& j = F(face, (adjacentFacesi[g]+1)%3);
            const int& k = F(face, (adjacentFacesi[g]+2)%3);
            
            t_V33 Jeij, Jeik;
            crossproduct_matrix((V.row(j)-V.row(i)).eval(), Jeij);
            crossproduct_matrix((V.row(k)-V.row(i)).eval(), Jeik);
            const t_V33 dPdi = -Jeij+Jeik, dPdj = -Jeik, dPdk = Jeij;
            
            const t_V3t factors = factorvi+factorvj+factorvk1+factorvk2+factorv;
            add_energyGrad(thread, i, factors*dPdi, j, factors*dPdj, k, factors*dPdk);
            assert(energyGrad.row(i)==energyGrad.row(i) && energyGrad.row(j)==energyGrad.row(j) && energyGrad.row(k)==energyGrad.row(k));
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
        handle_vertex_grad(vert, 0);
    accumulate(0);
#else
    //PARALLEL VERSION
    igl::parallel_for(V.rows(), prep_loop, handle_vertex_grad, accumulate);
#endif
    
}


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy>
IGL_INLINE void max_hinge_energy(
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
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3t;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 2> t_E;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 1> t_Fv;
    typedef Eigen::Matrix<t_F_i, 2, 1> t_F2;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergy::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    
    const auto macos = [] (const t_V_s& x)->t_V_s {
        return x > 1. ? acos(1.) : (x < -1. ? acos(-1.) : acos(x));
    };
    
    //Precomputation
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
    
    //Energy
    energy = t_energy(V.rows());
    t_E normalIndices(V.rows(), 2);
    t_E normalSigns(V.rows(), 2);
    t_E maxIndices(V.rows(), 2);
    t_V us(V.rows(), 3);
    const auto handle_vertex_energy = [&] (const t_F_i& vert) {
        //Boundary vertices contribute no energy
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || 1.+energy(vert)==1. || VF[vert].size()<4) {
#else
        if(isB[vert] || 1.+energy(vert)==1.) {
#endif
            energy(vert) = 0;
            return;
        }
        
        const t_V3& Nv = vertexNormals.row(vert);
        const std::vector<indexType>& adjacentFaces = VF[vert];
        
        energy(vert) = INFTY;
        
        //Loop through all pairs
        const auto handle_pair = [&](const int& i, const int& j, const t_V_s& signi, const t_V_s& signj) {
            const t_V3& Nfi = faceNormals.row(adjacentFaces[i]), &Nfj = faceNormals.row(adjacentFaces[j]);
            //if(1.+Nv.cross(Nfi).norm()==1. || 1.+Nv.cross(Nfj).norm()==1.)
            if(Nv.cross(Nfi).norm() < 1e-6 || Nv.cross(Nfj).norm() < 1e-6)
                return;
            const t_V3 Nfwi = signi*Nv.cross(Nfi).cross(Nv).normalized()*macos(Nv.cross(Nfi));
            const t_V3 Nfwj = signj*Nv.cross(Nfj).cross(Nv).normalized()*macos(Nv.cross(Nfj));
            const t_V3 base = Nfwi-Nfwj;
            //if(1.+base.norm()==1.)
            if(base.norm() < 1e-6)
                return;
            const t_V3 u = Nv.cross(base.normalized());
            
            t_F2 localMaxface = t_F2(-1,-1);
            t_V_s localEnergy = -1;
            for(int k1=0; k1<adjacentFaces.size(); ++k1) {
                const t_V3& Nfk1 = faceNormals.row(adjacentFaces[k1]);
                const t_V3 Nfwk1 = 1.+Nv.cross(Nfk1).norm()==1. ? t_V3::Zero() : (Nv.cross(Nfk1).cross(Nv).normalized()*macos(Nv.cross(Nfk1))).eval();
                for(int k2=k1; k2<adjacentFaces.size(); ++k2) {
                    const t_V3& Nfk2 = faceNormals.row(adjacentFaces[k2]);
                    const t_V3 Nfwk2 = 1.+Nv.cross(Nfk2).norm()==1. ? t_V3::Zero() : (Nv.cross(Nfk2).cross(Nv).normalized()*macos(Nv.cross(Nfk2))).eval();
                    const t_V_s udN = u.dot(Nfwk1-Nfwk2);
                    if(localEnergy < udN*udN) {
                        localEnergy = udN*udN;
                        localMaxface << k1, k2;
                    }
                }
            }
            
            if(localEnergy >= 0 && localEnergy < energy(vert)) {
                energy(vert) = localEnergy;
                normalIndices.row(vert) << i, j;
                normalSigns.row(vert) << signi, signj;
                maxIndices.row(vert) = localMaxface;
                us.row(vert) = u;
            }
        };
        
        for(int i=0; i<adjacentFaces.size(); ++i) {
            for(int j=i+1; j<adjacentFaces.size(); ++j) {
                handle_pair(i, j, 1, 1);
                handle_pair(i, j, 1, -1);
            }
            handle_pair(i, i, 1, -1);
        }
        
        if(energy(vert) == INFTY)
            energy(vert) = 0;
    };
    
#ifndef PARALLEL_COMPUTATION
    //SERIAL VERSION
    for(int vert=0; vert<V.rows(); ++vert)
        handle_vertex_energy(vert);
#else
    //PARALLEL VERSION
    igl::parallel_for(V.rows(), handle_vertex_energy);
#endif
    
}

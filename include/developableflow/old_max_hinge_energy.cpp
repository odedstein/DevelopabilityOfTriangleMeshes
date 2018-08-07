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


#include "old_max_hinge_energy.h"

#include <tools/kopp.h>

#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/per_face_normals.h>
#include <igl/squared_edge_lengths.h>

#include <set>

#define INFTY std::numeric_limits<double>::infinity()

template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void old_max_hinge_energy_and_grad(
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
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vvec;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3t;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 1> t_Fvec;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergy::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    
    //Helper function
    const auto orthogonalVector = [] (const t_V3& u) {
        t_V3 un = u.normalized();
        t_V3 ret;
        if(un[2]-1. == -2.) {
            ret << 0., -1., 0.;
        } else {
            const t_V_s a = 1./(1. + un[2]);
            const t_V_s b = -un[0]*un[1]*a;
            ret << 1. - un[0]*un[0]*a, b, -un[0];
        }
        return ret;
    };
    
    //Precomputation
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
    t_Vvec doubleAreas;
    igl::doublearea(doublecast(V), F, doubleAreas);
    
    //Actually compute energy
    energy = t_energy(V.rows());
    t_V umw = t_V::Zero(V.rows(), 3);
    t_V spmw(V.rows(), 3);
    t_F fpmw(V.rows(), 3);
    t_Fvec maxFaces(V.rows());
    for(int vert=0; vert<V.rows(); ++vert) {
        //Boundary vertices contribute no energy
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || VF[vert].size()<4)
#else
        if(isB[vert])
#endif
        {
            energy(vert) = 0;
            umw.row(vert) << 0., 0., 0.;
            fpmw.row(vert) << -1, -1, -1;
            spmw.row(vert) << 0., 0., 0.;
            maxFaces(vert) = -1;
            continue;
        }
        
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        //Do we have a small number of adjacent triangles?
        if(adjacentFaces.size() == 1) {
            umw.row(vert) << (V.row(F(adjacentFaces[0], (adjacentFacesi[0]+1)%3)) - V.row(vert)).normalized();
            fpmw.row(vert) << adjacentFaces[0], adjacentFaces[0], adjacentFaces[0];
            spmw.row(vert) << 1., 1., 1.;
            maxFaces(vert) = adjacentFaces[0];
            const t_V_s udotn = umw.row(vert).dot(faceNormals.row(maxFaces(vert)));
            energy(vert) = 0.5*udotn*udotn;
            continue;
        } else if(adjacentFaces.size() == 2) {
            umw.row(vert) << faceNormals.row(adjacentFaces[0]).cross(faceNormals.row(adjacentFaces[1])).normalized();
            fpmw.row(vert) << adjacentFaces[0], adjacentFaces[1], adjacentFaces[1];
            spmw.row(vert) << 1., 1., 1.;
            maxFaces(vert) = adjacentFaces[0];
            const t_V_s udotn = umw.row(vert).dot(faceNormals.row(maxFaces(vert)));
            energy(vert) = 0.5*udotn*udotn;
            continue;
        }
        
        //Loop through all face triplets
        t_energy_s& minEnergy = energy(vert);
        minEnergy = INFTY;
        for(int i=0; i<adjacentFaces.size(); ++i) {
            const t_F_i& facei = adjacentFaces[i];
            const t_V3& Ni = faceNormals.row(facei);
            for(int j=0; j<i; ++j) {
                const t_F_i& facej = adjacentFaces[j];
                const t_V3& Nj = faceNormals.row(facej);
                for(int k=0; k<j; ++k) {
                    const t_F_i& facek = adjacentFaces[k];
                    const t_V3& Nk = faceNormals.row(facek);
                    
                    t_V3 u;
                    bool uComputed = false;
                    
                    //Do we have a degenerate triangle configuration?
                    const t_V_s dotij = Ni.dot(Nj);
                    const t_V_s dotjk = Nj.dot(Nk);
                    const t_V_s dotki = Nk.dot(Ni);
                    if(dotij > 1.-1e-7 && dotjk > 1e-7) {
                    //if(1.+dotij >= 2. && 1.+dotjk >= 2.) {
                        u = orthogonalVector(Ni);
                        assert(u==u);
                        uComputed = true;
                        //std::cout << "total degeneracy. Returned norm vec " << u.transpose() << " to vec " << Ni.transpose() << ", dp: " << u.dot(Ni) << std::endl;
                    } else if(dotij > 1.-1e-7) {
                    //} else if(1.+dotij >= 2.) {
                        u = Ni.cross(Nk).normalized();
                        assert(u==u);
                        uComputed = true;
                    } else if(dotjk > 1.-1e-7) {
                    //} else if(1.+dotjk >= 2.) {
                        u = Nj.cross(Ni).normalized();
                        assert(u==u);
                        uComputed = true;
                    } else if(dotki > 1.-1e-7) {
                    //} else if(1.+dotki >= 2.) {
                        u = Nk.cross(Nj).normalized();
                        assert(u==u);
                        uComputed = true;
                    }
                    
                    //Loop through all sign combinations
                    for(t_V_s p : {/*-1.,*/ 1.}) {
                        for(t_V_s q : {-1., 1.}) {
                            for(t_V_s r : {-1., 1.}) {
                                if(!uComputed) {
                                    u = (q*Nj-p*Ni).cross(r*Nk-p*Ni);
                                    const t_V_s unorm = u.norm();
                                    if(u.norm() < 1e-7)
                                        continue;
                                    u /= unorm;
                                    assert(u==u);
                                }
                                
                                //Loop through adj faces to find maximum
                                t_energy_s Eijk = 0.;
                                t_F_i maxFace = adjacentFaces[0];
                                for(int f=0; f<adjacentFaces.size(); ++f) {
                                    const t_F_i& face = adjacentFaces[f];
                                    const t_V3& N = faceNormals.row(face);
                                    
                                    const t_V_s udotn = u.dot(N);
                                    const t_V_s localE = 0.5*udotn*udotn;
                                    if(localE > Eijk) {
                                        Eijk = localE;
                                        maxFace = face;
                                    }
                                }
                                
                                //Are we a per-vertex minimum?
                                if(Eijk < minEnergy) {
                                    minEnergy = Eijk;
                                    umw.row(vert) = u;
                                    fpmw.row(vert) << facei, facej, facek;
                                    spmw.row(vert) << p, q, r;
                                    maxFaces(vert) = maxFace;
                                }
                            }
                        }
                    }
                    
                }
                
            }
        }
        assert(umw.row(vert) == umw.row(vert));
        assert(umw.row(vert).norm() > 0);
        assert(maxFaces(vert) >= 0);
    }
    
    //Compute gradients
    energyGrad = t_energyGrad::Zero(V.rows(), 3);
    for(int vert=0; vert<V.rows(); ++vert) {
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || VF[vert].size()<4)
#else
        if(isB[vert])
#endif
            continue;
        
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        const t_V3& u = umw.row(vert);
        const t_F_i& maxFace = maxFaces(vert);
        const t_V3& N = faceNormals.row(maxFace);
        const t_V_s& dA = doubleAreas(maxFace);
        
        //Contributions from max face
        const t_V_s factor1 = u.dot(N);
        for(int j=0; j<3; ++j) {
            const t_V3 edge = V.row(F(maxFace, (j+2)%3)) - V.row(F(maxFace, (j+1)%3));
            energyGrad.row(F(maxFace, j)) += factor1/dA*u.dot(edge.cross(N))*N;
            assert(factor1/dA*u.dot(edge.cross(N))*N == factor1/dA*u.dot(edge.cross(N))*N);
        }
        
        t_V3 Np[3], ep[3];
        t_V_s dAp[3];
        for(int p=0; p<3; ++p) {
            Np[p] = spmw(vert,p)*faceNormals.row(fpmw(vert,p));
            dAp[p] = doubleAreas(fpmw(vert,p));
        }
        for(int p=0; p<3; ++p) {
            ep[p] = Np[(p+2)%3] - Np[(p+1)%3];
        }
        const t_V_s a2 = ep[0].squaredNorm(), b2 = ep[1].squaredNorm(), c2 = ep[2].squaredNorm();
        const t_V_s dAabc = 0.5*sqrt(4.*a2*b2 - (a2+b2-c2)*(a2+b2-c2));
        if(std::isnan(dAabc) || dAabc < 1e-3)
        //if(std::isnan(dAabc) || 1.+dAabc == 1.)
            continue;
        
        //Contributions from max face triplet
        for(int p=0; p<3; ++p) {
            const t_F_i& face = fpmw(vert,p);
            const t_V_s factor2 = factor1*ep[p].cross(u).dot(N)/dAp[p]/dAabc;
            for(int j=0; j<3; ++j) {
                const t_V3 edge = ( V.row(F(face, (j+2)%3)) - V.row(F(face, (j+1)%3)) ) * spmw(vert,p);
                energyGrad.row(F(face, j)) += factor2*u.dot(edge.cross(Np[p]))*Np[p];
                assert(factor2*u.dot(edge.cross(Np[p]))*Np[p] == factor2*u.dot(edge.cross(Np[p]))*Np[p]);
            }
         }
    }
    
}


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy>
IGL_INLINE void old_max_hinge_energy(
                                 const Eigen::PlainObjectBase<derivedV>& V,
                                 const Eigen::PlainObjectBase<derivedF>& F,
                                 const std::vector<std::vector<indexType> >& VF,
                                 const std::vector<std::vector<indexType> >& VFi,
                                 const std::vector<bool>& isB,
                                 Eigen::PlainObjectBase<derivedEnergy>& energy)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    
    //Helper function
    const auto orthogonalVector = [] (const t_V3& u) {
        t_V3 un = u.normalized();
        t_V3 ret;
        if(un[2]+1 == 1) {
            ret << 0., -1., 0.;
        } else {
            const t_V_s a = 1./(1. + un[2]);
            const t_V_s b = -un[0]*un[1]*a;
            ret << 1. - un[0]*un[0]*a, b, -un[0];
        }
        return ret;
    };
    
    //Precomputation
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(doublecast(V), F, faceNormals);
    
    //Actually compute energy
    energy = t_energy(V.rows());
    for(int vert=0; vert<V.rows(); ++vert) {
        //Boundary vertices contribute no energy
#ifdef IGNORE_VALENCE_3
        if(isB[vert] || VF[vert].size()<4)
#else
        if(isB[vert])
#endif
        {
            energy(vert) = 0;
            continue;
        }
        
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        //Do we have a small number of adjacent triangles?
        if(adjacentFaces.size() == 1) {
            t_V3 u = (V.row(F(adjacentFaces[0], (adjacentFacesi[0]+1)%3)) - V.row(vert)).normalized();
            const t_V_s udotn = u.dot(faceNormals.row(adjacentFaces[0]));
            energy(vert) = 0.5*udotn*udotn;
            continue;
        } else if(adjacentFaces.size() == 2) {
            t_V3 u = faceNormals.row(adjacentFaces[0]).cross(faceNormals.row(adjacentFaces[1])).normalized();
            const t_V_s udotn = u.dot(faceNormals.row(adjacentFaces[0]));
            energy(vert) = 0.5*udotn*udotn;
            continue;
        }
        
        //Loop through all face triplets
        t_energy_s& minEnergy = energy(vert);
        minEnergy = INFTY;
        for(int i=0; i<adjacentFaces.size(); ++i) {
            const t_F_i& facei = adjacentFaces[i];
            const t_V3& Ni = faceNormals.row(facei);
            for(int j=i+1; j<adjacentFaces.size(); ++j) {
                const t_F_i& facej = adjacentFaces[j];
                const t_V3& Nj = faceNormals.row(facej);
                for(int k=j+1; k<adjacentFaces.size(); ++k) {
                    const t_F_i& facek = adjacentFaces[k];
                    const t_V3& Nk = faceNormals.row(facek);
                    
                    t_V3 u;
                    bool uComputed = false;
                    
                    //Do we have a degenerate triangle configuration?
                    const t_V_s dotij = Ni.dot(Nj);
                    const t_V_s dotjk = Nj.dot(Nk);
                    const t_V_s dotki = Nk.dot(Ni);
                    if(dotij > 1.-1e-7 && dotjk > 1e-7) {
                        u = orthogonalVector(Ni);
                        uComputed = true;
                    } else if(dotij > 1.-1e-7) {
                        u = Ni.cross(Nk).normalized();
                        uComputed = true;
                    } else if(dotjk > 1.-1e-7) {
                        u = Nj.cross(Ni).normalized();
                        uComputed = true;
                    } else if(dotki > 1.-1e-7) {
                        u = Nk.cross(Nj).normalized();
                        uComputed = true;
                    }
                    
                    //Loop through all sign combinations
                    for(t_V_s p : {-1., 1.}) {
                        for(t_V_s q : {-1., 1.}) {
                            for(t_V_s r : {-1., 1.}) {
                                if(!uComputed) {
                                    u = (q*Nj-p*Ni);
                                    const t_V_s unorm = u.norm();
                                    if(u.norm() < 1e-7)
                                        continue;
                                    u /= unorm;
                                }
                                
                                //Loop through adj faces to find maximum
                                t_energy_s Eijk = 0.;
                                t_F_i maxFace = adjacentFaces[0];
                                for(int f=0; f<adjacentFaces.size(); ++f) {
                                    const t_F_i& face = adjacentFaces[f];
                                    const t_V3& N = faceNormals.row(face);
                                    
                                    const t_V_s udotn = u.dot(N);
                                    const t_V_s localE = 0.5*udotn*udotn;
                                    if(localE > Eijk) {
                                        Eijk = localE;
                                        maxFace = face;
                                    }
                                }
                                
                                //Are we a per-vertex minimum?
                                if(Eijk < minEnergy) {
                                    minEnergy = Eijk;
                                }
                            }
                        }
                    }
                    
                }
                
            }
        }
    }
    
}

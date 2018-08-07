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


#include "maxhingepairs_energy.h"

#include <tools/kopp.h>
#include <tools/triangle_dN.h>
#include <tools/triangle_dTheta.h>

#include <igl/internal_angles.h>
#include <igl/per_face_normals.h>
#include <igl/squared_edge_lengths.h>

#include <vector>
#include <list>

#define INFTY std::numeric_limits<double>::infinity()

#undef TWOSIDES_MAXIMUM


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void maxhingepairs_energy_and_grad(
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
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 2> t_ind;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergy::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    
    //Precomputation
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
    
    
    //Compute energy and edges that form the partition
    energy = t_energy(V.rows());
    t_ind maxNormalIndices(V.rows(), 2);
#ifdef TWOSIDES_MAXIMUM
    t_ind secondMaxNormalIndices(V.rows(), 2);
#else
    t_ind& secondMaxNormalIndices = maxNormalIndices;
#endif
    for(int vert=0; vert<V.rows(); ++vert) {
        if(isB[vert] || VF[vert].size() == 3) { //There can be no hinge if the neighborhood has only size 3
            energy(vert) = 0.;
            continue;
        }
        
        t_V_s& currentEnergy = energy(vert);
        currentEnergy = INFTY;
        const std::vector<indexType>& adjacentFaces = VF[vert];
        
        for(int edge1=0; edge1<adjacentFaces.size()-2; ++edge1) {
            for(int edge2=edge1+2; edge1==0 ? edge2<adjacentFaces.size()-1 : edge2<adjacentFaces.size(); ++edge2) {
                
                t_V_s maxPairEnergy = -1;
                t_ind localMaxNormalIndices(1,2);
                localMaxNormalIndices << edge1, edge2-1;
#ifdef TWOSIDES_MAXIMUM
                t_V_s secondMaxPairEnergy = -1;
                t_ind secondLocalMaxNormalIndices(1,2);
                secondLocalMaxNormalIndices << edge2, (edge1-1+adjacentFaces.size())%adjacentFaces.size();
#else
                t_V_s& secondMaxPairEnergy = maxPairEnergy;
                t_ind& secondLocalMaxNormalIndices = localMaxNormalIndices;
#endif
                
                const auto handleNormalPair = [&] (const t_F_i& n1, const t_F_i& n2, t_V_s& en, t_ind& ind) {
                    const t_V3t& normal1 = faceNormals.row(adjacentFaces[n1]);
                    const t_V3t& normal2 = faceNormals.row(adjacentFaces[n2]);
                    const t_V_s pairEnergy = 0.5*(normal1 - normal2).squaredNorm();
                    
                    if(pairEnergy > en) {
                        en = pairEnergy;
                        ind << n1, n2;
                    }
                };
                
                //Handle normal pairs from edge1 to edge2
                for(int n1=edge1; n1<edge2; ++n1) {
                    for(int n2=n1+1; n2<edge2; ++n2)
                    handleNormalPair(n1, n2, maxPairEnergy, localMaxNormalIndices);
                }
                
                //Handle normal pairs from edge2 to edge1
                for(int n1=edge2; n1<adjacentFaces.size(); ++n1) {
                    for(int n2=n1+1; n2<adjacentFaces.size(); ++n2)
                    handleNormalPair(n1, n2, secondMaxPairEnergy, secondLocalMaxNormalIndices);
                    for(int n2=0; n2<edge1; ++n2)
                    handleNormalPair(n1, n2, secondMaxPairEnergy, secondLocalMaxNormalIndices);
                }
                for(int n1=0; n1<edge1; ++n1) {
                    for(int n2=n1+1; n2<edge1; ++n2)
                    handleNormalPair(n1, n2, secondMaxPairEnergy, secondLocalMaxNormalIndices);
                }
                
                assert(maxPairEnergy>=0 && secondMaxPairEnergy>=0);
                
                //Is this the best configuration so far?
#ifdef TWOSIDES_MAXIMUM
                if(maxPairEnergy+secondMaxPairEnergy < currentEnergy) {
                    currentEnergy = maxPairEnergy+secondMaxPairEnergy;
                    maxNormalIndices.row(vert) << localMaxNormalIndices;
                    secondMaxNormalIndices.row(vert) << secondLocalMaxNormalIndices;
                }
#else
                if(maxPairEnergy < currentEnergy) {
                    currentEnergy = maxPairEnergy;
                    maxNormalIndices.row(vert) << localMaxNormalIndices;
                }
#endif
                
            }
        }
        
    }
    
    
    //Compute energy gradient
    energyGrad = t_energyGrad::Zero(V.rows(), 3);
    for(int vert=0; vert<V.rows(); ++vert) {
        if(isB[vert] || VF[vert].size() == 3) //There can be no hinge if the neighborhood has only size 3
        continue;
        if(energy(vert)+1. == 1.)
            continue;
        
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        const auto process_gradient = [&] (const t_ind& maxNormInd) {
            const t_F_i& n1 = maxNormInd(vert, 0);
            const t_F_i& n2 = maxNormInd(vert, 1);
            
            const int& face1 = adjacentFaces[n1];
            const t_V3& normal1 = faceNormals.row(face1);
            const int& face2 = adjacentFaces[n2];
            const t_V3& normal2 = faceNormals.row(face2);
            
            //Get the vertices nomenclature right
            const int& i1 = vert;
            const int& j1 = F(face1, (adjacentFacesi[n1]+1)%3);
            const int& k1 = F(face1, (adjacentFacesi[n1]+2)%3);
            const int& i2 = vert;
            const int& j2 = F(face2, (adjacentFacesi[n2]+1)%3);
            const int& k2 = F(face2, (adjacentFacesi[n2]+2)%3);
            
            //Compute edges
            t_V3 ejk1 = V.row(k1) - V.row(j1);
            t_V3 eki1 = V.row(i1) - V.row(k1);
            t_V3 eij1 = V.row(j1) - V.row(i1);
            t_V3 ejk2 = V.row(k2) - V.row(j2);
            t_V3 eki2 = V.row(i2) - V.row(k2);
            t_V3 eij2 = V.row(j2) - V.row(i2);
            
            //Double area
            t_V_s dA1 = eij1.cross(eki1).norm();
            t_V_s dA2 = eij2.cross(eki2).norm();
            
            //Compute normal grads
            t_V33 dN1di1, dN1dj1, dN1dk1;
            triangle_dN(eij1, ejk1, eki1, normal1, dA1, dN1di1, dN1dj1, dN1dk1);
            t_V33 dN2di2, dN2dj2, dN2dk2;
            triangle_dN(eij2, ejk2, eki2, normal2, dA2, dN2di2, dN2dj2, dN2dk2);
            
            //Actually add gradients
            energyGrad.row(i1) += (normal1 - normal2).transpose()*dN1di1;
            energyGrad.row(j1) += (normal1 - normal2).transpose()*dN1dj1;
            energyGrad.row(k1) += (normal1 - normal2).transpose()*dN1dk1;
            energyGrad.row(i2) -= (normal1 - normal2).transpose()*dN2di2;
            energyGrad.row(j2) -= (normal1 - normal2).transpose()*dN2dj2;
            energyGrad.row(k2) -= (normal1 - normal2).transpose()*dN2dk2;
        };
        
        process_gradient(maxNormalIndices);
#ifdef TWOSIDES_MAXIMUM
        process_gradient(secondMaxNormalIndices);
#endif
        
    }
    
    
}


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedMinCurvatureDirs>
IGL_INLINE void maxhingepairs_energy(
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
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3t;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedMinCurvatureDirs::Scalar t_mincurvaturedirs_s;
    typedef Eigen::Matrix<t_mincurvaturedirs_s, Eigen::Dynamic, 3> t_MinCurvatureDirs;
    
    //Precomputation
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
    
    
    //Compute energy
    energy = t_energy(V.rows());
    for(int vert=0; vert<V.rows(); ++vert) {
        if(isB[vert] || VF[vert].size() == 3) { //There can be no hinge if the neighborhood has only size 3
            energy(vert) = 0.;
            continue;
        }
        
        t_V_s& currentEnergy = energy(vert);
        currentEnergy = INFTY;
        const std::vector<indexType>& adjacentFaces = VF[vert];
        
        for(int edge1=0; edge1<adjacentFaces.size()-2; ++edge1) {
            for(int edge2=edge1+2; edge1==0 ? edge2<adjacentFaces.size()-1 : edge2<adjacentFaces.size(); ++edge2) {
                
                t_V_s maxPairEnergy = -1;
#ifdef TWOSIDES_MAXIMUM
                t_V_s secondMaxPairEnergy = -1;
#else
                t_V_s& secondMaxPairEnergy = maxPairEnergy;
#endif
                
                const auto handleNormalPair = [&] (const t_F_i& n1, const t_F_i& n2, t_V_s& en) {
                    const t_V3t& normal1 = faceNormals.row(adjacentFaces[n1]);
                    const t_V3t& normal2 = faceNormals.row(adjacentFaces[n2]);
                    const t_V_s pairEnergy = 0.5*(normal1 - normal2).squaredNorm();
                    
                    if(pairEnergy > maxPairEnergy) {
                        maxPairEnergy = pairEnergy;
                    }
                };
                
                //Handle normal pairs from edge1 to edge2
                for(int n1=edge1; n1<edge2; ++n1) {
                    for(int n2=n1+1; n2<edge2; ++n2)
                        handleNormalPair(n1, n2, maxPairEnergy);
                }
                
                //Handle normal pairs from edge2 to edge1
                for(int n1=edge2; n1<adjacentFaces.size(); ++n1) {
                    for(int n2=n1+1; n2<adjacentFaces.size(); ++n2)
                        handleNormalPair(n1, n2, secondMaxPairEnergy);
                    for(int n2=0; n2<edge1; ++n2)
                        handleNormalPair(n1, n2, secondMaxPairEnergy);
                }
                for(int n1=0; n1<edge1; ++n1) {
                    for(int n2=n1+1; n2<edge1; ++n2)
                        handleNormalPair(n1, n2, secondMaxPairEnergy);
                }
                
                assert(maxPairEnergy>=0 && secondMaxPairEnergy>=0);
                
                //Is this the best configuration so far?
#ifdef TWOSIDES_MAXIMUM
                if(maxPairEnergy+secondMaxPairEnergy < currentEnergy) {
                    currentEnergy = maxPairEnergy+secondMaxPairEnergy;
                }
#else
                if(maxPairEnergy < currentEnergy) {
                    currentEnergy = maxPairEnergy;
                }
#endif
                
            }
        }
        
    }
    
    
}

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


#include "hingepairs_energy.h"

#include <tools/kopp.h>
#include <tools/triangle_dN.h>
#include <tools/triangle_dTheta.h>
#include <tools/unitvec_derivative.h>

#include <igl/internal_angles.h>
#include <igl/per_face_normals.h>
#include <igl/squared_edge_lengths.h>

#include <vector>
#include <list>

#define INFTY std::numeric_limits<double>::infinity()

#undef WEIGH_BY_TIPANGLES
#define NORMALIZE_BY_PAIRNUM

template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad>
IGL_INLINE void hingepairs_energy_and_grad(
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
#ifdef WEIGH_BY_TIPANGLES
    t_V angles(V.rows(), 3);
    igl::internal_angles(V, F, angles);
#endif
    
    
    //Compute energy and edges that form the partition
    energy = t_energy(V.rows());
    t_ind partitionIndices(V.rows(), 2);
#ifdef NORMALIZE_BY_PAIRNUM
    t_Vv pairNumbers(V.rows());
#endif
    for(int vert=0; vert<V.rows(); ++vert) {
        if(isB[vert] || VF[vert].size() <= 3) { //There can be no hinge if the neighborhood has only size 3
            energy(vert) = 0.;
            continue;
        }
        
        t_V_s& currentEnergy = energy(vert);
        currentEnergy = INFTY;
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        for(int edge1=0; edge1<adjacentFaces.size()-2; ++edge1) {
            for(int edge2=edge1+2; edge1==0 ? edge2<adjacentFaces.size()-1 : edge2<adjacentFaces.size(); ++edge2) {
                
#ifdef NORMALIZE_BY_PAIRNUM
                t_V_s nPairs = 0.;
#endif
                const auto handleNormalPair = [&] (const t_F_i& n1, const t_F_i& n2) {
#ifdef NORMALIZE_BY_PAIRNUM
                    nPairs += 1.;
#endif
                    
                    const t_V3t& normal1 = faceNormals.row(adjacentFaces[n1]);
                    const t_V3 e1 = V.row(F(adjacentFaces[n1],(adjacentFacesi[n1]+2)%3)) - V.row(F(adjacentFaces[n1],(adjacentFacesi[n1]+1)%3));
                    const t_V3t& normal2 = faceNormals.row(adjacentFaces[n2]);
                    const t_V3 e2 = V.row(F(adjacentFaces[n2],(adjacentFacesi[n2]+2)%3)) - V.row(F(adjacentFaces[n2],(adjacentFacesi[n2]+1)%3));
#ifdef WEIGH_BY_TIPANGLES
                    return (normal1 - normal2).squaredNorm()*angles(adjacentFaces[n1], adjacentFacesi[n1])*angles(adjacentFaces[n2], adjacentFacesi[n2]);
#else
                    return (normal1 - normal2).squaredNorm();
#endif
                };
                
                //ORIGINAL: HANDLE ALL PAIRS
                //Handle normal pairs from edge1 to edge2
                t_V_s edge1sum = 0;
                for(int n1=edge1; n1<edge2; ++n1) {
                    for(int n2=n1+1; n2<edge2; ++n2)
                        edge1sum += handleNormalPair(n1, n2);
                }
                
                //Handle normal pairs from edge2 to edge1
                t_V_s edge2sum = 0;
                for(int n1=edge2; n1<adjacentFaces.size(); ++n1) {
                    for(int n2=n1+1; n2<adjacentFaces.size(); ++n2)
                        edge2sum += handleNormalPair(n1, n2);
                    for(int n2=0; n2<edge1; ++n2)
                        edge2sum += handleNormalPair(n1, n2);
                }
                for(int n1=0; n1<edge1; ++n1) {
                    for(int n2=n1+1; n2<edge1; ++n2)
                        edge2sum += handleNormalPair(n1, n2);
                }
                
                
                //ALTERNATIVE: HANDLE ONLY ADJACENTS
                /*t_V_s edge1sum = 0;
                for(int n1=edge1; n1<edge2-1; ++n1) {
                    int n2=n1+1;
                    edge1sum += handleNormalPair(n1, n2);
                }
                
                t_V_s edge2sum = 0;
                for(int n1=edge2; n1<adjacentFaces.size(); ++n1) {
                    int n2=n1+1;
                    if(n1==adjacentFaces.size()-1) {
                        if(edge1>0)
                            n2=0;
                        else
                            break;
                    }
                    edge2sum += handleNormalPair(n1, n2);
                }
                for(int n1=0; n1<edge1-1; ++n1) {
                    int n2=n1+1;
                    edge2sum += handleNormalPair(n1, n2);
                }*/
                
                
                t_V_s localEnergy = 0.5*(edge1sum + edge2sum);
#ifdef NORMALIZE_BY_PAIRNUM
                localEnergy /= nPairs;
#endif
                if(localEnergy < currentEnergy) {
                    currentEnergy = localEnergy;
                    partitionIndices.row(vert) << edge1, edge2;
#ifdef NORMALIZE_BY_PAIRNUM
                    pairNumbers(vert) = nPairs;
#endif
                }
                
            }
        }
        
    }
    
    
    //Compute energy gradient
    energyGrad = t_energyGrad::Zero(V.rows(), 3);
    for(int vert=0; vert<V.rows(); ++vert) {
        if(isB[vert] || VF[vert].size() <= 3) //There can be no hinge if the neighborhood has only size 3
            continue;
        //if(energy(vert) < 1e-6)
        //    continue;
        
        const std::vector<indexType>& adjacentFaces = VF[vert];
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        const t_F_i& edge1 = partitionIndices(vert, 0);
        const t_F_i& edge2 = partitionIndices(vert, 1);
        const t_F_i& p1tipind = F(adjacentFaces[edge1],(adjacentFacesi[edge1]+1)%3);
        const t_V3 p1 = (V.row(p1tipind) - V.row(vert)).normalized();
        const t_F_i& p2tipind = F(adjacentFaces[edge2],(adjacentFacesi[edge2]+1)%3);
        const t_V3 p2 = (V.row(p2tipind) - V.row(vert)).normalized();
        
        const auto sign = [] (const t_V_s& val)->int {
            if(1. == 1.+val)
            //if(fabs(val) < 1e-8)
                return 0;
            return (0<val) - (val<0);
        };
        
        const auto gradientAddFunc = [&] (const t_F_i& n1, const t_F_i& n2) {
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
            
#ifdef WEIGH_BY_TIPANGLES
            //Compute angle grads
            const t_V_s& theta1 = angles(adjacentFaces[n1], adjacentFacesi[n1]), &theta2 = angles(adjacentFaces[n2], adjacentFacesi[n2]);
            t_V3t dTheta1di, dTheta1dj, dTheta1dk;
            triangle_dTheta(eij1, eki1, normal1, dTheta1di, dTheta1dj, dTheta1dk);
            t_V3t dTheta2di, dTheta2dj, dTheta2dk;
            triangle_dTheta(eij2, eki2, normal2, dTheta2di, dTheta2dj, dTheta2dk);
#endif
            
#ifdef NORMALIZE_BY_PAIRNUM
            t_V_s nc = 1./pairNumbers(vert);
#else
            t_V_s nc = 1.;
#endif
            
            //Actually add gradients
#ifdef WEIGH_BY_TIPANGLES
            energyGrad.row(i1) += nc*(normal1 - normal2).transpose()*dN1di1*theta1*theta2;
            energyGrad.row(j1) += nc*(normal1 - normal2).transpose()*dN1dj1*theta1*theta2;
            energyGrad.row(k1) += nc*(normal1 - normal2).transpose()*dN1dk1*theta1*theta2;
            energyGrad.row(i2) -= nc*(normal1 - normal2).transpose()*dN2di2*theta1*theta2;
            energyGrad.row(j2) -= nc*(normal1 - normal2).transpose()*dN2dj2*theta1*theta2;
            energyGrad.row(k2) -= nc*(normal1 - normal2).transpose()*dN2dk2*theta1*theta2;
            
            energyGrad.row(i1) += nc*(normal1 - normal2).squaredNorm()*dTheta1di*theta2;
            energyGrad.row(j1) += nc*(normal1 - normal2).squaredNorm()*dTheta1dj*theta2;
            energyGrad.row(k1) += nc*(normal1 - normal2).squaredNorm()*dTheta1dk*theta2;
            
            energyGrad.row(i2) += nc*(normal1 - normal2).squaredNorm()*theta1*dTheta2di;
            energyGrad.row(j2) += nc*(normal1 - normal2).squaredNorm()*theta1*dTheta2dj;
            energyGrad.row(k2) += nc*(normal1 - normal2).squaredNorm()*theta1*dTheta2dk;
            
#else
            energyGrad.row(i1) += nc*(normal1 - normal2).transpose()*dN1di1;
            energyGrad.row(j1) += nc*(normal1 - normal2).transpose()*dN1dj1;
            energyGrad.row(k1) += nc*(normal1 - normal2).transpose()*dN1dk1;
            energyGrad.row(i2) -= nc*(normal1 - normal2).transpose()*dN2di2;
            energyGrad.row(j2) -= nc*(normal1 - normal2).transpose()*dN2dj2;
            energyGrad.row(k2) -= nc*(normal1 - normal2).transpose()*dN2dk2;
#endif
        };
        
        //ORIGINAL: HANDLE ALL PAIRS
        //Handle normal pairs from edge1 to edge2
        for(int n1=edge1; n1<edge2; ++n1) {
            for(int n2=n1+1; n2<edge2; ++n2)
                gradientAddFunc(n1, n2);
        }
        
        //Handle normal pairs from edge2 to edge1
        for(int n1=edge2; n1<adjacentFaces.size(); ++n1) {
            for(int n2=n1+1; n2<adjacentFaces.size(); ++n2)
                gradientAddFunc(n1, n2);
            for(int n2=0; n2<edge1; ++n2)
                gradientAddFunc(n1, n2);
        }
        for(int n1=0; n1<edge1; ++n1) {
            for(int n2=n1+1; n2<edge1; ++n2)
                gradientAddFunc(n1, n2);
        }
        
        
        //ALTERNATIVE: HANDLE ONLY ADJACENTS
        /*for(int n1=edge1; n1<edge2-1; ++n1) {
            int n2=n1+1;
            gradientAddFunc(n1, n2);
        }
        
        for(int n1=edge2; n1<adjacentFaces.size()-1; ++n1) {
            int n2=n1+1;
            if(n1==adjacentFaces.size()-1) {
                if(edge1>0)
                    n2=0;
                else
                    break;
            }
            gradientAddFunc(n1, n2);
        }
        for(int n1=0; n1<edge1-1; ++n1) {
            int n2=n1+1;
            gradientAddFunc(n1, n2);
        }*/
        
    }
    
    
}


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy>
IGL_INLINE void hingepairs_energy(
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
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3t;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    
    //Precomputation
    t_V faceNormals(F.rows(), 3);
    igl::per_face_normals(V, F, faceNormals);
#ifdef WEIGH_BY_TIPANGLES
    t_V angles(V.rows(), 3);
    igl::internal_angles(V, F, angles);
#endif
    
    
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
        const std::vector<indexType>& adjacentFacesi = VFi[vert];
        
        for(int edge1=0; edge1<adjacentFaces.size()-2; ++edge1) {
            for(int edge2=edge1+2; edge1==0 ? edge2<adjacentFaces.size()-1 : edge2<adjacentFaces.size(); ++edge2) {
                
#ifdef NORMALIZE_BY_PAIRNUM
                t_V_s nPairs = 0.;
#endif
                const auto handleNormalPair = [&] (const t_F_i& n1, const t_F_i& n2) {
#ifdef NORMALIZE_BY_PAIRNUM
                    nPairs += 1.;
#endif
                    const t_V3t& normal1 = faceNormals.row(adjacentFaces[n1]);
                    const t_V3t& normal2 = faceNormals.row(adjacentFaces[n2]);
#ifdef WEIGH_BY_TIPANGLES
                    return (normal1 - normal2).squaredNorm()*angles(adjacentFaces[n1], adjacentFacesi[n1])*angles(adjacentFaces[n2], adjacentFacesi[n2]);
#else
                    return (normal1 - normal2).squaredNorm();
#endif
                };
                
                //ORIGINAL: HANDLE ALL PAIRS
                //Handle normal pairs from edge1 to edge2
                t_V_s edge1sum = 0;
                for(int n1=edge1; n1<edge2; ++n1) {
                    for(int n2=n1+1; n2<edge2; ++n2)
                        edge1sum += handleNormalPair(n1, n2);
                }
                
                //Handle normal pairs from edge2 to edge1
                t_V_s edge2sum = 0;
                for(int n1=edge2; n1<adjacentFaces.size(); ++n1) {
                    for(int n2=n1+1; n2<adjacentFaces.size(); ++n2)
                        edge2sum += handleNormalPair(n1, n2);
                    for(int n2=0; n2<edge1; ++n2)
                        edge2sum += handleNormalPair(n1, n2);
                }
                for(int n1=0; n1<edge1; ++n1) {
                    for(int n2=n1+1; n2<edge1; ++n2)
                        edge2sum += handleNormalPair(n1, n2);
                }
                
                
                //ALTERNATIVE: HANDLE ONLY ADJACENTS
                /*t_V_s edge1sum = 0;
                for(int n1=edge1; n1<edge2-1; ++n1) {
                    int n2=n1+1;
                    edge1sum += handleNormalPair(n1, n2);
                }
                
                t_V_s edge2sum = 0;
                for(int n1=edge2; n1<adjacentFaces.size()-1; ++n1) {
                    int n2=n1+1;
                    if(n1==adjacentFaces.size()-1) {
                        if(edge1>0)
                            n2=0;
                        else
                            break;
                    }
                    edge2sum += handleNormalPair(n1, n2);
                }
                for(int n1=0; n1<edge1-1; ++n1) {
                    int n2=n1+1;
                    edge2sum += handleNormalPair(n1, n2);
                }*/
                
                t_V_s localEnergy = 0.5*(edge1sum + edge2sum);
#ifdef NORMALIZE_BY_PAIRNUM
                localEnergy /= nPairs;
#endif
                if(localEnergy < currentEnergy)
                    currentEnergy = localEnergy;
                
            }
        }
        
    }
    
}

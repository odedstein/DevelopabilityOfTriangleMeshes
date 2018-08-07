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


#include "measure_once_cut_twice.h"

#include <developableflow/compute_cut.h>
#include <developableflow/flatten_cut.h>
#include <developableflow/hinge_energy.h>

#include <igl/doublearea.h>
#include <igl/edge_lengths.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/boundary_loop.h>

#include <vector>
#include <set>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#define CUT_PERCENTAGE 0.1 //0.95
#define CLOSE_VERTEX_THRESHOLD 0.0005
#define ERROR_TOO_LARGE 1.3

#define MEASUREONCE

#define INFTY std::numeric_limits<double>::infinity()



template <typename derivedV, typename derivedE, typename indexType, typename thresholdType, typename derivedCut, typename derivedRetV, typename derivedRetE, typename derivedRetErr>
IGL_INLINE void measure_once_cut_twice(const Eigen::PlainObjectBase<derivedV>& V,
                                       const Eigen::PlainObjectBase<derivedE>& F,
                                       const Eigen::PlainObjectBase<derivedE>& E,
                                       const Eigen::PlainObjectBase<derivedE>& edgesC,
                                       const Eigen::PlainObjectBase<derivedE>& TT,
                                       const Eigen::PlainObjectBase<derivedE>& TTi,
                                       const std::vector<std::vector<indexType> >& VF,
                                       const std::vector<std::vector<indexType> >& VFi,
                                       const std::vector<bool>& isB,
                                       const thresholdType& cutThreshold,
                                       Eigen::PlainObjectBase<derivedCut>& cut,
                                       Eigen::PlainObjectBase<derivedRetV>& retV,
                                       Eigen::PlainObjectBase<derivedRetE>& retF,
                                       Eigen::PlainObjectBase<derivedRetErr>& error)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef typename derivedE::Scalar t_E_i;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 2> t_E;
    typedef Eigen::Matrix<t_E_i, 1, 2> t_E2;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 3> t_F;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 1> t_Fv;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> t_Ebool;
    typedef typename derivedRetV::Scalar t_retV_s;
    typedef Eigen::Matrix<t_retV_s, Eigen::Dynamic, 3> t_retV;
    typedef typename derivedRetE::Scalar t_retE_i;
    typedef Eigen::Matrix<t_retE_i, Eigen::Dynamic, 3> t_retF;
    typedef typename derivedCut::Scalar t_derivedCut_i;
    typedef Eigen::Matrix<t_derivedCut_i, Eigen::Dynamic, 1> t_derivedCut;
    
    
    //Compute hinge energy (for cost analysis)
    t_Vv hingeEnergy;
    hinge_energy(V, F, VF, VFi, isB, hingeEnergy);
    
    //Sort by energy
    std::vector<std::pair<t_V_s, t_E_i> > costVector;
    costVector.reserve(V.rows());
    for(int vert=0; vert<V.rows(); ++vert)
        costVector.push_back(std::make_pair(hingeEnergy(vert), vert));
    std::sort(costVector.begin(), costVector.end());
    
    //Build puncture list
    std::vector<t_E_i> punctureList;
    for(int v=floor(CUT_PERCENTAGE*V.rows()); v<V.rows(); ++v)
        if(costVector[v].first > cutThreshold)
            punctureList.push_back(costVector[v].second);
    
    compute_cut_erickson(V, F, E, edgesC, TT, TTi, VF, VFi, isB, punctureList, cut);
    flatten_cut(V, F, E, edgesC, TT, TTi, VF, VFi, isB, cut, retV, retF, error);
    
#ifndef MEASUREONCE
    
    //If an error at a face is very large, puncture all its vertices
    //std::set<t_E_i> punctureSet(punctureList.begin(), punctureList.end());
     //for(int face=0; face<error.size(); ++face)
     //if(error(face) > ERROR_TOO_LARGE)
     //for(int j=0; j<3; ++j)
     //punctureSet.insert(F(face,j));
     //punctureList = std::vector<t_E_i>(punctureSet.begin(), punctureSet.end());
     
     //compute_cut_erickson(V, F, E, edgesC, TT, TTi, VF, VFi, isB, punctureList, cut);
    //flatten_cut(V, F, E, edgesC, TT, TTi, VF, VFi, isB, cut, retV, retF, error);
    
    
    //Build retV to V index map
    t_Fv retVtoVmap(retV.rows());
    t_Fv VtoRetVMap(V.rows());
    for(int f=0; f<retF.rows(); ++f) {
        for(int j=0; j<3; ++j) {
            const t_E_i& oldv = F(f,j);
            const t_retE_i& newv = retF(f,j);
            retVtoVmap(newv) = oldv;
            VtoRetVMap(oldv) = newv;
        }
    }
    
    //Define the close vertices as boundary, so they will not be cut again
    std::vector<bool> forbiddenToCut(V.rows(), false);
    for(int v1=0; v1<retV.rows(); ++v1)
        for(int v2=v1+1; v2<retV.rows(); ++v2)
            if(retVtoVmap(v1)==retVtoVmap(v2) && (retV.row(v1)-retV.row(v2)).norm() < CLOSE_VERTEX_THRESHOLD*igl::bounding_box_diagonal(retV))
                for(const int& f : VF[retVtoVmap(v1)])
                    for(int j=0; j<3; ++j)
                        forbiddenToCut[F(f,j)] = true;
    
    /*std::vector<t_E_i> oldCutList;
    igl::boundary_loop(retF, oldCutList);
    punctureList.clear();
    for(const t_E_i& v : oldCutList) {
        const t_E_i& vert = retVtoVmap(v);
        if(!forbiddenToCut[vert])
            punctureList.push_back(vert);
    }*/
    
    std::vector<t_E_i> bdryLoop;
    igl::boundary_loop(retF, bdryLoop);
    std::vector<t_E_i> oldPunctureList = punctureList;
    punctureList.clear();
    for(const t_E_i& vert : oldPunctureList) {
        const t_E_i& v = VtoRetVMap(vert);
        if(forbiddenToCut[vert]) {
            //Walk along the boundary of the flattened mesh until you find a vertex that is not forbidden to cut
            auto iter = std::find(bdryLoop.begin(), bdryLoop.end(), v);
            assert(iter != bdryLoop.end());
            for(++iter; *iter!=vert; ++iter) {
                if(iter == bdryLoop.end())
                    iter = bdryLoop.begin();
                const t_E_i& vertToPush = retVtoVmap(*iter);
                if(!forbiddenToCut[vertToPush]) {
                    punctureList.push_back(vertToPush);
                    break;
                }
            }
        } else {
            punctureList.push_back(vert);
        }
    }
    
    
    compute_cut_erickson(V, F, E, edgesC, TT, TTi, VF, VFi, isB, punctureList, cut);
    flatten_cut(V, F, E, edgesC, TT, TTi, VF, VFi, isB, cut, retV, retF, error);
#endif
    
}

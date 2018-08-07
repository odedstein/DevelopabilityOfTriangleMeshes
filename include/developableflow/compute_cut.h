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


#ifndef DEVELOPABLEFLOW_COMPUTE_CUT_H
#define DEVELOPABLEFLOW_COMPUTE_CUT_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>

//Returns 0 on success, error code otherwise

//Dynamic thresholding, there is no predetermined threshold
template <typename derivedV, typename derivedE, typename indexType, typename derivedCut>
IGL_INLINE int compute_cut_erickson(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                    const Eigen::PlainObjectBase<derivedE>& F, //Faces
                                    const Eigen::PlainObjectBase<derivedE>& E, //Edges
                                    const Eigen::PlainObjectBase<derivedE>& edgesC, //EMAP
                                    const Eigen::PlainObjectBase<derivedE>& TT, //TT from triangle_triangle_adjacency
                                    const Eigen::PlainObjectBase<derivedE>& TTi, //TTi from triangle_triangle_adjacency
                                    const std::vector<std::vector<indexType> >& VF, //VF from vertex_triangle_adjacency
                                    const std::vector<std::vector<indexType> >& VFi, //VFi from vertex_triangle_adjacency
                                    const std::vector<bool>& isB,//Is a vertex a bdry
                                    Eigen::PlainObjectBase<derivedCut>& cut); //indices to cut edges return val

template <typename derivedV, typename derivedE, typename derivedCostScalar, typename indexType, typename derivedCut>
IGL_INLINE int compute_cut_erickson(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                    const Eigen::PlainObjectBase<derivedE>& F, //Faces
                                    const Eigen::PlainObjectBase<derivedE>& E, //Edges
                                    const Eigen::PlainObjectBase<derivedE>& edgesC, //EMAP
                                    const Eigen::PlainObjectBase<derivedE>& TT, //TT from triangle_triangle_adjacency
                                    const Eigen::PlainObjectBase<derivedE>& TTi, //TTi from triangle_triangle_adjacency
                                    const std::vector<std::vector<indexType> >& VF, //VF from vertex_triangle_adjacency
                                    const std::vector<std::vector<indexType> >& VFi, //VFi from vertex_triangle_adjacency
                                    const std::vector<bool>& isB,//Is a vertex a bdry
                                    const derivedCostScalar& costThreshold, //vertices above this absolutely have to be included in the cut
                                    Eigen::PlainObjectBase<derivedCut>& cut); //indices to cut edges return val

template <typename derivedV, typename derivedE, typename indexType, typename derivedCut>
IGL_INLINE int compute_cut_erickson(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                    const Eigen::PlainObjectBase<derivedE>& F, //Faces
                                    const Eigen::PlainObjectBase<derivedE>& E, //Edges
                                    const Eigen::PlainObjectBase<derivedE>& edgesC, //EMAP
                                    const Eigen::PlainObjectBase<derivedE>& TT, //TT from triangle_triangle_adjacency
                                    const Eigen::PlainObjectBase<derivedE>& TTi, //TTi from triangle_triangle_adjacency
                                    const std::vector<std::vector<indexType> >& VF, //VF from vertex_triangle_adjacency
                                    const std::vector<std::vector<indexType> >& VFi, //VFi from vertex_triangle_adjacency
                                    const std::vector<bool>& isB,//Is a vertex a bdry
                                    const std::vector<typename derivedE::Scalar>& punctureList, //A list of the vertices that have to be included in the cuts
                                    Eigen::PlainObjectBase<derivedCut>& cut); //indices to cut edges return val

/*template <typename derivedV, typename derivedF, typename derivedCost, typename derivedCostScalar, typename derivedSegmentIDs>
IGL_INLINE int compute_cut_cgal(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                const Eigen::PlainObjectBase<derivedF>& F, //Faces
                                const Eigen::PlainObjectBase<derivedCost>& cost, //per-vertex cost used to compute cut
                                const derivedCostScalar& costThreshold, //vertices above this absolutely have to be included in the cut
                                Eigen::PlainObjectBase<derivedSegmentIDs>& segmentIds); //A vector containing the segment id for each face*/



//Utility function: Given a mesh, and a cut, return a mesh where the edges have been cut so the resulting mesh can be used with a flattening algorithm
template <typename derivedV, typename derivedE, typename indexType, typename derivedCut, typename derivedRetV, typename derivedRetE>
IGL_INLINE void cut_mesh(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                        const Eigen::PlainObjectBase<derivedE>& F, //Faces
                        const Eigen::PlainObjectBase<derivedE>& E, //Edges
                        const Eigen::PlainObjectBase<derivedE>& edgesC, //EMAP
                        const Eigen::PlainObjectBase<derivedE>& TT, //TT from triangle_triangle_adjacency
                        const Eigen::PlainObjectBase<derivedE>& TTi, //TTi from triangle_triangle_adjacency
                        const std::vector<std::vector<indexType> >& VF, //VF from vertex_triangle_adjacency
                        const std::vector<std::vector<indexType> >& VFi, //VFi from vertex_triangle_adjacency
                        const std::vector<bool>& isB,//Is a vertex a bdry
                        const Eigen::PlainObjectBase<derivedCut>& cut, //indices to cut edges
                        Eigen::PlainObjectBase<derivedV>& retV, //return value V
                        Eigen::PlainObjectBase<derivedE>& retF); //return value F


#ifndef IGL_STATIC_LIBRARY
#  include "compute_cut.cpp"
#endif

#endif

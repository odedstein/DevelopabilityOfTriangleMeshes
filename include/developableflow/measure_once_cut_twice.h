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


#ifndef DEVELOPABLEFLOW_MEASURE_ONCE_CUT_TWICE_H
#define DEVELOPABLEFLOW_MEASURE_ONCE_CUT_TWICE_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>

//Returns 0 on success, error code otherwise

template <typename derivedV, typename derivedE, typename indexType, typename thresholdType, typename derivedCut, typename derivedRetV, typename derivedRetE, typename derivedRetErr>
IGL_INLINE void measure_once_cut_twice(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                       const Eigen::PlainObjectBase<derivedE>& F, //Faces
                                       const Eigen::PlainObjectBase<derivedE>& E, //Edges
                                       const Eigen::PlainObjectBase<derivedE>& edgesC, //EMAP
                                       const Eigen::PlainObjectBase<derivedE>& TT, //TT from triangle_triangle_adjacency
                                       const Eigen::PlainObjectBase<derivedE>& TTi, //TTi from triangle_triangle_adjacency
                                       const std::vector<std::vector<indexType> >& VF, //VF from vertex_triangle_adjacency
                                       const std::vector<std::vector<indexType> >& VFi, //VFi from vertex_triangle_adjacency
                                       const std::vector<bool>& isB, //isB from is_border_vertex
                                       const thresholdType& cutThreshold, //The threshold value used to start the cut with
                                       Eigen::PlainObjectBase<derivedCut>& cut, //A list of edges that make up the cut, indexed into E
                                       Eigen::PlainObjectBase<derivedRetV>& retV, //return mesh vertices
                                       Eigen::PlainObjectBase<derivedRetE>& retF, //return mesh faces
                                       Eigen::PlainObjectBase<derivedRetErr>& error); //flattening error


#ifndef IGL_STATIC_LIBRARY
#  include "measure_once_cut_twice.cpp"
#endif

#endif

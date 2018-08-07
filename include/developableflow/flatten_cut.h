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


#ifndef DEVELOPABLEFLOW_FLATTEN_CUT_H
#define DEVELOPABLEFLOW_FLATTEN_CUT_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>

//Returns 0 on success, error code otherwise

template <typename derivedV, typename derivedE, typename indexType, typename derivedCut, typename derivedRetV, typename derivedRetE>
IGL_INLINE int flatten_cut(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                           const Eigen::PlainObjectBase<derivedE>& F, //Faces
                           const Eigen::PlainObjectBase<derivedE>& E, //Edges
                           const Eigen::PlainObjectBase<derivedE>& edgesC, //EMAP
                           const Eigen::PlainObjectBase<derivedE>& TT, //TT from triangle_triangle_adjacency
                           const Eigen::PlainObjectBase<derivedE>& TTi, //TTi from triangle_triangle_adjacency
                           const std::vector<std::vector<indexType> >& VF, //VF from vertex_triangle_adjacency
                           const std::vector<std::vector<indexType> >& VFi, //VFi from vertex_triangle_adjacency
                           const std::vector<bool>& isB, //isB from is_border_vertex
                           const Eigen::PlainObjectBase<derivedCut>& cut, //indices to cut edges
                           Eigen::PlainObjectBase<derivedRetV>& retV, //return mesh vertices
                           Eigen::PlainObjectBase<derivedRetE>& retF); //return mesh faces


//Version that outputs error

template <typename derivedV, typename derivedE, typename indexType, typename derivedCut, typename derivedRetV, typename derivedRetE, typename derivedRetErr>
IGL_INLINE int flatten_cut(const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                           const Eigen::PlainObjectBase<derivedE>& F, //Faces
                           const Eigen::PlainObjectBase<derivedE>& E, //Edges
                           const Eigen::PlainObjectBase<derivedE>& edgesC, //EMAP
                           const Eigen::PlainObjectBase<derivedE>& TT, //TT from triangle_triangle_adjacency
                           const Eigen::PlainObjectBase<derivedE>& TTi, //TTi from triangle_triangle_adjacency
                           const std::vector<std::vector<indexType> >& VF, //VF from vertex_triangle_adjacency
                           const std::vector<std::vector<indexType> >& VFi, //VFi from vertex_triangle_adjacency
                           const std::vector<bool>& isB, //isB from is_border_vertex
                           const Eigen::PlainObjectBase<derivedCut>& cut, //indices to cut edges
                           Eigen::PlainObjectBase<derivedRetV>& retV, //return mesh vertices
                           Eigen::PlainObjectBase<derivedRetE>& retF, //return mesh faces
                           Eigen::PlainObjectBase<derivedRetErr>& error); //conformal error


#ifndef IGL_STATIC_LIBRARY
#  include "flatten_cut.cpp"
#endif

#endif

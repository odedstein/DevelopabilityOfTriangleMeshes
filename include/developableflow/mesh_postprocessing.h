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


#ifndef DEVELOPABLEFLOW_MESH_POSTPROCESSING_H
#define DEVELOPABLEFLOW_MESH_POSTPROCESSING_H

#include <igl/igl_inline.h>

#include <Eigen/Core>

//Perform the mesh postprocessing, such as removing small triangles etc.
//Returns 0 if no changes happened, returns 1 if a change to F happened


template <typename derivedV, typename derivedF, typename derivedEnergy>
IGL_INLINE int mesh_postprocessing(
                                   Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                   Eigen::PlainObjectBase<derivedF>& F, //Faces
                                   const Eigen::PlainObjectBase<derivedF>& E, //Edges list
                                   const Eigen::PlainObjectBase<derivedF>& edgesC, //EMAP
                                   const Eigen::PlainObjectBase<derivedF>& TT, //TT from triangle_triangle_adjacency
                                   const Eigen::PlainObjectBase<derivedF>& TTi, //TTi from triangle_triangle_adjacency
                                   const std::vector<std::vector<typename OMatrixXi::Scalar> >& VF, //Vertex-face adjacency
                                   const bool& removeUnreferenced=true); //Should we remove unreferenced?



#ifndef IGL_STATIC_LIBRARY
#  include "mesh_postprocessing.cpp"
#endif

#endif

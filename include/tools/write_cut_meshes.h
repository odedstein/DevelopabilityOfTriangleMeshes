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


#ifndef DEVELOPABLEFLOW_WRITE_CUT_MESHES_H
#define DEVELOPABLEFLOW_WRITE_CUT_MESHES_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <vector>
#include <string>


template <typename derivedV, typename derivedF, typename derivedFlatV, typename derivedFlatF>
IGL_INLINE void write_cut_meshes(
                                 const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                 const Eigen::PlainObjectBase<derivedF>& F, //Faces
                                 const Eigen::PlainObjectBase<derivedFlatV>& flatV, //Flattened vertices
                                 const Eigen::PlainObjectBase<derivedFlatF>& flatF, //Flattened faces (rows must correspond to rows in F)
                                 const std::string& filename1 = "", //file to write the output 1 to
                                 const std::string& filename2 = "", //file to write the output 2 to
                                 const std::vector<std::pair<int,int> >& protocol = std::vector<std::pair<int,int> >(), //Protocol
                                 const int& maxTime = 0); //The current timestep number

template <typename derivedV, typename derivedF, typename derivedFlatV, typename derivedFlatF, typename derivedOrigV, typename derivedOrigF>
IGL_INLINE void write_cut_meshes(
                                 const Eigen::PlainObjectBase<derivedV>& V, //Vertices
                                 const Eigen::PlainObjectBase<derivedF>& F, //Faces
                                 const Eigen::PlainObjectBase<derivedFlatV>& flatV, //Flattened vertices
                                 const Eigen::PlainObjectBase<derivedFlatF>& flatF, //Flattened faces (rows must correspond to rows in F)
                                 const Eigen::PlainObjectBase<derivedOrigV>& origV, //Original mesh vertices
                                 const Eigen::PlainObjectBase<derivedOrigF>& origF, //Original mesh faces
                                 const std::string& filename1 = "", //file to write the output 1 to
                                 const std::string& filename2 = "", //file to write the output 2 to
                                 const std::vector<std::pair<int,int> >& protocol = std::vector<std::pair<int,int> >(), //Protocol
                                 const int& maxTime = 0);



#ifndef IGL_STATIC_LIBRARY
#  include "write_cut_meshes.cpp"
#endif

#endif

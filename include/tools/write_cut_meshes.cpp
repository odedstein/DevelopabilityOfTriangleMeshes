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


#include "write_cut_meshes.h"

#include <iostream>
#include <igl/writeOBJ.h>


template <typename derivedV, typename derivedF, typename derivedFlatV, typename derivedFlatF>
IGL_INLINE void write_cut_meshes(
                                 const Eigen::PlainObjectBase<derivedV>& V,
                                 const Eigen::PlainObjectBase<derivedF>& F,
                                 const Eigen::PlainObjectBase<derivedFlatV>& flatV,
                                 const Eigen::PlainObjectBase<derivedFlatF>& flatF,
                                 const std::string& filename1,
                                 const std::string& filename2,
                                 const std::vector<std::pair<int,int> >& protocol,
                                 const int& maxTime)
{
    Eigen::PlainObjectBase<derivedV> dummyV;
    Eigen::PlainObjectBase<derivedF> dummyF;
    write_cut_meshes(V, F, flatV, flatF, dummyV, dummyF, filename1, filename2, protocol, maxTime);
}



template <typename derivedV, typename derivedF, typename derivedFlatV, typename derivedFlatF, typename derivedOrigV, typename derivedOrigF>
IGL_INLINE void write_cut_meshes(
                                 const Eigen::PlainObjectBase<derivedV>& V,
                                 const Eigen::PlainObjectBase<derivedF>& F,
                                 const Eigen::PlainObjectBase<derivedFlatV>& flatV,
                                 const Eigen::PlainObjectBase<derivedFlatF>& flatF,
                                 const Eigen::PlainObjectBase<derivedOrigV>& origV,
                                 const Eigen::PlainObjectBase<derivedOrigF>& origF,
                                 const std::string& filename1,
                                 const std::string& filename2,
                                 const std::vector<std::pair<int,int> >& protocol,
                                 const int& maxTime)
{
    
    //Write cut as mesh
    if(!filename1.empty())
        igl::writeOBJ(filename1, flatV, flatF);
    
    //Write cut as old mesh plus texture coords
    if(!filename2.empty()) {
        assert(F.rows() == flatF.rows() && "The flattened mesh and the original mesh must have the same number of faces");
        
        std::ofstream s(filename2);
        s.precision(17);
        
        s << "# Protocol of actions (translations of numbers to keys: http://www.glfw.org/docs/latest/group__keys.html)" << std::endl;
        for(int p=0; p<protocol.size(); ++p)
            s << "#" << protocol[p].first << " " << protocol[p].second << std::endl;
        if(maxTime>0) {
            s << "#" << maxTime << " " << -1 << std::endl;
        }
        
        for(int v=0; v<V.rows(); ++v)
            s << "v " << V(v,0) << " " << V(v,1) << " " << V(v,2) << std::endl;
        
        for(int v=0; v<flatV.rows(); ++v)
            s << "vt " << flatV(v,0) << " " << flatV(v,1) << std::endl;
        
        if(origF.rows()>0 && origF==F) {
            //Do we also print an original mesh?
            for(int v=0; v<V.rows(); ++v)
                s << "vo " << origV(v,0) << " " << origV(v,1) << " " << origV(v,2) << std::endl;
            for(int f=0; f<F.rows(); ++f)
                s << "f " << F(f,0)+1 << "/" << flatF(f,0)+1 << "/" << origF(f,0)+1 << " " << F(f,1)+1 << "/" << flatF(f,1)+1 << "/" << origF(f,1)+1 << " " << F(f,2)+1 << "/" << flatF(f,2)+1 << "/" << origF(f,2)+1 << std::endl;
        } else {
            for(int f=0; f<F.rows(); ++f)
                s << "f " << F(f,0)+1 << "/" << flatF(f,0)+1 << " " << F(f,1)+1 << "/" << flatF(f,1)+1 << " " << F(f,2)+1 << "/" << flatF(f,2)+1 << std::endl;
        }
    }
    
}

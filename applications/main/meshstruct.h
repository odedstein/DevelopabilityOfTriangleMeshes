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


#ifndef DEVELOPABLEFLOW_MESHSTRUCT_H
#define DEVELOPABLEFLOW_MESHSTRUCT_H

#include <vector>
//#include <list>

#include <igl/all_edges.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/unique_simplices.h>
#include <igl/is_border_vertex.h>

namespace Developables {
    
    struct Mesh {
        OMatrixXs origV;
        OMatrixXi origF;
        OMatrixXs V;
        OMatrixXi F;
        OMatrixXi allE;
        OMatrixXi E;
        OMatrixXi edgesA;
        OMatrixXi edgesC;
        OMatrixXi TT; //triangle-triangle adjacency
        OMatrixXi TTi; //triangle-triangle adjacencyi
        std::vector<std::vector<typename OMatrixXi::Scalar> > VF; //vertex-triangle adjacency
        std::vector<std::vector<typename OMatrixXi::Scalar> > VFi; //vertex-triangle adjacencyi
        std::vector<bool> isB; //is border vertex
        
        void update()
        {
            igl::all_edges(F, allE);
            igl::unique_simplices(allE, E, edgesA, edgesC);
            igl::triangle_triangle_adjacency(F, TT, TTi);
            isB = igl::is_border_vertex(V, F);
            
            //Build sorted vertex triangle adjacency matrix
            igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
            for(int i=0; i<VF.size(); ++i) {
                if(isB[i])
                    continue;
                
                std::vector<typename OMatrixXi::Scalar>& newVF = VF[i];
                std::vector<typename OMatrixXi::Scalar>& newVFi = VFi[i];
                
                //We keep the first face intact, then we rotate over the others
                for(int ind=1; ind < newVF.size(); ++ind) {
                    int nextface = TT(newVF[ind-1], (newVFi[ind-1]+2)%3);
                    newVF[ind] = nextface;
                    for(int j=0; j<3; ++j) {
                        if(F(nextface,j)==i) {
                            newVFi[ind] = j;
                            break;
                        }
                    }
                    
                }
                
            }
            
            //Redo isB
            //for(int i=0; i<V.rows(); ++i) {
            //    if(i!=230)
            //        isB[i]=true;
            //}
        }
        
        Mesh()
        {}
        
        Mesh(const OMatrixXs& iV, const OMatrixXi& iF)
        {
            V = iV;
            F = iF;
            origV = iV;
            origF = iF;
            update();
        }
    };
    
}


#endif

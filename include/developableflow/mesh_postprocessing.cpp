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


#include "mesh_postprocessing.h"

#include <igl/internal_angles.h>
#include <igl/collapse_edge.h>
#include <igl/collapse_small_triangles.h>
#include <igl/edge_flaps.h>
#include <igl/remove_unreferenced.h>
#include <igl/colon.h>
#include <igl/per_face_normals.h>
#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>

#include <tools/triangle_cosdihedral_angle.h>


#define ANGLE_COLLAPSE_THRESHOLD 0.01
#define MIN_COS_DIHEDRAL_ANGLE -0.9999

#undef COLLAPSE_FACES
#define FACE_COLLAPSE_THRESHOLD 2e-9 //1e-8


template <typename derivedV, typename derivedF>
IGL_INLINE int mesh_postprocessing(
                                   Eigen::PlainObjectBase<derivedV>& V,
                                   Eigen::PlainObjectBase<derivedF>& F,
                                   const Eigen::PlainObjectBase<derivedF>& E,
                                   const Eigen::PlainObjectBase<derivedF>& edgesC,
                                   const Eigen::PlainObjectBase<derivedF>& TT,
                                   const Eigen::PlainObjectBase<derivedF>& TTi,
                                   const std::vector<std::vector<typename OMatrixXi::Scalar> >& VF,
                                   const bool& removeUnreferenced)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3t;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef Eigen::Matrix<t_F_i, 1, 3> t_F3t;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, Eigen::Dynamic> t_Fdyn;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 1> t_FV;
    typedef Eigen::Matrix<t_F_i, 1, 2> t_F2;
    
    int retVal = 0;
    
    //Precomputation
    t_V edgesSq;
    igl::squared_edge_lengths(V, F, edgesSq);
    t_V tipAngles;
    igl::internal_angles_using_squared_edge_lengths(edgesSq, tipAngles);
    
    //Flip edges
    int edgesFlipped = 0;
    
    t_F TTm = TT, TTim = TTi;
    const auto flip_edge = [&F, &E, &TTm, &TTim, &edgesFlipped] (const t_F_i& face, const t_F_i& j) {
        const t_F_i& adjFace = TTm(face, (j+1)%3);
        const t_F_i& adjJ = (TTim(face, (j+1)%3)+2)%3;
        
        static Eigen::Vector2i lastFacesFlipped(-1, -1);
        if(face==lastFacesFlipped(0) && adjFace==lastFacesFlipped(1))
            return;
        if(face==lastFacesFlipped(1) && adjFace==lastFacesFlipped(0))
            return;
        lastFacesFlipped << face, adjFace;
        
        const t_F_i v1 = F(face, j);
        const t_F_i v2 = F(face, (j+1)%3);
        assert(v2 == F(adjFace, (adjJ+2)%3) && "These should be the same vertex. There is some problem with the he flip");
        const t_F_i v3 = F(adjFace, adjJ);
        const t_F_i v4 = F(adjFace, (adjJ+1)%3);
        assert(v4 == F(face, (j+2)%3) && "These should be the same vertex. There is some problem with the he flip");
        
        //If this edge already exists, don't flip
        for(int e=0; e<E.rows(); ++e) {
            if(E(e,0)==v1 && E(e,1)==v3)
                return;
            if(E(e,1)==v1 && E(e,0)==v3)
                return;
        }
        
        F.row(face) << v1, v2, v3;
        F.row(adjFace) << v3, v4, v1;
        
        t_F3t ttmface = TTm.row(face);
        t_F3t ttmadjface = TTm.row(adjFace);
        t_F3t ttimface = TTim.row(face);
        t_F3t ttimadjface = TTim.row(adjFace);
        TTm.row(face) << ttmface(j), ttmadjface((adjJ+2)%3), adjFace;
        TTm.row(adjFace) << ttmadjface(j), ttmface((j+2)%3), adjFace;
        TTim.row(face) << ttimface(j), ttimadjface((adjJ+2)%3), 2;
        TTim.row(adjFace) << ttimadjface(j), ttimface((j+2)%3), 2;
        
        TTm(TTm(face, j), TTim(face, j)) = face;
        TTim(TTm(face, j), TTim(face, j)) = 0;
        TTm(TTm(face, (j+1)%3), TTim(face, (j+1)%3)) = face;
        TTim(TTm(face, (j+1)%3), TTim(face, (j+1)%3)) = 1;
        TTm(TTm(adjFace, j), TTim(adjFace, adjJ)) = adjFace;
        TTim(TTm(adjFace, j), TTim(adjFace, adjJ)) = 0;
        TTm(TTm(adjFace, (j+1)%3), TTim(adjFace, (adjJ+1)%3)) = adjFace;
        TTim(TTm(adjFace, (j+1)%3), TTim(adjFace, (adjJ+1)%3)) = 1;
        
        ++edgesFlipped;
    };
    
    const auto flip_largest_edge = [&flip_edge, &V, &F] (const t_F_i& face) {
        t_F_i largestJ = -1;
        t_V_s largestLength = 0;
        for(int j=0; j<3; ++j) {
            const t_V_s length = ( V.row(F(face,(j+2)%3)) - V.row(F(face,(j+1)%3)) ).norm();
            if(length > largestLength) {
                largestLength = length;
                largestJ = j;
            }
        }
        
        flip_edge(face, largestJ);
    };
    
    for(int i=0; i<F.rows(); ++i) {
        for(int j=0; j<3; ++j) {
            const t_V_s& angle = tipAngles(i,j);
            
            //If a vertex only has 3 neighbors or we have a boundary edge, we really don't want to flip anything
            if(angle < ANGLE_COLLAPSE_THRESHOLD && VF[F(i,j)].size()>3 && TTm(i,(j+1)%3)>=0) {
                //If the edge is rather large, flip it, otherwise it will be collapsed later
                const t_V3 ei = V.row(F(i,(j+2)%3)) - V.row(F(i,(j+1)%3));
                const t_V3 ej = V.row(F(i,j)) - V.row(F(i,(j+2)%3));
                const t_V3 ek = V.row(F(i,(j+1)%3)) - V.row(F(i,j));
                if(ei.norm() > 0.1*(ej.norm()+ek.norm()) ) {
                    //Flip the largest edge in this triangle and then go on to the next triangle
                    flip_largest_edge(i);
                    break;
                }
            }
        }
    }
    
    if(edgesFlipped>0)
        retVal = 1;
    
    
    
    
    //Collapse short edges (angle-based)
    std::set<t_F_i> edgesToRemove;
    if(retVal == 0) {
        for(int i=0; i<F.rows(); ++i) {
            for(int j=0; j<3; ++j) {
                const t_V_s& angle = tipAngles(i,j);
                const t_F_i& edge = edgesC(j*F.rows() + i);
                
                //If a vertex only has 3 neighbors or we have a is a boundary edge, we really don't want to collapse anything
                if(angle < ANGLE_COLLAPSE_THRESHOLD && VF[F(i,j)].size()>3 && TTm(i,(j+1)%3)>=0) {
                    edgesToRemove.insert(edge);
                }
            }
        }
        
        if(edgesToRemove.size() > 0) {
            t_V cachedV = V;
            t_F cachedF = F;
            
            t_Fdyn EF, EI;
            igl::edge_flaps(F, E, edgesC, EF, EI);
            
            Eigen::MatrixXd dV = doublecast(V);
            t_Fdyn oldF = F;
            t_FV EMAP = edgesC;
            t_Fdyn nE = E;
            
            for(typename std::set<t_F_i>::iterator iter = edgesToRemove.begin(); iter != edgesToRemove.end(); ++iter) {
                const t_F_i& edge = *iter;
                const Eigen::RowVector3d& p = dV.row(E(edge,0));
                
                igl::collapse_edge(edge, p, dV, oldF, nE, EMAP, EF, EI);
            }
            
            //remove all IGL_COLLAPSE_EDGE_NULL faces
            F = t_F(oldF.rows(), 3);
            int m = 0;
            for(int f = 0; f<oldF.rows(); ++f)
            {
                if(oldF(f,0) != IGL_COLLAPSE_EDGE_NULL || oldF(f,1) != IGL_COLLAPSE_EDGE_NULL || oldF(f,2) != IGL_COLLAPSE_EDGE_NULL) {
                    F.row(m) = oldF.row(f);
                    ++m;
                }
            }
            F.conservativeResize(m, 3);
            
            //Move collapsed edge points to average
            for(typename std::set<t_F_i>::iterator iter = edgesToRemove.begin(); iter != edgesToRemove.end(); ++iter) {
                const t_F_i& edge = *iter;
                const t_F_i& ind0 = E(edge,0);
                const t_F_i& ind1 = E(edge,1);
                const t_V3& v0 = V.row(ind0);
                const t_V3& v1 = V.row(ind1);
                
                //Simple average
                t_V3t p = 0.5*(v0 + v1);
                
                V.row(E(edge,0)) = p;
                V.row(E(edge,1)) = p;
            }
            
            retVal = 1;
        }
    }
    
    
    
    //Remove bad faces (code based on collapse_small_triangles.h)
    int num_face_collapses = 0;
  #ifdef COLLAPSE_FACES
    if(retVal == 0) {
        t_V l = edgesSq.array().sqrt();
        t_Vv doubleAreas;
        igl::doublearea(l, doubleAreas);
        t_V normals;
        igl::per_face_normals(V, F, normals);
        t_Vv dihedralAngles;
        triangle_cosdihedral_angle(edgesC, TT, normals, dihedralAngles);
        Eigen::VectorXi FIM = igl::colon<t_F_i>(0, V.rows()-1);
        int num_edge_collapses = 0;
        for(int f=0; f<F.rows(); ++f) {
            bool dihedralAngleTooSmall = false;
            for(int j=0; j<3; ++j) {
                t_F_i edgeIndex = edgesC(j*F.rows()+f);
                if(dihedralAngles(edgeIndex) < MIN_COS_DIHEDRAL_ANGLE) {
                    dihedralAngleTooSmall = true;
                    break;
                }
            }
            if(doubleAreas(f) < FACE_COLLAPSE_THRESHOLD || dihedralAngleTooSmall) {
                t_V_s minl = 0;
                int minli = -1;
                // Find shortest edge
                for(int e=0; e<3; ++e) {
                    if(minli==-1 || l(f,e)<minl) {
                        minli = e;
                        minl = l(f,e);
                    }
                }
                t_V_s maxl = 0;
                int maxli = -1;
                // Find longest edge
                for(int e=0; e<3; ++e) {
                    if(maxli==-1 || l(f,e)>maxl) {
                        maxli = e;
                        maxl = l(f,e);
                    }
                }
                // Be sure that min and max aren't the same
                maxli = (minli==maxli?(minli+1)%3:maxli);
                
                // Collapse min edge maintaining max edge: i-->j
                int i = maxli;
                int j = ((minli+1)%3 == maxli ? (minli+2)%3: (minli+1)%3);
                assert(i != minli);
                assert(j != minli);
                assert(i != j);
                FIM(F(f,i)) = FIM(F(f,j));
                ++num_edge_collapses;
            }
        }
        
        // Reindex faces
        t_F rF = F;
        // Loop over triangles
        for(int f=0; f<rF.rows(); ++f) {
            for(int i=0; i<rF.cols(); ++i) {
                rF(f,i) = FIM(rF(f,i));
            }
        }
        
        F.resizeLike(rF);
        // Only keep uncollapsed faces
        {
            int ff = 0;
            // Loop over triangles
            for(int f=0; f<rF.rows(); ++f) {
                bool collapsed = false;
                // Check if any indices are the same
                for(int i=0; i<rF.cols(); ++i) {
                    for(int j=i+1; j<rF.cols(); ++j) {
                        if(rF(f,i)==rF(f,j)) {
                            collapsed = true;
                            num_face_collapses++;
                            break;
                        }
                    }
                }
                if(!collapsed) {
                    F.row(ff++) = rF.row(f);
                }
            }
            // Use conservative resize
            F.conservativeResize(ff,F.cols());
        }
        
        if(num_face_collapses > 0) {
            
            retVal = 1;
        }
    }
#endif
    
    //If any processing has happened, resolve errors that could have happened
    if(retVal==1) {
        std::cout << edgesToRemove.size() << " edges collapsed, " << edgesFlipped << " edges flipped, " << num_face_collapses << " faces collapsed." << std::endl;
        
        if(removeUnreferenced) {
            t_FV _1;
            igl::remove_unreferenced(t_V(V), t_F(F), V, F, _1);
        }
    }
    
    
    return retVal;
}

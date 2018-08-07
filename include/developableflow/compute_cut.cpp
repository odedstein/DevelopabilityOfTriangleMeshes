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


#include "compute_cut.h"

#include <igl/per_face_normals.h>

#include <developableflow/hinge_energy.h>
#include <developableflow/hingepairs_energy.h>

#include <vector>
#include <list>
#include <set>
#include <utility>


#define INFTY std::numeric_limits<double>::infinity()


template <typename derivedV, typename derivedE, typename indexType, typename derivedCut>
IGL_INLINE int compute_cut_erickson(const Eigen::PlainObjectBase<derivedV>& V,
                                    const Eigen::PlainObjectBase<derivedE>& F,
                                    const Eigen::PlainObjectBase<derivedE>& E,
                                    const Eigen::PlainObjectBase<derivedE>& edgesC,
                                    const Eigen::PlainObjectBase<derivedE>& TT,
                                    const Eigen::PlainObjectBase<derivedE>& TTi,
                                    const std::vector<std::vector<indexType> >& VF,
                                    const std::vector<std::vector<indexType> >& VFi,
                                    const std::vector<bool>& isB,
                                    Eigen::PlainObjectBase<derivedCut>& cut)
{
    return compute_cut_erickson(V, F, E, edgesC, TT, TTi, VF, VFi, isB, 0.001, cut);
}


template <typename derivedV, typename derivedE, typename derivedCostScalar, typename indexType, typename derivedCut>
IGL_INLINE int compute_cut_erickson(const Eigen::PlainObjectBase<derivedV>& V,
                                    const Eigen::PlainObjectBase<derivedE>& F,
                                    const Eigen::PlainObjectBase<derivedE>& E,
                                    const Eigen::PlainObjectBase<derivedE>& edgesC,
                                    const Eigen::PlainObjectBase<derivedE>& TT,
                                    const Eigen::PlainObjectBase<derivedE>& TTi,
                                    const std::vector<std::vector<indexType> >& VF,
                                    const std::vector<std::vector<indexType> >& VFi,
                                    const std::vector<bool>& isB,
                                    const derivedCostScalar& costThreshold,
                                    Eigen::PlainObjectBase<derivedCut>& cut)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef typename derivedE::Scalar t_E_i;
    
    //Compute hinge energy (for cost analysis)
    t_Vv hingeEnergy;
    hinge_energy(V, F, VF, VFi, isB, hingeEnergy);
    
    //Build puncture list
    std::vector<t_E_i> punctureList;
    for(int vert=0; vert<V.rows(); ++vert) {
        if(!isB[vert] && hingeEnergy(vert) > costThreshold)
            punctureList.push_back(vert);
    }
    
    return compute_cut_erickson(V, F, E, edgesC, TT, TTi, VF, VFi, isB, punctureList, cut);
}


template <typename derivedV, typename derivedE, typename indexType, typename derivedCut>
IGL_INLINE int compute_cut_erickson(const Eigen::PlainObjectBase<derivedV>& V,
                                    const Eigen::PlainObjectBase<derivedE>& F,
                                    const Eigen::PlainObjectBase<derivedE>& E,
                                    const Eigen::PlainObjectBase<derivedE>& edgesC,
                                    const Eigen::PlainObjectBase<derivedE>& TT,
                                    const Eigen::PlainObjectBase<derivedE>& TTi,
                                    const std::vector<std::vector<indexType> >& VF,
                                    const std::vector<std::vector<indexType> >& VFi,
                                    const std::vector<bool>& isB,
                                    const std::vector<typename derivedE::Scalar>& punctureList,
                                    Eigen::PlainObjectBase<derivedCut>& cut)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef typename derivedE::Scalar t_E_i;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 2> t_E;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 3> t_F;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 1> t_Fv;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> t_Ebool;
    typedef typename derivedCut::Scalar t_cut_i;
    typedef Eigen::Matrix<t_cut_i, Eigen::Dynamic, 1> t_cut;
    
    int retVal = 0;
    
    
    //Compute the hinge energy for the min curvature vecs
    t_Vv hingeEnergy;
    t_V minCurvatureVecs;
    hinge_energy(V, F, VF, VFi, isB, hingeEnergy, minCurvatureVecs);
    
    //Constants
    const t_E_i npunctures = punctureList.size();
    
    
    //Halfedge data structure
    struct Halfedge {
        t_E_i face;
        t_E_i j;
        t_E_i normalj() const {return (j+2)%3;}
        Halfedge(const t_E_i& i_f, const t_E_i& i_j) : face(i_f), j(i_j) {}
        Halfedge() {}
        Halfedge next() const {return Halfedge(face, (j+1)%3);}
        bool operator==(const Halfedge& he) const {return face==he.face && j==he.j;}
    };
    typedef std::vector<Halfedge> Halfedges;
    const auto heedge = [&F, &edgesC](const Halfedge& he) {
        return edgesC(he.face + F.rows()*he.normalj());
    };
    const auto hetail = [&F](const Halfedge& he) {
        return F(he.face, he.j);
    };
    const auto hehead = [&F](const Halfedge& he) {
        return F(he.face, (he.j+1)%3);
    };
    const auto heflip = [&TT, &TTi](const Halfedge& he) {
        return Halfedge(TT(he.face,he.j), TTi(he.face,he.j));
    };
    const auto hesatvert = [&VF, &VFi](const t_E_i& vert) {
        const std::vector<indexType>& adjFaces = VF[vert];
        const std::vector<indexType>& adjFacesi = VFi[vert];
        Halfedges halfedges(adjFaces.size());
        for(int i=0; i<halfedges.size(); ++i)
            halfedges[i] = Halfedge(adjFaces[i], adjFacesi[i]);
        return halfedges;
    };
    const auto hecost = [&isB, &hehead, &hetail, &minCurvatureVecs, &V](const Halfedge& he) {
        const t_E_i& v = hetail(he);
        const t_E_i& w = hehead(he);
        if(isB[v] && isB[w])
            return INFTY;
        t_V3 edgevec = V.row(w) - V.row(v);
        const t_V_s enorm = edgevec.norm(); edgevec /= enorm;
        return minCurvatureVecs.row(v).cross(edgevec).norm() + enorm;
        //return enorm;
    };
    const auto sharpness = [&hecost](const Halfedges& string) {
        t_V_s sharp = 0.;
        for(int i=0; i<string.size(); ++i)
            sharp += hecost(string[i]);
        return sharp;
    };

    
    //Build adjacency matrix for puncture list. Cost of each edge is determined by dijkstra's shortest path
    std::vector<std::vector<Halfedges> > adjacency(npunctures);
    for(int puncture1=0; puncture1<npunctures; ++puncture1) {
        const t_E_i& p1ind = punctureList[puncture1];
        std::vector<Halfedges>& localAdjacency = adjacency[puncture1];
        std::vector<bool> localPsVisited(npunctures);
        for(int i=0; i<npunctures; ++i)
            localPsVisited[i] = false;
        
        //Use Dijkstra to find the shortest path from puncture1 to the localPs. Based on libigl's dijkstra.h
        t_Vv minDistance = t_Vv::Constant(V.rows(), INFTY);
        minDistance(p1ind) = 0.;
        t_Fv prev = t_Fv::Constant(V.rows(), -1);
        std::set<std::pair<t_V_s, t_E_i> > vertexQueue;
        vertexQueue.insert(std::make_pair(minDistance(p1ind), p1ind));
        
        while(!vertexQueue.empty()) {
            const t_V_s dist = vertexQueue.begin()->first;
            const t_E_i u = vertexQueue.begin()->second;
            vertexQueue.erase(vertexQueue.begin());
            
            bool allPsVisited = true;
            for(int i=0; i<npunctures; ++i) {
                if(punctureList[i]==u)
                    localPsVisited[i] = true;
                if(!localPsVisited[i]) {
                    allPsVisited = false;
                    break;
                }
            }
            if(allPsVisited)
                break;
            
            //Visit each halfedge exiting u
            for(Halfedge he : hesatvert(u)) {
                const t_E_i& v = hehead(he);
                const t_V_s distThroughU = dist + hecost(he);
                if(distThroughU < minDistance(v)) {
                    vertexQueue.erase(std::make_pair(minDistance(v), v));
                    minDistance(v) = distThroughU;
                    prev(v) = u;
                    vertexQueue.insert(std::make_pair(minDistance(v), v));
                }
            }
        }
        
        //Build the local adjacency list
        localAdjacency = std::vector<Halfedges>(npunctures);
        for(int target=0; target<npunctures; ++target) {
            if(puncture1 == target)
                continue;
            
            Halfedges reversePath;
            t_E_i currentVert = punctureList[target];
            while(currentVert != p1ind) {
                t_E_i nextVert = prev(currentVert);
                for(Halfedge he : hesatvert(nextVert)) {
                    if(hehead(he) == currentVert) {
                        reversePath.push_back(he);
                        break;
                    }
                }
                currentVert = nextVert;
            }
            localAdjacency[target] = Halfedges(reversePath.rbegin(), reversePath.rend());
        }
    }
    
    //Using the adjacency matrix from earlier, compute a minimum spanning tree. Uses Prim's algorithm (Wikipedia)
    std::set<t_E_i> cutEdges;
    t_Vv pC = t_Vv::Constant(npunctures, INFTY);
    t_Fv pE = t_Fv::Constant(npunctures, -1);
    std::set<std::pair<t_V_s, t_E_i> > pQ;
    for(int i=0; i<npunctures; ++i)
        pQ.insert(std::make_pair(pC(i), i));
    
    while(!pQ.empty()) {
        const t_E_i v = pQ.begin()->second;
        pQ.erase(pQ.begin());
        
        if(pE(v) != -1) {
            //This edge will be included in the final cut. Mark these edges.
            for(Halfedge h : adjacency[pE(v)][v])
                cutEdges.insert(heedge(h));
        }
        
        for(int w=0; w<npunctures; ++w) {
            if(v==w)
                continue;
            if(pQ.find(std::make_pair(pC(w), w)) != pQ.end()) {
                const t_V_s weight = sharpness(adjacency[v][w]);
                if(weight < pC(w)) {
                    pQ.erase(std::make_pair(pC(w), w));
                    pC(w) = weight;
                    pE(w) = v;
                    pQ.insert(std::make_pair(pC(w), w));
                }
            }
        }
    }
    
    cut = t_cut(cutEdges.size());
    int ind = 0;
    for(t_E_i edge : cutEdges)
        cut(ind++) = edge;
    
    
    return retVal;
}


/*#include <igl/copyleft/cgal/extract_feature.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>


template <typename derivedV, typename derivedF, typename derivedCost, typename derivedCostScalar, typename derivedSegmentIDs>
IGL_INLINE int compute_cut_cgal(const Eigen::PlainObjectBase<derivedV>& V,
                                const Eigen::PlainObjectBase<derivedF>& F,
                                const Eigen::PlainObjectBase<derivedCost>& cost,
                                const derivedCostScalar& costThreshold,
                                Eigen::PlainObjectBase<derivedSegmentIDs>& segmentIds)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedCost::Scalar t_cost_s;
    typedef Eigen::Matrix<t_cost_s, Eigen::Dynamic, 1> t_cost;
    typedef typename derivedSegmentIDs::Scalar t_segmentid_i;
    typedef Eigen::Matrix<t_segmentid_i, Eigen::Dynamic, 1> t_SegmentID;
    
    int retVal = 0;
    
    //Code from CGAL doc, http://doc.cgal.org/latest/Surface_mesh_segmentation/index.html
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
    
    //Transform to polyhedron
    Polyhedron polymesh;
    igl::copyleft::cgal::mesh_to_polyhedron(V, F, polymesh);
    
    // assign id field for each facet
    t_F_i facet_id = 0;
    for(Polyhedron::Facet_iterator facet_it = polymesh.facets_begin(); facet_it != polymesh.facets_end(); ++facet_it)
        facet_it->id() = facet_id++;
    
    // create a property-map for SDF values
    typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
    Facet_double_map internal_map;
    boost::associative_property_map<Facet_double_map> sdf_property_map(internal_map);
    // compute SDF values
    const std::size_t number_of_rays = 25;  // cast 25 rays per facet
    const double cone_angle = 2.0 / 3.0 * CGAL_PI; // set cone opening-angle
    CGAL::sdf_values(polymesh, sdf_property_map, cone_angle, number_of_rays, false);
    std::pair<double, double> min_max_sdf = CGAL::sdf_values_postprocessing(polymesh, sdf_property_map);
    
    // create a property-map for segment-ids
    typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);
    const std::size_t number_of_clusters = 5;
    const double smoothing_lambda = 0.1;
    std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(polymesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);
    
    //Write segments into eigen vector
    segmentIds = t_SegmentID(F.rows());
    std::cout << "Number of segments: " << number_of_segments << std::endl;
    for(Polyhedron::Facet_const_iterator facet_it = polymesh.facets_begin(); facet_it != polymesh.facets_end(); ++facet_it) {
        segmentIds(facet_it->id()) = segment_property_map[facet_it];
    }
    std::cout << std::endl;
    
    return retVal;
}*/



template <typename derivedV, typename derivedE, typename indexType, typename derivedCut, typename derivedRetV, typename derivedRetE>
IGL_INLINE void cut_mesh(const Eigen::PlainObjectBase<derivedV>& V,
                        const Eigen::PlainObjectBase<derivedE>& F,
                        const Eigen::PlainObjectBase<derivedE>& E,
                        const Eigen::PlainObjectBase<derivedE>& edgesC,
                        const Eigen::PlainObjectBase<derivedE>& TT,
                        const Eigen::PlainObjectBase<derivedE>& TTi,
                        const std::vector<std::vector<indexType> >& VF,
                        const std::vector<std::vector<indexType> >& VFi,
                        const std::vector<bool>& isB,
                        const Eigen::PlainObjectBase<derivedCut>& cut,
                        Eigen::PlainObjectBase<derivedRetV>& retV,
                        Eigen::PlainObjectBase<derivedRetE>& retF)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef typename derivedE::Scalar t_E_i;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 2> t_E;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 3> t_F;
    typedef Eigen::Matrix<t_E_i, Eigen::Dynamic, 1> t_Fv;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> t_Ebool;
    typedef typename derivedCut::Scalar t_cut_i;
    typedef Eigen::Matrix<t_cut_i, Eigen::Dynamic, 1> t_cut;
    typedef typename derivedRetV::Scalar t_retV_s;
    typedef Eigen::Matrix<t_retV_s, Eigen::Dynamic, 3> t_retV;
    typedef typename derivedRetE::Scalar t_retE_i;
    typedef Eigen::Matrix<t_retE_i, Eigen::Dynamic, 3> t_retF;
    
    
    //Halfedge data structure
    struct Halfedge {
        t_E_i face;
        t_E_i j;
        t_E_i normalj() const {return (j+2)%3;}
        Halfedge(const t_E_i& i_f, const t_E_i& i_j) : face(i_f), j(i_j) {}
        Halfedge() {}
        Halfedge next() const {return Halfedge(face, (j+1)%3);}
        bool operator==(const Halfedge& he) const {return face==he.face && j==he.j;}
        bool operator!=(const Halfedge& he) const {return face!=he.face || j!=he.j;}
    };
    typedef std::vector<Halfedge> Halfedges;
    const auto heedge = [&F, &edgesC](const Halfedge& he) {
        return edgesC(he.face + F.rows()*he.normalj());
    };
    const auto hetail = [&F](const Halfedge& he) {
        return F(he.face, he.j);
    };
    const auto hehead = [&F](const Halfedge& he) {
        return F(he.face, (he.j+1)%3);
    };
    const auto heflip = [&TT, &TTi](const Halfedge& he) {
        return Halfedge(TT(he.face,he.j), TTi(he.face,he.j));
    };
    const auto hesatvert = [&VF, &VFi](const t_E_i& vert) {
        const std::vector<indexType>& adjFaces = VF[vert];
        const std::vector<indexType>& adjFacesi = VFi[vert];
        Halfedges halfedges(adjFaces.size());
        for(int i=0; i<halfedges.size(); ++i)
            halfedges[i] = Halfedge(adjFaces[i], adjFacesi[i]);
        return halfedges;
    };
    
    
    //Compute an edge cut bool
    t_Ebool cutBool = t_Ebool::Constant(E.rows(), false);
    for(int i=0; i<cut.size(); ++i)
        cutBool(cut(i)) = true;
    
    std::vector<bool> bdryOrCut = isB;
    //Assign cut dofs
    retF.resize(F.rows(), 3);
    int nDof = 0;
    for(int vert=0; vert<V.rows(); ++vert) {
        Halfedges hesAtVert = hesatvert(vert);
        
        //Find a cut edge if it exists
        Halfedge he = hesAtVert[0];
        for(Halfedge h : hesAtVert) {
            if(cutBool(heedge(h))) {
                he = h;
                bdryOrCut[vert] = true;
                break;
            }
        }
        
        //Walk back one
        he = heflip(he).next();
        Halfedge h = he;
        do {
            retF(h.face, h.j) = nDof;
            if(cutBool(heedge(h)) && heflip(h).next()!=he)
                ++nDof;
            h = heflip(h).next();
        } while(h != he);
        ++nDof;
    }
    
    retV.resize(nDof, 3);
    for(int face=0; face<F.rows(); ++face)
        for(int j=0; j<3; ++j)
            retV.row(retF(face,j)) = V.row(F(face,j));
    
            
}

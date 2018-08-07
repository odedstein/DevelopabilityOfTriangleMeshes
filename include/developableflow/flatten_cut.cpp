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


#include "flatten_cut.h"

#include <igl/doublearea.h>
#include <igl/edge_lengths.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include <igl/lscm.h>
#include <igl/boundary_loop.h>
#include <igl/bounding_box_diagonal.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>


#define LSCM_FLATTEN 1
#define SCP_FLATTEN 2
#define BFF_FLATTEN 3
#define BLENDER_FLATTEN 4
#define FLATTENING_METHOD SCP_FLATTEN

#if FLATTENING_METHOD==BFF_FLATTEN
using namespace std;
#include <DenseMatrix.h>
#include <Mesh.h>
#include <Bff.h>
#endif

#if FLATTENING_METHOD==BLENDER_FLATTEN

#endif


template <typename derivedV, typename derivedE, typename indexType, typename derivedCut, typename derivedRetV, typename derivedRetE>
IGL_INLINE int flatten_cut(const Eigen::PlainObjectBase<derivedV>& V,
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
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 2> t_V2;
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
    typedef Eigen::Matrix<t_retV_s, Eigen::Dynamic, Eigen::Dynamic> t_retVd;
    typedef typename derivedRetE::Scalar t_retE_i;
    typedef Eigen::Matrix<t_retE_i, Eigen::Dynamic, 3> t_retF;
    
    
    typedef Eigen::SparseMatrix<t_V_s> t_sparse;
    typedef Eigen::Triplet<t_V_s> t_triplet;
    typedef std::complex<t_V_s> t_c;
    typedef Eigen::SparseMatrix<t_c> t_sparse_c;
    typedef Eigen::Triplet<t_c> t_triplet_c;
    typedef Eigen::Matrix<t_c, Eigen::Dynamic, 1> t_vec_c;
    const t_c I(0, 1);
    
    
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
    
    
    int retVal = 0;
    
    //Precompute areas and edge lengths
    t_V edgesSq;
    igl::squared_edge_lengths(V, F, edgesSq);
    t_V edges = edgesSq.array().sqrt();
    t_V angles;
    igl::internal_angles_using_squared_edge_lengths(edgesSq, angles);
    t_Vv doubleAreas;
    igl::doublearea(edges, doubleAreas);
    
    //Compute an edge cut bool
    t_Ebool cutBool = t_Ebool::Constant(E.rows(), false);
    for(int i=0; i<cut.size(); ++i)
        cutBool(cut(i)) = true;
    
    std::vector<bool> bdryOrCut = isB;
    //Assign cut dofs
    t_F cutIndices(F.rows(), 3);
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
            cutIndices(h.face, h.j) = nDof;
            if(cutBool(heedge(h)) && heflip(h).next()!=he)
                ++nDof;
            h = heflip(h).next();
        } while(h != he);
        ++nDof;
    }
    t_Fv dofToVert(nDof);
    for(int face=0; face<F.rows(); ++face) {
        for(int j=0; j<3; ++j) {
            dofToVert(cutIndices(face,j)) = F(face,j);
        }
    }
    
    //Create cut mesh
    t_V cutV(nDof, 3);
    for(int face=0; face<F.rows(); ++face)
        for(int j=0; j<3; ++j)
            cutV.row(cutIndices(face,j)) = V.row(F(face,j));
    
    
#if FLATTENING_METHOD==LSCM_FLATTEN
    
    Eigen::VectorXi bnd, b(2,1);
    igl::boundary_loop(cutIndices,bnd);
    b(0) = bnd(0);
    Eigen::MatrixXd bcr(2,3), bc(2,2);
    bcr.row(0) = cutV.row(b(0));
    for(int i=1; i<bnd.size(); ++i) {
        b(1) = bnd(i);
        bcr.row(1) = cutV.row(b(1));
        if((bcr.row(1)-bcr.row(0)).norm() > 1e-4)
            break;
    }
    bc.row(0).setZero();
    bc.row(1) << 0, (bcr.row(1)-bcr.row(0)).norm() / igl::bounding_box_diagonal(V);
    Eigen::MatrixXd retVd;
    igl::lscm(Eigen::MatrixXd(cutV), Eigen::MatrixXi(cutIndices), b, bc, retVd);
    retVd.conservativeResize(retVd.rows(), 3);
    retVd.col(2).setZero();
    retF = cutIndices;
    retV = retVd;
    
#elif FLATTENING_METHOD==SCP_FLATTEN
    
    //Compute SCP matrices
    t_sparse_c A(nDof, nDof), B(nDof, nDof);
    std::vector<t_triplet_c> tripletListA, tripletListB;
    for(int face=0; face<F.rows(); ++face) {
        const t_V_s A = 0.5*doubleAreas(face);
        for(int j=0; j<3; ++j) {
            int i = cutIndices(face,j);
            int i1 = cutIndices(face,(j+1)%3);
            int i2 = cutIndices(face,(j+2)%3);
            
            //Compute cotan
            const t_V_s& u = edges(face,j);
            const t_V_s& v = edges(face,(j+1)%3);
            const t_V_s& w = edges(face,(j+2)%3);
            const t_V_s diam = u*v*w/(2.*A);
            const t_V_s sinu = u/diam;
            const t_V_s cosu = (v*v+w*w-u*u)/(2.*v*w);
            const t_V_s cotanh = cosu/sinu;
            
            tripletListA.emplace_back(i1, i2, -0.25*cotanh);
            tripletListA.emplace_back(i2, i1, -0.25*cotanh);
            tripletListA.emplace_back(i1, i1, 0.25*cotanh);
            tripletListA.emplace_back(i2, i2, 0.25*cotanh);
            if(cutBool(edgesC(face + F.rows()*j))) {
                tripletListA.emplace_back(i1, i2, 0.25*I);
                tripletListA.emplace_back(i2, i1, -0.25*I);
            }
            
            tripletListB.emplace_back(i, i, 1./3.*A);
        }
    }
    //For semipositive definiteness
    for(int i=0; i<nDof; ++i)
        tripletListA.emplace_back(i, i, 1e-8);
    A.setFromTriplets(tripletListA.begin(), tripletListA.end());
    Eigen::SimplicialLDLT<t_sparse_c> solver(A);
    B.setFromTriplets(tripletListB.begin(), tripletListB.end());
    
    //Find smallest eigenvalue
    t_vec_c y = t_vec_c::Random(nDof);
    y /= (t_V_s) nDof;
    
    const int MAX_ITER = 100; //10*V.rows();
    const t_V_s EPS = 1e-7;
    t_V_s residual;
    for(int i=0; i<MAX_ITER; ++i) {
        t_vec_c w = B*y;
        y = solver.solve(w);
        
        t_c mean = y.sum() / ((t_V_s) y.size());
        y.array() -= mean;
        y /= sqrt(y.dot(B*y));
        
        t_V_s residual = (A*y - y.dot(A*y)*y).norm();
        if(residual < EPS)
            break;
    }
    
    //Create new mesh
    retV = t_retV(nDof, 3);
    for(int i=0; i<nDof; ++i) {
        retV(i,0) = real(y(i));
        retV(i,1) = imag(y(i));
        retV(i,2) = 0;
    }
    retF = t_retF(F.rows(), 3);
    for(int face=0; face<F.rows(); ++face) {
        for(int j=0; j<3; ++j) {
            retF(face, j) = cutIndices(face, j);
        }
    }
    
#elif FLATTENING_METHOD==BFF_FLATTEN
    
    //Create bff mesh
    std::string tmpfile = std::tmpnam(nullptr);
    igl::writeOBJ(tmpfile, cutV, cutIndices);
    Mesh mesh;
    mesh.read(tmpfile);
    
    //Use Bff library
    BFF bff(mesh);
    DenseMatrix u(bff.data->bN);
    bff.flatten(u, true);
    
    //Extract data from bff to retV, retF
    retV.resize(mesh.vertices.size(), 3);
    t_retF prRetF(mesh.faces.size(), 3);
    for(int f=0; f<mesh.faces.size(); ++f) {
        HalfEdgeIter he = mesh.faces[f].he;
        for(int j=0; j<3; ++j) {
            if(!mesh.faces[f].isReal())
                prRetF(f,j) = -1;
            else
                prRetF(f,j) = he->vertex->index;
            
            const Vector& uv = he->wedge()->uv;
            retV.row(he->next->next->vertex->index) << 0.5*(uv.x+1.), 0.5*(uv.y+1.), 0;
            
            he = he->next;
        }
    }
    
    //Remove the not-real face
    retF.resize(prRetF.rows(), 3);
    int m=0;
    for(int f=0; f<prRetF.rows(); ++f) {
        if(prRetF(f,0)!=-1) {
            retF.row(m) = prRetF.row(f);
            ++m;
        }
    }
    
    //Delete the temp file from earlier
    remove(tmpfile.c_str());
    
#elif FLATTENING_METHOD==BLENDER_FLATTEN
    
    //Write mesh to tmp file and import to blender
    std::string tmpMesh = std::tmpnam(nullptr);
    igl::writeOBJ(tmpMesh, cutV, cutIndices);
    
    //Write blender python script to a file
    std::string tmpPython = std::tmpnam(nullptr);
    std::string tmpOutput = std::tmpnam(nullptr);
    std::ofstream out(tmpPython);
    out << R"(
import bpy
import os
    
for obj in bpy.data.objects:
    if obj.type == 'MESH':
        obj.select = True
        bpy.context.scene.objects.active = obj
        bpy.ops.object.delete()

inputFile = ')"
    << tmpMesh <<
    R"('
outputFile = ')"
    << tmpOutput <<
    R"('
bpy.ops.import_scene.obj(filepath=inputFile)
objName = os.path.splitext(os.path.basename(inputFile))[0]
obj = bpy.context.selected_objects[0]
obj.select = True
bpy.context.scene.objects.active = obj
bpy.ops.object.mode_set(mode='EDIT')
bpy.ops.uv.smart_project(angle_limit=66.0, island_margin=0.01, use_aspect=True, stretch_to_bounds=False)
bpy.ops.object.mode_set(mode='OBJECT')
bpy.ops.export_scene.obj(filepath=outputFile)
)";
    out.close();
    
    //Run blender
    std::system((std::string("/Applications/blender.app/Contents/MacOS/blender --background --python ") + tmpPython + std::string(" >/dev/null 2>/dev/null")).c_str());
    std::cout << "out: " << tmpOutput << std::endl;
    
    //Read in file
    t_V dummy1, dummy2;
    t_F dummy3, dummy4;
    t_V2 retVd;
    igl::readOBJ(tmpOutput, dummy1, retVd, dummy2, dummy3, retF, dummy4);
    retV.resize(retVd.rows(), 3);
    retV.block(0,0,retVd.rows(),retVd.cols()) = retVd;
    retV.col(2).setZero();
    
    //Delete the temp files from earlier
    std::remove(tmpMesh.c_str());
    std::remove(tmpPython.c_str());
    //std::remove(tmpOutput.c_str());
    
#endif
    
    return retVal;
}


template <typename derivedV, typename derivedE, typename indexType, typename derivedCut, typename derivedRetV, typename derivedRetE, typename derivedRetErr>
IGL_INLINE int flatten_cut(const Eigen::PlainObjectBase<derivedV>& V,
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
                           Eigen::PlainObjectBase<derivedRetE>& retF,
                           Eigen::PlainObjectBase<derivedRetErr>& error)
{
    
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef Eigen::Matrix<t_V_s, 3, 1> t_V3;
    typedef Eigen::Matrix<t_V_s, 3, 3> t_V33;
    typedef typename derivedRetErr::Scalar t_error_s;
    typedef Eigen::Matrix<t_error_s, Eigen::Dynamic, 1> t_error;
    
    int retVal = flatten_cut(V, F, E, edgesC, TT, TTi, VF, VFi, isB, cut, retV, retF);
    
    //Compute conformal error
    auto quasi_conformal_dist = [](t_V3 p1, t_V3 p2, t_V3 p3, t_V3 q1, t_V3 q2, t_V3 q3) {
        //compute edge vectors
        const t_V3 u1 = p2 - p1;
        const t_V3 u2 = p3 - p1;
        const t_V3 v1 = q2 - q1;
        const t_V3 v2 = q3 - q1;
        
        //compute orthonormal bases
        const t_V3 e1 = u1 / u1.norm();
        t_V3 e2 = u2 - u2.dot(e1)*e1;
        e2 /= e2.norm();
        
        const t_V3 f1 = v1 / v1.norm();
        t_V3 f2 = v2 - v2.dot(f1)*f1;
        f2 /= f2.norm();
        
        // project onto bases
        p1 = t_V3( 0., 0., 0. );
        p2 = t_V3( u1.dot(e1), u1.dot(e2), 0. );
        p3 = t_V3( u2.dot(e1), u2.dot(e2), 0. );
        
        q1 = t_V3( 0., 0., 0. );
        q2 = t_V3( v1.dot(f1), v1.dot(f2), 0. );
        q3 = t_V3( v2.dot(f1), v2.dot(f2), 0. );
        
        t_V_s A = (p2-p1).cross(p3-p1).norm();
        
        t_V3 Ss = (q1*(p2(1)-p3(1))+q2*(p3(1)-p1(1))+q3*(p1(1)-p2(1)))/(2.*A);
        t_V3 St = (q1*(p3(0)-p2(0))+q2*(p1(0)-p3(0))+q3*(p2(0)-p1(0)))/(2.*A);
        t_V_s a = Ss.dot(Ss);
        t_V_s b = Ss.dot(St);
        t_V_s c = St.dot(St);
        t_V_s Gamma = sqrt(.5*((a+c)+sqrt(pow(a-c, 2)+4.*b*b)));
        t_V_s gamma = sqrt(.5*((a+c)-sqrt(pow(a-c, 2)+4.*b*b)));
        
        if( Gamma < gamma ) {
            t_V_s og = Gamma;
            Gamma = gamma;
            gamma = og;
        }
        
        return Gamma/gamma;
    };
    
    /*auto svd_error = [](t_V3 p1, t_V3 p2, t_V3 p3, t_V3 q1, t_V3 q2, t_V3 q3) {
        t_V33 p, q;
        p.col(0) = p2-p1;
        p.col(1) = p3-p1;
        p.col(2) = p.col(0).cross(p.col(1)).normalized();
        q.col(0) = q2-q1;
        q.col(1) = q3-q1;
        q.col(2) = q.col(0).cross(q.col(1)).normalized();
        
        t_V33 mat = q*p.inverse();
        Eigen::JacobiSVD<t_V33> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
        t_V3 singularValues = svd.singularValues().array().abs();
        for(int i=0; i<singularValues.size(); ++i)
            if(singularValues(i)<1)
                singularValues(i) = 1./singularValues(i);
        return singularValues.sum();
        
    };
    
    auto area_error = [](t_V3 p1, t_V3 p2, t_V3 p3, t_V3 q1, t_V3 q2, t_V3 q3) {
        const t_V_s a0 = (p2-p1).cross(p3-p1).norm();
        const t_V_s a1 = (q2-q1).cross(q3-q1).norm();
        return fabs(a0-a1)/a0;
    };*/
    
    error = t_error(F.rows());
    for(int face=0; face<F.rows(); ++face) {
        const t_V3& p0 = V.row(F(face,0));
        const t_V3& p1 = V.row(F(face,1));
        const t_V3& p2 = V.row(F(face,2));
        const t_V3& q0 = retV.row(retF(face,0));
        const t_V3& q1 = retV.row(retF(face,1));
        const t_V3& q2 = retV.row(retF(face,2));
        
        //Errors
        t_error_s conformalError = quasi_conformal_dist(p0, p1, p2, q0, q1, q2);
        //t_error_s svdError = svd_error(p0, p1, p2, q0, q1, q2);
        //t_error_s areaError = area_error(p0, p1, p2, q0, q1, q2);
        error(face) = conformalError;
    }
    
    return retVal;
    
}

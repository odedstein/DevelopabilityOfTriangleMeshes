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


#undef PARALLEL_COMPUTATION
#define IGNORE_VALENCE_3
#undef REQUIRE_MATLAB

#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <thread>

#include "types.h"
#include "meshstruct.h"
#include "timestepstruct.h"
#include "handle_configfile.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <viewer/OViewer.h>
#include <igl/viewer/Viewer.h>

#include <developableflow/energy_selector.h>
#include <developableflow/curvature_energy.h>
#include <developableflow/timestep.h>
#include <developableflow/mesh_postprocessing.h>
#include <developableflow/measure_once_cut_twice.h>

#include <tools/triangle_cosdihedral_angle.h>
#include <tools/qc_color.h>
#include <tools/write_energy.h>
#include <tools/write_cut_meshes.h>
#include <tools/numerical_gradient.h>
#include <tools/perturb.h>

#include <igl/read_triangle_mesh.h>
#include <igl/upsample.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/hausdorff.h>
#include <igl/parula.h>
#include <igl/colormap.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>


Developables::Mesh m; //Mesh struct
Timestep t; //Timestep struct
std::vector<std::string> titles; //Plot titles
std::vector<Eigen::VectorXd> lines; //Plot lines

std::vector<std::pair<int, int> > protocol; //Protocol

bool animating = false;
bool meshPosChanged = true;
bool meshChanged = false;
bool recordingVideo = false;

//Video variables
double lastEnergy = -1;
int frameCount = 0;

int energyColoring = 0;
bool dihedralColoring = true;
bool remeshingEnabled = true;
int renderMode = 0;
Linesearch linesearchMode = LINESEARCH_NONE;
StepType stepType = STEP_TYPE_GRADDESC;
EnergyType energyMode = ENERGY_TYPE_HINGE;
Scalar cutThreshold = 0.001; //Vertices with more energy than that are cut


void animation_step()
{
    //Perform a timestep
    const auto oldt = t; const auto oldm = m;
    int success = timestep(m.V, m.F, m.VF, m.VFi, m.isB, t.t, t.p, t.energy, t.energyGrad, linesearchMode, stepType, energyMode);
    if(success<0) {
        std::stringstream stream;
        //stream << "energy_at_step_" << t.totalSteps << ".mat";
        //write_energy(oldm.V, oldm.F, oldm.VF, oldm.VFi, oldm.isB, oldt.p, stream.str(), 1e-7, energyMode);
        std::cout << "Line search failed with error code " << success << " at step " << t.totalSteps << std::endl;
        //animating = false;
    }
    t.post_step_processing();
    
    //Do postprocessing
    if(remeshingEnabled) {
        int change = mesh_postprocessing(m.V, m.F, m.E, m.edgesC, m.TT, m.TTi, m.VF, false);
        if(change==1) { //Structural change happened
            std::cout << "A structural change happened to the mesh." << std::endl;
            m.origF = m.F;
            OVectorXi _1;
            igl::remove_unreferenced(OMatrixXs(m.V), OMatrixXi(m.F), m.V, m.F, _1);
            igl::remove_unreferenced(OMatrixXs(m.origV), OMatrixXi(m.origF), m.origV, m.origF, _1);
            m.update();
            t.invalidate();
            meshChanged = true;
        }
    }
    meshPosChanged = true;
    
    //Update labels
    double totalEnergy = toDouble(t.energy.sum());
    lines[0].conservativeResize(lines[0].size()+1);
    lines[0](lines[0].size()-1) = totalEnergy;
    std::stringstream title0;
    title0 << "Energy (" << totalEnergy << ")";
    titles[0] = title0.str();
    double energyGradL2 = t.energyGrad.rows()==0 ? 0 : doublecast(t.energyGrad).norm()/sqrt((double) t.energyGrad.rows());
    lines[1].conservativeResize(lines[1].size()+1);
    lines[1](lines[1].size()-1) = energyGradL2;
    std::stringstream title1;
    title1 << "Scaled L2 gradient norm (" << energyGradL2 << ")";
    titles[1] = title1.str();
    std::stringstream title2;
    title2 << "Timestep: " << toDouble(t.t) << ", total time: " << toDouble(t.totalT) << ", total steps: " << t.totalSteps;
    titles[2] = title2.str();
    
    //Record video
    if(recordingVideo) {
        if(lastEnergy < 0 || totalEnergy/lastEnergy < 0.99) {
            if(totalEnergy>0)
                lastEnergy = totalEnergy;
            std::stringstream str;
            str << "./frame_" << std::setw(8) << std::setfill('0') << frameCount++ << ".obj";
            igl::writeOBJ(str.str(), doublecast(m.V), intcast(m.F));
        }
    }
}

void predraw_callback(igl::OViewer& viewer, bool modHappened)
{
    //Execute instructions
    execute_instructions(t.totalSteps);
    
    if(animating)
        animation_step();
    
    if(meshChanged) {
        
        Eigen::MatrixXd V = doublecast(m.V);
        Eigen::MatrixXi F = intcast(m.F);
        Eigen::MatrixXi E = intcast(m.E);
        
        //viewer.modify_V(V);
        //viewer.modify_F(F, E);
        viewer.set_mesh(V, F, E);
        
        meshChanged = false;
    }
    
    if(meshPosChanged) {
        
        Eigen::MatrixXd V = doublecast(m.V);
        Eigen::MatrixXi F = intcast(m.F);
        static Eigen::MatrixXd origV = doublecast(m.origV);
        static Eigen::MatrixXi origF = intcast(m.origF);
        
        //Update mesh
        viewer.modify_V(V);
        
        //Dihedral angle coloring
        Eigen::MatrixXd normals;
        igl::per_face_normals(V, F, normals);
        if(dihedralColoring) {
            Eigen::MatrixXd dihedralAngles;
            triangle_cosdihedral_angle(intcast(m.edgesC), intcast(m.TT), normals, dihedralAngles);
            viewer.set_linetransparency(5.*dihedralAngles.array().acos()/M_PI);
        } else {
            viewer.set_linetransparency(Eigen::VectorXd::Ones(m.E.rows()));
        }
        
        //Color vertices by energy, threshold, or gradient dotted with normal
        if(energyColoring == 0 || t.energy.rows()==0) {
            viewer.set_colors(Eigen::Vector3d(0., 118./255., 189./255.));
        } else if(energyColoring == 1  && t.energy.rows() > 0) {
            Eigen::VectorXd energy = doublecast(t.energy);
            Eigen::VectorXd vertcolors(3*F.rows());
            for(int i=0; i<F.rows(); ++i) {
                for(int j=0; j<3; ++j) {
                    vertcolors(3*i+j) = energy(F(i,j));
                }
            }
            
            Eigen::MatrixXd colors;
            igl::parula(vertcolors, true, colors);
            viewer.set_colors(colors);
        } else if(energyColoring == 2 && t.energyGrad.rows() > 0) {
            Eigen::MatrixXd energyGrad = doublecast(t.energyGrad);
            Eigen::MatrixXd perVertexNormals;
            igl::per_vertex_normals(V, F, normals, perVertexNormals);
            Eigen::VectorXd vertcolors(3*F.rows());
            for(int i=0; i<F.rows(); ++i) {
                for(int j=0; j<3; ++j) {
                    vertcolors(3*i+j) = perVertexNormals.row(F(i,j)).dot(energyGrad.row(F(i,j)));
                }
            }
            
            Eigen::MatrixXd colors;
            igl::parula(vertcolors, true, colors);
            viewer.set_colors(colors);
        } else if(energyColoring == 3  && t.energy.rows() > 0) {
            Eigen::VectorXd energy = doublecast(t.energy);
            Eigen::MatrixXd colors(3*F.rows(), 3);
            for(int i=0; i<F.rows(); ++i) {
                for(int j=0; j<3; ++j) {
                    if(energy(F(i,j)) > cutThreshold)
                        colors.row(3*i+j) << 0., 118./255., 189./255.;
                    else
                        colors.row(3*i+j) << 0.8, 0.8, 0.8;
                }
            }
            viewer.set_colors(colors);
        } else if(energyColoring == 4) {
            Eigen::VectorXd gaussEnergy;
            curvature_energy(m.V, m.F, m.VF, m.VFi, m.isB, gaussEnergy);
            std::cout << "Total angle defect: " << gaussEnergy.sum() << std::endl;
            
            Eigen::VectorXd vertcolors(3*F.rows());
            for(int i=0; i<F.rows(); ++i) {
                for(int j=0; j<3; ++j) {
                    vertcolors(3*i+j) = gaussEnergy(F(i,j));
                }
            }
            Eigen::MatrixXd colors;
            igl::parula(vertcolors, true, colors);
            viewer.set_colors(colors);
        } else if(energyColoring == 5) {
            Eigen::MatrixXd N;
            igl::per_face_normals(m.V, m.F, N);
            Eigen::VectorXd A;
            igl::doublearea(m.V, m.F, A);
            Eigen::VectorXd err(N.rows());
            for(int face=0; face<N.rows(); ++face) {
                double weightSum = 0;
                Eigen::RowVector3d normalSum(0,0,0);
                for(int j=0; j<3; ++j) {
                    int neigh = m.TT(face,j);
                    if(neigh>=0) {
                        weightSum += A(neigh);
                        normalSum += A(neigh)*N.row(neigh);
                    }
                }
                if(weightSum>0)
                    err(face) = (normalSum/weightSum - N.row(face)).norm();
                else
                    err(face) = -1;
            }
            Eigen::MatrixXd errorColor;
            igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, err, 0, err.maxCoeff(), errorColor);
            viewer.set_colors(errorColor);
        }
        
        //Statistics about angles
        Eigen::VectorXd doubleareas;
        igl::doublearea(V, F, doubleareas);
        std::stringstream title3;
        title3 << "Smallest triangle area: " << doubleareas.minCoeff();
        Eigen::MatrixXd tipangles;
        igl::internal_angles(V, F, tipangles);
        title3 << ", smallest triangle angle: " << tipangles.minCoeff();
        titles[3] = title3.str();
        
        //Hausdorff distance
        double dist;
        igl::hausdorff(V, F, origV, origF, dist);
        std::stringstream title4;
        title4 << "Hausdorff dist to original: " << dist;
        titles[4] = title4.str();
        
        
        //Update plot
        viewer.set_plots(titles, lines);
        
        meshPosChanged = false;
    }
}

void keyboard_callback(igl::OViewer& viewer, int key, int scancode, int action, int modifier)
{
    if(key<0)
        exit(-1);
    
    if(action == GLFW_PRESS)
        protocol.push_back(std::make_pair(t.totalSteps, key));
    
    if(action == GLFW_PRESS && key == GLFW_KEY_W) {
        //Save mesh
        igl::writeOBJ("./out.obj", doublecast(m.V), intcast(m.F));
        //Output protocol to stdout
        std::cout << "# Protocol of actions (translations of numbers to keys: http://www.glfw.org/docs/latest/group__keys.html)" << std::endl;
        for(int p=0; p<protocol.size(); ++p)
            std::cout << "#" << protocol[p].first << " " << protocol[p].second << std::endl;
        std::cout << "#" << t.totalSteps << " " << -1 << std::endl;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_V) {
        //Toggle video recording
        recordingVideo = !recordingVideo;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_P) {
        //Toggle plotting
        viewer.set_plotting_enabled(!viewer.get_plotting_enabled());
    } else if(action == GLFW_PRESS && key == GLFW_KEY_L) {
        //One pass of Laplacian smoothing
        Eigen::SparseMatrix<Scalar> L, M;
        igl::cotmatrix(m.V, m.F, L);
        igl::massmatrix(m.V, m.F, igl::MASSMATRIX_TYPE_DEFAULT, M);
        Eigen::SparseMatrix<Scalar> lhs = M - 1e-5/igl::bounding_box_diagonal(m.V)*L;
        OMatrixXs rhs = M*m.V;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar> > solver(lhs);
        m.V = solver.solve(rhs);
        
        meshChanged = true;
        std::cout << "brap" << std::endl;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_D) {
        //Toggle dihedral angle coloring
        dihedralColoring = !dihedralColoring;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_S) {
        //In-place subdivision of mesh
        igl::upsample(OMatrixXs(m.V), OMatrixXi(m.F), m.V, m.F);
        igl::upsample(OMatrixXs(m.origV), OMatrixXi(m.origF), m.origV, m.origF);
        m.update();
        t.invalidate();
        meshChanged = true;
        meshPosChanged = true;
        //This invalidates the energy
        lastEnergy = -1;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_F) {
        //Cut and flatten mesh
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXd error;
        OVectorXi cut;
        measure_once_cut_twice(m.V, m.F, m.E, m.edgesC, m.TT, m.TTi, m.VF, m.VFi, m.isB, cutThreshold, cut, V, F, error);
        std::cout << "Flattening error, max: " << error.maxCoeff() << ", avg: " << error.mean() << std::endl;
        Eigen::VectorXd linet = Eigen::VectorXd::Zero(m.E.rows());
        for(int i=0; i<cut.size(); ++i)
            linet(cut(i)) = 1;
        viewer.set_linetransparency(linet);
        Eigen::MatrixXd errorColor;
        qc_color_matrix(error, errorColor);
        //igl::parula(error, true, errorColor);
        viewer.set_colors(errorColor);
        
        write_cut_meshes(m.V, m.F, V, F, m.origV, m.origF, "./flattened_out.obj", "./mesh_with_flattened_texcoord.obj", protocol, t.totalSteps);
    } else if(action == GLFW_PRESS && key == GLFW_KEY_C) {
        //Toggle gradient coloring and energy coloring (0 none, 1 energy, 2 gradient, 3 threshold, 4 angle defect, 6 smoothness)
        energyColoring = (energyColoring+1)%6;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_R) {
        //Toggle render mode (0 standard, 1 metal, 2 cel)
        renderMode = (renderMode+1)%3;
        if(renderMode==0) {
            viewer.set_colors(Eigen::Vector3d(0., 118./255., 189./255.));
            viewer.set_display_mode(igl::OViewer::FillWireframe);
            viewer.set_shininess(35.);
            viewer.set_lighting_factor(0.8);
            viewer.set_cel_shading(false);
        } else if(renderMode==1) {
            Eigen::MatrixXd c(m.F.rows(), 3);
            for(int i=0; i<m.F.rows(); ++i)
                c.row(i) = Eigen::RowVector3d(0.85, 0.8, 0.8);
            viewer.set_colors(0.1*c, c, 0.7*c);
            viewer.set_shininess(100.);
            viewer.set_lighting_factor(1.);
            viewer.set_display_mode(igl::OViewer::Fill);
            viewer.set_cel_shading(false);
        } else if(renderMode==2) {
            viewer.set_colors(Eigen::Vector3d(0., 118./255., 189./255.));
            viewer.set_display_mode(igl::OViewer::FillWireframe);
            viewer.set_shininess(35.);
            viewer.set_lighting_factor(0.8);
            viewer.set_cel_shading(true);
        }
    } else if(action == GLFW_PRESS && key == GLFW_KEY_M) {
        //Toggle linesearch mode
        linesearchMode = (Linesearch) ( (linesearchMode+1)%LINESEARCH_NUMS );
        if(linesearchMode==LINESEARCH_NONE)
            stepType = STEP_TYPE_GRADDESC;
        else if(linesearchMode==LINESEARCH_OVERTON)
            //stepType = STEP_TYPE_GRADDESC;
            stepType = STEP_TYPE_LBFGS;
        else if(linesearchMode==LINESEARCH_BACKTRACK)
            stepType = STEP_TYPE_GRADDESC;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_E) {
        //Toggle energy
        energyMode = (EnergyType) ( (energyMode+1)%ENERGY_TYPE_NUMS );
        t.invalidate();
        lastEnergy = -1;
        std::cout << "Using energy: " << energyMode << std::endl;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_TAB) {
        //Perform one animation step
        animation_step();
    } else if(action == GLFW_PRESS && key == GLFW_KEY_SPACE) {
        //Animate
        animating = !animating;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_UP) {
        t.t = exp(log(t.t)+1.);
        std::stringstream title2;
        title2 << "Timestep: " << toDouble(t.t) << ", total time: " << toDouble(t.totalT) << ", total steps: " << t.totalSteps;
        titles[2] = title2.str();
        viewer.set_plottitles(titles);
    } else if(action == GLFW_PRESS && key == GLFW_KEY_DOWN) {
        t.t = exp(log(t.t)-1.);
        std::stringstream title2;
        title2 << "Timestep: " << toDouble(t.t) << ", total time: " << toDouble(t.totalT) << ", total steps: " << t.totalSteps;
        titles[2] = title2.str();
        viewer.set_plottitles(titles);
    } else if(action == GLFW_PRESS && key == GLFW_KEY_N) {
        remeshingEnabled = !remeshingEnabled;
        std::cout << "remeshing enabled: " << remeshingEnabled << std::endl;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_Q) {
        //Plot the energy (with resolution as function argument) along the search direction
        std::stringstream stream;
        stream << "energy_at_step_" << t.totalSteps << ".mat";
        write_energy(m.V, m.F, m.VF, m.VFi, m.isB, t.p, stream.str(), 1e-8, energyMode);
    } else if(action == GLFW_PRESS && key == GLFW_KEY_G) {
        //Numerical gradient test
        const std::function<Scalar(const Eigen::PlainObjectBase<OMatrixXs>&)> energyFct = [] (const Eigen::PlainObjectBase<OMatrixXs>& pos)->Scalar {
            OVectorXs energy;
            OMatrixXs dummy;
            energy_selector(energyMode, pos, m.F, m.VF, m.VFi, m.isB, energy, dummy);
            return energy.sum();
        };
        OMatrixXs numericalGrad;
        numerical_gradient(energyFct, m.V, 1e-8, numericalGrad);
        OMatrixXs dummy, energyGrad;
        energy_selector(energyMode, m.V, m.F, m.VF, m.VFi, m.isB, dummy, energyGrad);
        
        std::cout << "The difference between exact and numerical grad is: " << (numericalGrad - energyGrad).norm() << ", numerical norm: " << numericalGrad.norm() << ", exact norm: " << energyGrad.norm() << std::endl;
    } else if(action == GLFW_PRESS && key == GLFW_KEY_J) {
        perturb(m.V, m.VF, m.isB, 0.01*igl::bounding_box_diagonal(m.V));
        meshPosChanged = true;
    }
}


int main(int argc, char *argv[])
{
    //Init viewer
    igl::OViewer viewer;
    
    //Initialize the scalar data type
    types_startup();
    Eigen::initParallel();
    
    
    //Read in mesh
    std::string path;
    if(argc>1) {
        path = std::string(argv[1]);
    } else {
        path = "../assets/bunny.obj";
    }
    if(path.compare(path.length()-3,3,"txt")==0) {
        //There is a config file
        read_configfile(path, &viewer, keyboard_callback, m);
    } else {
        OMatrixXs V;
        OMatrixXi F;
        igl::read_triangle_mesh(path, V, F);
        m = Developables::Mesh(V, F);
    }
    
    
    //Plots
    titles = {"Energy", "Scaled L2 gradient norm", "Timestep: , total time: , total steps: ", "Smallest triangle area: , smallest triangle angle: ", ""};
    lines = {Eigen::VectorXd(1), Eigen::VectorXd(1)};
    lines[0] << 0.;
    lines[1] << 1.;
    
    
    //viewer
    viewer.set_mesh(doublecast(m.V), intcast(m.F), intcast(m.E));
    viewer.set_plotting_enabled(true);
    viewer.set_plots(titles, lines);
    viewer.set_face_based_coloring(true);
    viewer.predraw_callback = predraw_callback;
    viewer.keyboard_callback = keyboard_callback;
    
    viewer.launch();
    
    return EXIT_SUCCESS;
}

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


#include "OViewer.h"

#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <iostream>

#include <igl/opengl/gl.h>
#include <igl/quat_to_mat.h>
#include <igl/barycenter.h>
#include <igl/frustum.h>
#include <igl/ortho.h>
#include <igl/edges.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/look_at.h>
#include <igl/trackball.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>

#include "shaders.txt"


namespace igl {
    
    //CLASS INITIALIZATION
    
    IGL_INLINE void OViewer::defaults()
    {
        backgroundColor = Eigen::Vector3f(0.9, 0.9, 0.9);
        wireColor = Eigen::Vector3f(0., 0., 0.);
        defaultColor = Eigen::Vector3d(0., 118./255., 189./255.);
        set_colors(defaultColor);
        defaultColorUsed = true;
        m_linetransp = Eigen::VectorXf::Constant(m_E_mod.rows(), 1.);
        
        pointsOutlineColor = Eigen::Vector3f(0., 0., 0.);
        pointsFillColor = Eigen::Vector3f(255./255., 255./255., 179./255.);
        pointsImportantColor = Eigen::Vector3f(141./255., 211./255., 199./255.);
        pointsSelectedColor = Eigen::Vector3f(225./255., 26./255., 39./255.);
        
        //cameraRotation = Eigen::Quaternion<float>::Identity();
        cameraRotation = Eigen::AngleAxis<float>(-0.1*M_PI, Eigen::Vector3f(1., 0., 0.));
        cameraTranslation = Eigen::Vector3f::Zero();
        
        camera_eye = Eigen::Vector3f(0., 0., 5.);
        camera_center = Eigen::Vector3f(0., 0., 0.);
        camera_up = Eigen::Vector3f(0., 1., 0.);
        
        lightPosition = Eigen::Vector3f(0.0, 0.3, 5.0);
        
        plotColors = {Eigen::Vector3f(31./255., 120./255., 180./255.),
                        Eigen::Vector3f(51./255., 160./255., 44./255.),
                        Eigen::Vector3f(227./255., 26./255., 28./255.),
                        Eigen::Vector3f(255./255., 127./255., 0./255.),
                        Eigen::Vector3f(106./255., 61./255., 154./255.),
                        Eigen::Vector3f(177./255., 89./255., 40./255.),
                        Eigen::Vector3f(166./255., 206./255., 227./255.),
                        Eigen::Vector3f(178./255., 223./255., 138./255.),
                        Eigen::Vector3f(251./255., 154./255., 153./255.),
                        Eigen::Vector3f(253./255., 191./255., 111./255.),
                        Eigen::Vector3f(202./255., 178./255., 214./255.),
                        Eigen::Vector3f(255./255., 255./255., 153./255.)};
        plotBGColor = Eigen::Vector3f(0.97, 0.97, 0.97);
    }
    
    IGL_INLINE OViewer::OViewer()
    {
        defaults();
    }
    
    IGL_INLINE OViewer::OViewer(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
    {
        set_mesh(V, F);
        defaults();
    }
    
    IGL_INLINE OViewer::OViewer(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& E)
    {
        set_mesh(V, F, E);
        defaults();
    }
    
    
    
    //GETTERS & SETTERS
    
    #include "OViewer_getters_setters.cxx"
    
    
    
    //GRAPHICS
    
    #include "OViewer_graphics.cxx"
    
}

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


#ifndef DEVELOPABLEFLOW_HANDLE_CONFIGFILE_H
#define DEVELOPABLEFLOW_HANDLE_CONFIGFILE_H

#include <iostream>
#include <fstream>

#include <viewer/OViewer.h>

#include <igl/read_triangle_mesh.h>


std::vector<std::pair<int, int> > instructions; //Instructions
igl::OViewer *viewerref;
void (*keyboard_fct)(igl::OViewer&, int, int, int, int);


void read_configfile(const std::string& path, igl::OViewer *viewer, void (*kfct)(igl::OViewer&, int, int, int, int), Developables::Mesh &m) {
    
    std::ifstream input(path);
    std::string meshpath;
    std::getline(input, meshpath);
    
    OMatrixXs V;
    OMatrixXi F;
    igl::read_triangle_mesh(meshpath, V, F);
    m = Developables::Mesh(V, F);
    
    int timestep, keypress;
    while(input >> timestep >> keypress) {
        instructions.push_back(std::make_pair(timestep,keypress));
    }
    
    std::reverse(instructions.begin(), instructions.end());
    
    keyboard_fct = kfct;
    viewerref = viewer;
}


void execute_instructions(const int& step) {
    if(instructions.size()>0) {
        while(instructions.back().first==step) {
            std::cout << "Executing instruction " << instructions.back().first << " " << instructions.back().second << std::endl;
            keyboard_fct(*viewerref, instructions.back().second, -1, GLFW_PRESS, -1);
            instructions.pop_back();
        }
    }
}


#endif

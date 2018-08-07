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


#ifndef DEVELOPABLEFLOW_TIMESTEPSTRUCT_H
#define DEVELOPABLEFLOW_TIMESTEPSTRUCT_H


struct Timestep {
    Scalar t = 1e-5; //0.1; //Timestep
    Scalar totalT = 0; //Total time
    int totalSteps = 0;
    OVectorXs energy; //Cached energy
    OMatrixXs energyGrad; //Cached grad
    OMatrixXs p; //Search direction
    
    void post_step_processing() //Increases all relevant counters
    {
        totalT += t;
        ++totalSteps;
    }
    
    void invalidate()
    {
        energy = OVectorXs();
        energyGrad = OMatrixXs();
        p = OMatrixXs();
        t = 1e-8;
    }
};


#endif

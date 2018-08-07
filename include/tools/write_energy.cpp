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


#include "write_energy.h"

#include <developableflow/energy_selector.h>

#ifdef REQUIRE_MATLAB
#include <igl/matlab/MatlabWorkspace.h>
#endif

#include <limits>
#include <deque>

#include <Eigen/Core>


template <typename derivedV, typename derivedF, typename indexType, typename derivedP, typename t_res>
IGL_INLINE void write_energy(
                             const Eigen::PlainObjectBase<derivedV>& V,
                             const Eigen::PlainObjectBase<derivedF>& F,
                             const std::vector<std::vector<indexType> >& VF,
                             const std::vector<std::vector<indexType> >& VFi,
                             const std::vector<bool>& isB,
                             const Eigen::PlainObjectBase<derivedP>& p,
                             const std::string& filename,
                             const t_res& resolution,
                             EnergyType energyType)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, Eigen::Dynamic> t_Vdyn;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef typename derivedP::Scalar t_p_s;
    typedef Eigen::Matrix<t_p_s, Eigen::Dynamic, 3> t_p;
    
    assert(p.rows() == V.rows() && "The search direction must be a valid step direction for your current mesh.");
    
#ifdef REQUIRE_MATLAB
    int nSamples = 10000;
    
    t_Vv alphas(nSamples);
    t_Vv energies(nSamples);
    t_Vv energyGrads(nSamples);
    t_Vv energyGradDotPs(nSamples);
    
    t_Vdyn energiesPerVert(nSamples, V.rows());
    t_Vdyn energyGradsPerVert(nSamples, 3*V.rows());
    t_Vdyn smallestEigenvecPerVert(nSamples, 3*V.rows());
    
    for(int i=0; i<nSamples; ++i) {
        t_V_s alpha = resolution*(i-0.5*nSamples);
        t_Vv en;
        t_V egrad;
        
        t_V newV = V + alpha*p;
        
        energy_selector(energyType, newV, F, VF, VFi, isB, en, egrad);
        
        alphas(i) = alpha;
        energies(i) = en.sum();
        energyGrads(i) = egrad.norm() / sqrt((t_V_s) egrad.rows());
        energyGradDotPs(i) = (egrad.array()*p.array()).sum();
        
        energiesPerVert.row(i) = en;
        t_Vv reshapedGrads = Eigen::Map<t_Vv>(egrad.data(), egrad.size());
        energyGradsPerVert.row(i) = reshapedGrads;
        
        t_V minCurvatureVecs;
        hinge_energy(newV, F, VF, VFi, isB, en, minCurvatureVecs);
        t_Vv reshapedCurvs = Eigen::Map<t_Vv>(minCurvatureVecs.data(), minCurvatureVecs.size());
        smallestEigenvecPerVert.row(i) = reshapedCurvs;
    }
    
    igl::matlab::MatlabWorkspace mw;
    mw.save(alphas, "alphas");
    mw.save(energies, "energies");
    mw.save(energyGrads, "energyGrads");
    mw.save(energyGradDotPs, "energyGradDotPs");
    mw.save(energiesPerVert, "energiesPerVert");
    mw.save(energyGradsPerVert, "energyGradsPerVert");
    mw.save(smallestEigenvecPerVert, "smallestEigenvecPerVert");
    mw.save(doublecast(V), "V");
    mw.save_index(intcast(F), "F");
    mw.write(filename);
#else
    std::cerr << "Please enable MATLAB for this feature." << std::endl;
#endif
}

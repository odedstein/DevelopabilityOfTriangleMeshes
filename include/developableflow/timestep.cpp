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


#include "timestep.h"

#include <limits>
#include <deque>

#include <Eigen/Core>

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>

//A Few compiler flags that control the line searches
#define MAX_LINESEARCH_TRIES 100
#define ARMIJO_C1 1e-4
#define WOLFE_C2 0.99
#define N_LBFGS_VECTORS 8
#define MAX_T 100. //1.
#define MIN_T 1e-12

#undef CLAMP
#undef NORMALIZE
#undef CONSTRAIN_TO_CYLINDER
#undef NORMAL_DAMPENING

//Infinity
#define INFTY std::numeric_limits<double>::infinity()


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad, typename derivedT, typename derivedP>
IGL_INLINE int hinge_timestep(
                              Eigen::PlainObjectBase<derivedV>& V,
                              const Eigen::PlainObjectBase<derivedF>& F,
                              const std::vector<std::vector<indexType> >& VF,
                              const std::vector<std::vector<indexType> >& VFi,
                              const std::vector<bool>& isB,
                              derivedT& t,
                              Eigen::PlainObjectBase<derivedP>& p,
                              Eigen::PlainObjectBase<derivedEnergy>& energy,
                              Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad,
                              Linesearch mode,
                              StepType type)
{
    return timestep_tool(V, F, VF, VFi, isB, t, p, energy, energyGrad, mode, type, ENERGY_TYPE_HINGE);
}

template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad, typename derivedT, typename derivedP>
IGL_INLINE int minwidth_timestep(
                                 Eigen::PlainObjectBase<derivedV>& V,
                                 const Eigen::PlainObjectBase<derivedF>& F,
                                 const std::vector<std::vector<indexType> >& VF,
                                 const std::vector<std::vector<indexType> >& VFi,
                                 const std::vector<bool>& isB,
                                 derivedT& t,
                                 Eigen::PlainObjectBase<derivedP>& p,
                                 Eigen::PlainObjectBase<derivedEnergy>& energy,
                                 Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad,
                                 Linesearch mode,
                                 StepType type)
{
    return timestep_tool(V, F, VF, VFi, isB, t, p, energy, energyGrad, mode, type, ENERGY_TYPE_MINWIDTH);
}


template <typename derivedV, typename derivedF, typename indexType, typename derivedEnergy, typename derivedEnergyGrad, typename derivedT, typename derivedP>
IGL_INLINE int timestep(
                        Eigen::PlainObjectBase<derivedV>& V,
                        const Eigen::PlainObjectBase<derivedF>& F,
                        const std::vector<std::vector<indexType> >& VF,
                        const std::vector<std::vector<indexType> >& VFi,
                        const std::vector<bool>& isB,
                        derivedT& t,
                        Eigen::PlainObjectBase<derivedP>& p,
                        Eigen::PlainObjectBase<derivedEnergy>& energy,
                        Eigen::PlainObjectBase<derivedEnergyGrad>& energyGrad,
                        Linesearch mode,
                        StepType type,
                        EnergyType energyType)
{
    typedef typename derivedV::Scalar t_V_s;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 3> t_V;
    typedef Eigen::Matrix<t_V_s, Eigen::Dynamic, 1> t_Vv;
    typedef Eigen::Matrix<t_V_s, 1, 3> t_V3;
    typedef typename derivedF::Scalar t_F_i;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 3> t_F;
    typedef Eigen::Matrix<t_F_i, Eigen::Dynamic, 1> t_Fv;
    typedef typename derivedP::Scalar t_p_s;
    typedef Eigen::Matrix<t_p_s, Eigen::Dynamic, 3> t_p;
    typedef Eigen::Matrix<t_p_s, 1, 3> t_p3;
    typedef typename derivedEnergy::Scalar t_energy_s;
    typedef Eigen::Matrix<t_energy_s, Eigen::Dynamic, 1> t_energy;
    typedef typename derivedEnergyGrad::Scalar t_energyGrad_s;
    typedef Eigen::Matrix<t_energyGrad_s, Eigen::Dynamic, 3> t_energyGrad;
    typedef Eigen::Matrix<t_energyGrad_s, 1, 3> t_energyGrad3;
    
    
    //The matrices-are-actually-vectors dot product. Will be used a lot.
    const auto dot = [] (t_p A, t_energyGrad B) { return (A.array()*B.array()).sum(); };
    
    int retVal = 0;
    
    bool resetHappened = false;
    if(p.rows()<1) {
        //The search direction is not initialized or was just reset. Initialize it with the gradient.
        energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
        p = -energyGrad;
        resetHappened = true;
    }
    
    //Cache old values
    t_V oldV = V;
    t_energyGrad oldGrad = energyGrad;
    t_energy_s oldTotalEnergy = energy.sum();
    
    if(mode==LINESEARCH_NONE) {
        //No linesearch happening. Just step towards the direction
        V += t*p;
#ifdef CONSTRAIN_TO_CYLINDER
        for(int i=0; i<V.rows(); ++i)
            V.block(i,0,1,2).normalize();
#endif
        
        //Calculate energy and gradient at the new spot
        energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
        
    } else if(mode==LINESEARCH_OVERTON) {
        //Try to cautiously increase t
        t = std::min(2.*t, MAX_T);
        
        t_p_s g0dotp = dot(p, oldGrad);
        derivedT tmin = 0, tmax = INFTY;
        int tries = 0;
        for(; tries<MAX_LINESEARCH_TRIES; ++tries) {
            V = oldV + t*p;
#ifdef CONSTRAIN_TO_CYLINDER
            for(int i=0; i<V.rows(); ++i)
                V.block(i,0,1,2).normalize();
#endif
            energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
            t_energy_s totalEnergy = energy.sum();
            
            if(totalEnergy > oldTotalEnergy + ARMIJO_C1*t*g0dotp) {
                tmax = t;
            } else if(dot(p, energyGrad) < WOLFE_C2*g0dotp) {
                tmin = t;
            } else {
                break;
            }
            
            if(tmax < INFTY) {
                t = 0.5*(tmin+tmax);
            } else {
                t = 2*tmin;
            }
        }
        
        if(tries==MAX_LINESEARCH_TRIES)
            retVal = -1;
        
        if(t > MAX_T) {
            t = MAX_T;
            V = oldV + t*p;
#ifdef CONSTRAIN_TO_CYLINDER
            for(int i=0; i<V.rows(); ++i)
                V.block(i,0,1,2).normalize();
#endif
            energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
        }
        if(t < MIN_T) {
            t = MIN_T;
            V = oldV + t*p;
#ifdef CONSTRAIN_TO_CYLINDER
            for(int i=0; i<V.rows(); ++i)
                V.block(i,0,1,2).normalize();
#endif
            energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
        }
    } else if(mode==LINESEARCH_BACKTRACK) {
        const t_V_s rho = 0.9;
        
        //Try to cautiously increase t
        t = std::min(2.*t, MAX_T);
        
        t_p_s g0dotp = dot(p, oldGrad);
        int tries = 0;
        for(; tries<MAX_LINESEARCH_TRIES; ++tries) {
            V = oldV + t*p;
#ifdef CONSTRAIN_TO_CYLINDER
            for(int i=0; i<V.rows(); ++i)
                V.block(i,0,1,2).normalize();
#endif
            energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
            t_energy_s totalEnergy = energy.sum();
            
            if(totalEnergy > oldTotalEnergy + ARMIJO_C1*t*g0dotp) {
                t *= rho;
            } else {
                break;
            }
        }
        
        if(tries==MAX_LINESEARCH_TRIES) {
            retVal = -1;
        }
        
        if(t > MAX_T) {
            t = MAX_T;
            V = oldV + t*p;
#ifdef CONSTRAIN_TO_CYLINDER
            for(int i=0; i<V.rows(); ++i)
                V.block(i,0,1,2).normalize();
#endif
            energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
        }
        if(t < MIN_T) {
            t = MIN_T;
            V = oldV + t*p;
#ifdef CONSTRAIN_TO_CYLINDER
            for(int i=0; i<V.rows(); ++i)
                V.block(i,0,1,2).normalize();
#endif
            energy_selector(energyType, V, F, VF, VFi, isB, energy, energyGrad);
        }
        
    }
    
    
    if(type==STEP_TYPE_GRADDESC) {
        //For gradient descent, the new search direction is just minus the energy gradient
        p = -energyGrad;
    } else if(type==STEP_TYPE_LBFGS) {
        //Cache vectors
        static std::deque<t_energyGrad> y;
        static std::deque<t_V> s;
        if(resetHappened) {
            s.clear();
            y.clear();
        }
        y.push_back(energyGrad - oldGrad);
        s.push_back(V-oldV);
        if(s.size() > N_LBFGS_VECTORS) {
            y.pop_front();
            s.pop_front();
        }
        int nvecs = s.size();
        
        //Compute next step
        t_energyGrad q = energyGrad;
        t_Vv alpha(nvecs), rho(nvecs);
        for(int i = nvecs-1; i>=0; --i) {
            rho(i) = 1./dot(y[i],s[i]);
            if(rho(i) != rho(i) || !std::isfinite(rho(i)))
                rho(i) = 0.;
            alpha(i) = rho(i)*dot(s[i],q);
            q -= alpha(i)*y[i];
        }
        t_energyGrad_s gammak = dot(s[nvecs-1], y[nvecs-1])/dot(y[nvecs-1], y[nvecs-1]);
        q *= gammak;
        for(int i=0; i<nvecs; ++i) {
            t_energyGrad_s beta = rho(i)*dot(y[i], q);
            q += s[i]*(alpha(i)-beta);
        }
        p = -q;
        if(p!=p || retVal == -1) {
            s.clear();
            y.clear();
            p = t_p();
            resetHappened = true;
        }
    }
    
    
#ifdef CLAMP
    t_p_s maxNorm = 0.;
    for(int i=0; i<p.rows(); ++i) {
        const t_p3& g = p.row(i);
        if(g.norm() > maxNorm)
            maxNorm = g.norm();
    }
    
    for(int i=0; i<p.rows(); ++i) {
        const t_p3& g = p.row(i);
        if(g.norm() > 0.5*maxNorm)
            p.row(i) = g*0.5*maxNorm/g.norm();
        
    }
#endif
    
#ifdef NORMALIZE
    for(int i=0; i<p.rows(); ++i) {
        p.row(i).normalize();
    }
#endif
    
#ifdef NORMAL_DAMPENING
    t_V normals;
    igl::per_face_normals(V, F, normals);
    t_V perVertexNormals;
    igl::per_vertex_normals(t_V(V), t_F(F), normals, perVertexNormals);
    for(int i=0; i<p.rows(); ++i) {
        const t_V3& n = perVertexNormals.row(i);
        const t_V3 normpart = p.row(i).dot(n)*n;
        const t_V3 tangpart = p.row(i) - normpart;
        p.row(i) = tangpart + normpart/1000;
        //p.row(i) = tangpart*1000 + normpart;
    }
#endif
    
    if(V != V) {
        retVal = -2;
        assert(false && "There are nans in the vertices return value");
    }
    
    return retVal;
}

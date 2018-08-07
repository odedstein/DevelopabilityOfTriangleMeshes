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


#include "qc_color.h"

static Eigen::Vector3d qc_HSV(double h, double s, double v)
{
    double r = 0., g = 0., b = 0.;
    
    if (s == 0)
    {
        r = v;
        g = v;
        b = v;
    }
    else
    {
        h = (h == 1 ? 0 : h) * 6;
        
        int i = (int) floor(h);
        
        double f = h - i;
        double p = v * (1 - s);
        double q = v * (1 - (s * f));
        double t = v * (1 - s * (1 - f));
        
        switch (i)
        {
            case 0:
                r = v;
                g = t;
                b = p;
                break;
                
            case 1:
                r = q;
                g = v;
                b = p;
                break;
                
            case 2:
                r = p;
                g = v;
                b = t;
                break;
                
            case 3:
                r = p;
                g = q;
                b = v;
                break;
                
            case 4:
                r = t;
                g = p;
                b = v;
                break;
                
            case 5:
                r = v;
                g = p;
                b = q;
                break;
                
            default:
                break;
        }
    }
    
    return Eigen::Vector3d(r, g, b);
}


template <typename derivedError>
IGL_INLINE Eigen::Vector3d qc_color(const derivedError& error)
{
    double qc = std::max(1., std::min(1.5, error));
    return qc_HSV((2.-4.*(qc-1.))/3., 0.7, 0.65);
}


template <typename derivedError, typename derivedColor>
IGL_INLINE void qc_color_matrix(const Eigen::PlainObjectBase<derivedError>& errors,
                                Eigen::PlainObjectBase<derivedColor>& colors)
{
    typedef typename derivedError::Scalar t_error_s;
    typedef Eigen::Matrix<t_error_s, Eigen::Dynamic, 1> t_error;
    typedef typename derivedError::Scalar t_color_s;
    typedef Eigen::Matrix<t_color_s, Eigen::Dynamic, 3> t_colors;
    
    colors = t_colors(errors.size(), 3);
    for(int i=0; i<errors.size(); ++i) {
        colors.row(i) = qc_color(errors(i));
    }
}


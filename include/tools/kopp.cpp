// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "kopp.h"




//EASY ACCESS FUNCTION, TEMPLATED

template<typename derivedScalar>
IGL_INLINE void eigendecomp(const Eigen::PlainObjectBase<derivedScalar>& A,
                            Eigen::PlainObjectBase<derivedScalar>& D,
                            Eigen::PlainObjectBase<derivedScalar>& V)
{
    typedef typename derivedScalar::Scalar t_s;
    typedef Eigen::Matrix<t_s, 3, 3> t_mat;
    typedef Eigen::Matrix<t_s, 3, 1> t_vec;
    
    t_s M[3][3];
    t_s Q[3][3];
    t_s w[3];
    
    M[0][0] = A(0,0); M[0][1] = A(0,1); M[0][2] = A(0,2);
    /*M[1][0] = A(1,0);*/ M[1][1] = A(1,1); M[1][2] = A(1,2);
    /*M[2][0] = A(2,0); M[2][1] = A(2,1);*/ M[2][2] = A(2,2);
    
    dsyevh3<t_s>(M, Q, w);
    //dsyevq3<t_s>(M, Q, w);
    
    D = t_mat();
    D(0,0) = w[0]; D(0,1) =  0.0; D(0,2) =  0.0;
    D(1,0) =  0.0; D(1,1) = w[1]; D(1,2) =  0.0;
    D(2,0) =  0.0; D(2,1) =  0.0; D(2,2) = w[2];
    
    V = t_mat();
    V(0,0) = Q[0][0]; V(0,1) = Q[0][1]; V(0,2) = Q[0][2];
    V(1,0) = Q[1][0]; V(1,1) = Q[1][1]; V(1,2) = Q[1][2];
    V(2,0) = Q[2][0]; V(2,1) = Q[2][1]; V(2,2) = Q[2][2];
    
    // sort by eigenvalue magntiude, smallest first
    if( fabs(D(0,0)) > fabs(D(1,1))) {
        t_s dbuf = D(0,0);
        D(0,0) = D(1,1);
        D(1,1) = dbuf;
        t_vec vbuf = V.col(0);
        V.col(0) = V.col(1);
        V.col(1) = V.col(0);
    }
    if( fabs(D(0,0)) > fabs(D(2,2))) {
        t_s dbuf = D(0,0);
        D(0,0) = D(2,2);
        D(2,2) = dbuf;
        t_vec vbuf = V.col(0);
        V.col(0) = V.col(2);
        V.col(2) = V.col(0);
    }
    if( fabs(D(1,1)) > fabs(D(2,2))) {
        t_s dbuf = D(1,1);
        D(1,1) = D(2,2);
        D(2,2) = dbuf;
        t_vec vbuf = V.col(1);
        V.col(1) = V.col(2);
        V.col(2) = V.col(1);
    }
    
}




// Constants
//#define M_SQRT3    1.73205080756887729357372529559455642811371944844722747802734375l   // sqrt(3)

// Macros
#define SQR(x)      ((x)*(x))                        // x^2


// ----------------------------------------------------------------------------
template<typename Scalar>
IGL_INLINE int dsyevc3(Scalar A[3][3], Scalar w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
{
    Scalar m, c1, c0;
    
    // Determine coefficients of characteristic poynomial. We write
    //       | a   d   f  |
    //  A =  | d*  b   e  |
    //       | f*  e*  c  |
    Scalar de = A[0][1] * A[1][2];                                    // d * e
    Scalar dd = SQR(A[0][1]);                                         // d^2
    Scalar ee = SQR(A[1][2]);                                         // e^2
    Scalar ff = SQR(A[0][2]);                                         // f^2
    m  = A[0][0] + A[1][1] + A[2][2];
    c1 = (A[0][0]*A[1][1] + A[0][0]*A[2][2] + A[1][1]*A[2][2])        // a*b + a*c + b*c - d^2 - e^2 - f^2
    - (dd + ee + ff);
    c0 = A[2][2]*dd + A[0][0]*ee + A[1][1]*ff - A[0][0]*A[1][1]*A[2][2]
    - Scalar(2.0) * A[0][2]*de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
    
    Scalar p, sqrt_p, q, c, s, phi;
    p = SQR(m) - Scalar(3.0)*c1;
    q = m*(p - (Scalar(3.0)/Scalar(2.0))*c1) - (Scalar(27.0)/Scalar(2.0))*c0;
    sqrt_p = sqrt(fabs(p));
    
    phi = Scalar(27.0) * ( Scalar(0.25)*SQR(c1)*(p - c1) + c0*(q + Scalar(27.0)/Scalar(4.0)*c0));
    phi = (Scalar(1.0)/Scalar(3.0)) * atan2(sqrt(fabs(phi)), q);
    
    c = sqrt_p*cos(phi);
    s = (Scalar(1.0)/sqrt(Scalar(3.0)))*sqrt_p*sin(phi);
    
    w[1]  = (Scalar(1.0)/Scalar(3.0))*(m - c);
    w[2]  = w[1] + s;
    w[0]  = w[1] + c;
    w[1] -= s;
    
    return 0;
}


// ----------------------------------------------------------------------------
template<typename Scalar>
IGL_INLINE int dsyevh3(Scalar A[3][3], Scalar Q[3][3], Scalar w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using Cardano's method for the eigenvalues and an analytical
// method based on vector cross products for the eigenvectors. However,
// if conditions are such that a large error in the results is to be
// expected, the routine falls back to using the slower, but more
// accurate QL algorithm. Only the diagonal and upper triangular parts of A need
// to contain meaningful values. Access to A is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
// Dependencies:
//   dsyevc3(), dsytrd3(), dsyevq3()
// ----------------------------------------------------------------------------
// Version history:
//   v1.1: Simplified fallback condition --> speed-up
//   v1.0: First released version
// ----------------------------------------------------------------------------
{
#ifndef EVALS_ONLY
    Scalar norm;          // Squared norm or inverse norm of current eigenvector
    //  Scalar n0, n1;        // Norm of first and second columns of A
    Scalar error;         // Estimated maximum roundoff error
    Scalar t, u;          // Intermediate storage
    int j;                // Loop counter
#endif
    
    // Calculate eigenvalues
    dsyevc3<Scalar>(A, w);
    
#ifndef EVALS_ONLY
    //  n0 = SQR(A[0][0]) + SQR(A[0][1]) + SQR(A[0][2]);
    //  n1 = SQR(A[0][1]) + SQR(A[1][1]) + SQR(A[1][2]);
    
    t = fabs(w[0]);
    if ((u=fabs(w[1])) > t)
        t = u;
    if ((u=fabs(w[2])) > t)
        t = u;
    if (t < 1.0l)
        u = t;
    else
        u = SQR(t);
    error = Scalar(256.0) * std::numeric_limits<Scalar>::epsilon() * SQR(u);
    //  error = 256.0l * LDBL_EPSILON * (n0 + u) * (n1 + u);
    
    Q[0][1] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
    Q[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0];
    Q[2][1] = SQR(A[0][1]);
    
    // Calculate first eigenvector by the formula
    //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
    Q[0][0] = Q[0][1] + A[0][2]*w[0];
    Q[1][0] = Q[1][1] + A[1][2]*w[0];
    Q[2][0] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Q[2][1];
    norm    = SQR(Q[0][0]) + SQR(Q[1][0]) + SQR(Q[2][0]);
    
    // If vectors are nearly linearly dependent, or if there might have
    // been large cancellations in the calculation of A[i][i] - w[0], fall
    // back to QL algorithm
    // Note that this simultaneously ensures that multiple eigenvalues do
    // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
    // i.e. all columns of A - w[0] * I are linearly dependent.
    if (norm <= error)
        return dsyevq3<Scalar>(A, Q, w);
    else                      // This is the standard branch
    {
        norm = sqrt(Scalar(1.0) / norm);
        for (j=0; j < 3; j++)
            Q[j][0] = Q[j][0] * norm;
    }
    
    // Calculate second eigenvector by the formula
    //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
    Q[0][1]  = Q[0][1] + A[0][2]*w[1];
    Q[1][1]  = Q[1][1] + A[1][2]*w[1];
    Q[2][1]  = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Q[2][1];
    norm     = SQR(Q[0][1]) + SQR(Q[1][1]) + SQR(Q[2][1]);
    if (norm <= error)
        return dsyevq3<Scalar>(A, Q, w);
    else
    {
        norm = sqrt(Scalar(1.0) / norm);
        for (j=0; j < 3; j++)
            Q[j][1] = Q[j][1] * norm;
    }
    
    // Calculate third eigenvector according to
    //   v[2] = v[0] x v[1]
    Q[0][2] = Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1];
    Q[1][2] = Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1];
    Q[2][2] = Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1];
#endif
    
    return 0;
}


// ----------------------------------------------------------------------------
template<typename Scalar>
IGL_INLINE int dsyevq3(Scalar A[3][3], Scalar Q[3][3], Scalar w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   dsytrd3()
// ----------------------------------------------------------------------------
{
    const int n = 3;
    Scalar e[3];                   // The third element is used only as temporary workspace
    Scalar g, r, p, f, b, s, c, t; // Intermediate storage
    int nIter;
    int m;
    
    // Transform A to real tridiagonal form by the Householder method
    dsytrd3<Scalar>(A, Q, w, e);
    
    // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
    // with the QL method
    //
    // Loop over all off-diagonal elements
    for (int l=0; l < n-1; l++)
    {
        nIter = 0;
        while (1)
        {
            // Check for convergence and exit iteration loop if off-diagonal
            // element e(l) is zero
            for (m=l; m <= n-2; m++)
            {
                g = fabs(w[m])+fabs(w[m+1]);
                if (fabs(e[m]) + g == g)
                    break;
            }
            if (m == l)
                break;
            
            if (nIter++ >= 30)
                return -1;
            
            // Calculate g = d_m - k
            g = (w[l+1] - w[l]) / (e[l] + e[l]);
            r = sqrt(SQR(g) + Scalar(1.0));
            if (g > 0)
                g = w[m] - w[l] + e[l]/(g + r);
            else
                g = w[m] - w[l] + e[l]/(g - r);
            
            s = c = Scalar(1.0);
            p = Scalar(0.0);
            for (int i=m-1; i >= l; i--)
            {
                f = s * e[i];
                b = c * e[i];
                if (fabs(f) > fabs(g))
                {
                    c      = g / f;
                    r      = sqrt(SQR(c) + Scalar(1.0));
                    e[i+1] = f * r;
                    c     *= (s = Scalar(1.0)/r);
                }
                else
                {
                    s      = f / g;
                    r      = sqrt(SQR(s) + Scalar(1.0));
                    e[i+1] = g * r;
                    s     *= (c = Scalar(1.0)/r);
                }
                
                g = w[i+1] - p;
                r = (w[i] - g)*s + Scalar(2.0)*c*b;
                p = s * r;
                w[i+1] = g + p;
                g = c*r - b;
                
                // Form eigenvectors
#ifndef EVALS_ONLY
                for (int k=0; k < n; k++)
                {
                    t = Q[k][i+1];
                    Q[k][i+1] = s*Q[k][i] + c*t;
                    Q[k][i]   = c*Q[k][i] - s*t;
                }
#endif
            }
            w[l] -= p;
            e[l]  = g;
            e[m]  = Scalar(0.0);
        }
    }
    
    return 0;
}


// ----------------------------------------------------------------------------
template<typename Scalar>
IGL_INLINE void dsytrd3(Scalar A[3][3], Scalar Q[3][3], Scalar d[3], Scalar e[2])
// ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
    const int n = 3;
    Scalar u[n], q[n];
    Scalar omega, f;
    Scalar K, h, g;
    
    // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
    for (int i=0; i < n; i++)
    {
        Q[i][i] = Scalar(1.0);
        for (int j=0; j < i; j++)
            Q[i][j] = Q[j][i] = Scalar(0.0);
    }
#endif
    
    // Bring first row and column to the desired form
    h = SQR(A[0][1]) + SQR(A[0][2]);
    if (A[0][1] > 0)
        g = -sqrt(h);
    else
        g = sqrt(h);
    e[0] = g;
    f    = g * A[0][1];
    u[1] = A[0][1] - g;
    u[2] = A[0][2];
    
    omega = h - f;
    if (omega > Scalar(0.0))
    {
        omega = Scalar(1.0) / omega;
        K     = Scalar(0.0);
        for (int i=1; i < n; i++)
        {
            f    = A[1][i] * u[1] + A[i][2] * u[2];
            q[i] = omega * f;                  // p
            K   += u[i] * f;                   // u* A u
        }
        K *= Scalar(0.5) * SQR(omega);
        
        for (int i=1; i < n; i++)
            q[i] = q[i] - K * u[i];
        
        d[0] = A[0][0];
        d[1] = A[1][1] - Scalar(2.0)*q[1]*u[1];
        d[2] = A[2][2] - Scalar(2.0)*q[2]*u[2];
        
        // Store inverse Householder transformation in Q
#ifndef EVALS_ONLY
        for (int j=1; j < n; j++)
        {
            f = omega * u[j];
            for (int i=1; i < n; i++)
                Q[i][j] = Q[i][j] - f*u[i];
        }
#endif
        
        // Calculate updated A[1][2] and store it in e[1]
        e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
    }
    else
    {
        for (int i=0; i < n; i++)
            d[i] = A[i][i];
        e[1] = A[1][2];
    }
}




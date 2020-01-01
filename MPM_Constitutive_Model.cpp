/*
 * @file MPM_CONSTITUTIVE_MODEL.h
 * @brief 
 *
 * CONFIDENTIAL INFORMATION: This software is the confidential and 
 * proprietary information of Walt Disney Animation Studios ("WDAS"). 
 * This software may not be used, disclosed, reproduced or distributed 
 * for any purpose without prior written authorization and license from WDAS. 
 * Reproduction of any section of this software must include this legend
 * and all copyright notices.
 * (C) 2011 Disney Enterprises, Inc. All rights reserved.
 *
 * @author Andrew Selle
 *
 */
#ifndef __MPM_CONSTITUTIVE_MODEL__
#define __MPM_CONSTITUTIVE_MODEL__




#include <PhysBAM_Tools/Math_Tools/Robust_Functions.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>


namespace PhysBAM{

float hardening_max=FLT_MAX;

template<class T>
inline T log1p_x_over_x(const T x) // log(x+1)/x
{if(fabs(x)<1e-18) return 1;return log1p(x)/x;}

template<class T>
inline T diff_log_over_diff(const T x,const T y) // (log(x)-log(y))/(x-y)
{return log1p_x_over_x((y-x)/x)/x;}

template<class T>
MATRIX<T,2> Times_Cofactor_Matrix_Derivative(const MATRIX<T,2>& F,const MATRIX<T,2>& dF)
{
    return MATRIX<T,2>(dF.x[3],-dF.x[2],-dF.x[1],dF.x[0]);
}

template<class T>
MATRIX<T,3> Times_Cofactor_Matrix_Derivative(const MATRIX<T,3>& F,const MATRIX<T,3>& dF)
{
    return MATRIX<T,3>(
        dF.x[4]*F.x[8]+F.x[4]*dF.x[8]-dF.x[5]*F.x[7]-F.x[5]*dF.x[7],
        dF.x[5]*F.x[6]+F.x[5]*dF.x[6]-dF.x[3]*F.x[8]-F.x[3]*dF.x[8],
        dF.x[3]*F.x[7]+F.x[3]*dF.x[7]-dF.x[4]*F.x[6]-F.x[4]*dF.x[6],
        dF.x[2]*F.x[7]+F.x[2]*dF.x[7]-dF.x[1]*F.x[8]-F.x[1]*dF.x[8],
        dF.x[0]*F.x[8]+F.x[0]*dF.x[8]-dF.x[2]*F.x[6]-F.x[2]*dF.x[6],
        dF.x[1]*F.x[6]+F.x[1]*dF.x[6]-dF.x[0]*F.x[7]-F.x[0]*dF.x[7],
        dF.x[1]*F.x[5]+F.x[1]*dF.x[5]-dF.x[2]*F.x[4]-F.x[2]*dF.x[4],
        dF.x[2]*F.x[3]+F.x[2]*dF.x[3]-dF.x[0]*F.x[5]-F.x[0]*dF.x[5],
        dF.x[0]*F.x[4]+F.x[0]*dF.x[4]-dF.x[1]*F.x[3]-F.x[1]*dF.x[3]);
}

template<class T>
MATRIX<T,2> dR(const MATRIX<T,2>& dF,const MATRIX<T,2>& R,const MATRIX<T,2>& F)
{
    MATRIX<T,2> M=R.Transpose_Times(dF);
    MATRIX<T,2> S=R.Transpose_Times(F);
    T a=(M.x[2]-M.x[1])/(S.x[3]+S.x[0]);
    MATRIX<T,2> RTdR(0,-a,a,0);
    return R*RTdR;
}

template<class T>
MATRIX<T,3> dR(const MATRIX<T,3>& dF,const MATRIX<T,3>& R,const MATRIX<T,3>& F)
{
    MATRIX<T,3> M=R.Transpose_Times(dF);
    MATRIX<T,3> S=R.Transpose_Times(F);

    MATRIX<T,3> K(S.x[4]+S.x[0],S.x[5],-S.x[2],S.x[5],S.x[8]+S.x[0],S.x[1],-S.x[2],S.x[1],S.x[8]+S.x[4]);
    VECTOR<T,3> L(M.x[1]-M.x[3],-M.x[6]+M.x[2],-M.x[7]+M.x[5]);
    VECTOR<T,3> RV=K.Solve_Linear_System(L);
    MATRIX<T,3> RTdR(0,RV.x,RV.y,-RV.x,0,RV.z,-RV.y,-RV.z,0);
    return R*RTdR;
}


}
#endif

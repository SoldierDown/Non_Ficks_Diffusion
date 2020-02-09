//!#####################################################################
//! \file MPM_Constitutive_Model.h
//!#####################################################################
// Class MPM_Constitutive_Model
//######################################################################
#ifndef __MPM_Constitutive_Model__
#define __MPM_Constitutive_Model__
#include <float.h>
#include <nova/Tools/Matrices/Matrix_2x2.h>
#include <nova/Tools/Matrices/Matrix_3x3.h>
#include <nova/Tools/Matrices/Diagonal_Matrix.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_2x2.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_3x3.h>
#include <nova/Tools/Log/Log.h>
#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/Tools/Vectors/Vector.h>
namespace Nova
{
template<class T>
inline T sqr(const T a)
{return a*a;}
template<class T>
inline T log1p_x_over_x(const T x) // log(x+1)/x
{if(fabs(x)<1e-18) return 1;return log1p(x)/x;}
template<class T>
inline T diff_log_over_diff(const T x,const T y) // (log(x)-log(y))/(x-y)
{return log1p_x_over_x((y-x)/x)/x;}
template<class T>
Matrix<T,2> Times_Cofactor_Matrix_Derivative(const Matrix<T,2>& F,const Matrix<T,2>& dF)
{
    return Matrix<T,2>(dF.x[3],-dF.x[2],-dF.x[1],dF.x[0]);
}
template<class T>
Matrix<T,3> Times_Cofactor_Matrix_Derivative(const Matrix<T,3>& F,const Matrix<T,3>& dF)
{
    return Matrix<T,3>(
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
Matrix<T,2> dR(const Matrix<T,2>& dF,const Matrix<T,2>& R,const Matrix<T,2>& F)
{
    Matrix<T,2> M=R.Transpose_Times(dF);
    Matrix<T,2> S=R.Transpose_Times(F);
    T a=(M.x[2]-M.x[1])/(S.x[3]+S.x[0]);
    Matrix<T,2> RTdR(0,-a,a,0);
    return R*RTdR;
}
template<class T>
Matrix<T,3> dR(const Matrix<T,3>& dF,const Matrix<T,3>& R,const Matrix<T,3>& F)
{
    Matrix<T,3> M=R.Transpose_Times(dF);
    Matrix<T,3> S=R.Transpose_Times(F);
    Matrix<T,3> K(S.x[4]+S.x[0],S.x[5],-S.x[2],S.x[5],S.x[8]+S.x[0],S.x[1],-S.x[2],S.x[1],S.x[8]+S.x[4]);
    Vector<T,3> L;
    L(0)=M.x[1]-M.x[3];
    L(1)=-M.x[6]+M.x[2];
    L(2)=-M.x[7]+M.x[5];
    Vector<T,3> RV=K.Solve_Linear_System(L);
    Matrix<T,3> RTdR(0,RV(0),RV(1),-RV(0),0,RV(2),-RV(1),-RV(2),0);
    return R*RTdR;
}
template<class T,int d>
class MPM_Constitutive_Model
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Mat                       = Matrix<T,d>;

public:
    float hardening_max=FLT_MAX;
    T lambda,mu,lambda0,mu0;
    Mat Fe;       
    Mat Ue,Ve;
    Diagonal_Matrix<T,d> Se;
    T Je;
    Mat Fp;       

    // only used by fixed corotated
    Mat Re,He;

    T eta;

    // snow plasticity
    bool plastic;
    T stretching_yield;
    T compression_yield;
    T hardening_factor;

    MPM_Constitutive_Model()
        :lambda(0),mu(0),lambda0(0),mu0(0),plastic(false),stretching_yield(FLT_MAX),compression_yield(-FLT_MAX),hardening_factor(0),eta(0)
    {
        Fe=Matrix<T,d>::Identity_Matrix();
        Fp=Matrix<T,d>::Identity_Matrix();
        Je=(T)1.;
        plastic=false;
    }
    ~MPM_Constitutive_Model(){}
//######################################################################
    void Compute_Lame_Parameters(const T E,const T nu);
    void Precompute();
    void Test();
    T    Psi() const;
    Mat  P() const;
    Mat  Times_dP_dF(const Mat& dF) const;
//######################################################################
}; 
} // namespace Nova



#endif
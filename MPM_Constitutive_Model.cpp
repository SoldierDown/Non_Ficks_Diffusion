
//!#####################################################################
//! \file MPM_Constitutive_Model.cpp
//!#####################################################################
#include "MPM_Constitutive_Model.h"
using namespace Nova;
//######################################################################
// Constructor
//######################################################################
// template<class T,int d> MPM_Constitutive_Model<T,d>::
// MPM_Constitutive_Model()

//######################################################################
// Destructor
//######################################################################
// template<class T,int d> MPM_Constitutive_Model<T,d>::
// ~MPM_Constitutive_Model()
// {
// }
//######################################################################
// Compute_Lame_Parameters
//######################################################################
template<class T,int d> void MPM_Constitutive_Model<T,d>::
Compute_Lame_Parameters(const T E,const T nu)
{
    mu0=E/((T)2.*((T)1.+nu));
    lambda0=(E*nu)/(((T)1.+nu)*((T)1.-(T)2.*nu));
}
//######################################################################
// Precompute
//######################################################################
template<class T,int d> void MPM_Constitutive_Model<T,d>::
Precompute()
{
    Ue=Matrix<T,d>();Ve=Matrix<T,d>();
    Fe.Fast_Singular_Value_Decomposition(Ue,Se,Ve);
    if(plastic){
        for(int i=0;i<d;++i) if(Se(i)>stretching_yield) Se(i)=stretching_yield;
        for(int i=0;i<d;++i) if(Se(i)<compression_yield) Se(i)=compression_yield;
        Fp=Ve*Se.Inverse()*Ue.Transposed()*Fe*Fp;
        Fe=Ue*Se*Ve.Transposed();
        T power=std::min(hardening_factor*(1-Fp.Determinant()),(T)hardening_max);
        lambda=exp(power)*lambda0;
        mu=exp(power)*mu0;}
    else{
        lambda=lambda0;
        mu=mu0;
        // printf("lambda: %f, mu: %f\n",lambda, mu);
        }
    Je=Se.Determinant();
    Re=Ue*Ve.Transposed();
    He=Fe.Cofactor_Matrix();
}
//######################################################################
// Psi
//######################################################################
template<class T,int d> T MPM_Constitutive_Model<T,d>::
Psi() const
{
    return mu*(Se.To_Vector()-1).Norm_Squared()+lambda/2*sqr(Je-1);
}
//######################################################################
// Precompute
//######################################################################
template<class T,int d> Matrix<T,d> MPM_Constitutive_Model<T,d>::
P() const
{
    return 2*mu*(Fe-Re)+lambda*(Je-1)*He;
}
//######################################################################
// Precompute
//######################################################################
template<class T,int d> Matrix<T,d> MPM_Constitutive_Model<T,d>::
Times_dP_dF(const Matrix<T,d>& dF) const
{
    Matrix<T,d> dHdF_times_dF=Times_Cofactor_Matrix_Derivative(Fe,dF);
    return 2*mu*(dF-dR(dF,Re,Fe))+lambda*(He.Times_Transpose(dF)).Trace()*He+lambda*(Je-1)*dHdF_times_dF;
}
//######################################################################
// Precompute
//######################################################################
template<class T,int d> void MPM_Constitutive_Model<T,d>::
Test()
{

}
//######################################################################
template class Nova::MPM_Constitutive_Model<float,2>;
template class Nova::MPM_Constitutive_Model<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::MPM_Constitutive_Model<double,2>;
template class Nova::MPM_Constitutive_Model<double,3>;
#endif



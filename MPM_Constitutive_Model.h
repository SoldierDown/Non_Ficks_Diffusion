//!#####################################################################
//! \file MPM_Constitutive_Model.h
//!#####################################################################
// Class MPM_Constitutive_Model
//######################################################################
#ifndef __MPM_Constitutive_Model__
#define __MPM_Constitutive_Model__
namespace Nova
{
template<class T,int d>
class MPM_CONSTITUTIVE_MODEL
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;

public:

    ///////////////////////////////////////////////////////////////
    int model; 
    enum MODEL_TYPE{
        FIXED_COROTATED=1,
        HENCKY=2};
    ///////////////////////////////////////////////////////////////

    T lambda,mu,lambda0,mu0;
    MATRIX<T,TV::m> Fe;       
    MATRIX<T,TV::m> Ue,Ve;
    DIAGONAL_MATRIX<T,TV::m> Se;
    T Je;
    MATRIX<T,TV::m> Fp;       

    // only used by fixed corotated
    MATRIX<T,TV::m> Re,He;

    // only used by hencky
    T failure_threshold;
    T hencky_psi;
    MATRIX<T,TV::m> hencky_P;
    DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,TV::m> hencky_dP_dF_diagonal;

    // snow plasticity
    bool plastic;
    T stretching_yield;
    T compression_yield;
    T hardening_factor;


    MPM_CONSTITUTIVE_MODEL(const T E=0,const T nu=0)
        :lambda(0),mu(0),lambda0(0),mu0(0),plastic(false),stretching_yield(FLT_MAX),
        compression_yield(-FLT_MAX),hardening_factor(0)
    {
        model=FIXED_COROTATED;
        Compute_Lame_Parameters(E,nu);
        Fe+=1;
        Fp+=1;
        Je=1;
        plastic=false;
    }
    ~MPM_CONSTITUTIVE_MODEL(){}
    void Compute_Lame_Parameters(const T E,const T nu)
    {
        mu0=E/(2*(1+nu));
        lambda0=(E*nu)/((1+nu)*(1-2*nu));

        if(model==HENCKY) failure_threshold=0;
    }
    void Precompute()
    {
        Ue=MATRIX<T,TV::m>();Ve=MATRIX<T,TV::m>();
        Fe.Fast_Singular_Value_Decomposition(Ue,Se,Ve);
        if(plastic){
            for(int i=1;i<=TV::m;i++) if(Se(i)>stretching_yield) Se(i)=stretching_yield;
            for(int i=1;i<=TV::m;i++) if(Se(i)<compression_yield) Se(i)=compression_yield;
            Fp=Ve*Se.Inverse()*Ue.Transposed()*Fe*Fp;
            Fe=Ue*Se*Ve.Transposed();
            T power=min(hardening_factor*(1-Fp.Determinant()),(T)hardening_max);
            lambda=exp(power)*lambda0;
            mu=exp(power)*mu0;}
        else{
            lambda=lambda0;
            mu=mu0;}
        Je=Se.Determinant();

        if(model==FIXED_COROTATED){
            Re=Ue*Ve.Transposed();
            He=Fe.Cofactor_Matrix();}
        else if(model==HENCKY){
            if(Se(TV::m)<=0){
                hencky_psi=FLT_MAX;
                hencky_P=MATRIX<T,TV::m>();}
            else{
                DIAGONAL_MATRIX<T,TV::m> log_Se=log(Se);
                T sum_sqr=sqr(log_Se.Trace());
                hencky_psi=mu*log_Se.Frobenius_Norm_Squared()+(T)0.5*lambda*sum_sqr;
                DIAGONAL_MATRIX<T,TV::m> F_Inv=Se.Inverse();
                DIAGONAL_MATRIX<T,TV::m> P_hat=((T)2*mu*log_Se+lambda*log_Se.Trace())*F_Inv;
                hencky_P=(Ue*P_hat).Times_Transpose(Ve);}
            Evaluate_Hencky_Isotropic_Stress_Derivative(Se,hencky_dP_dF_diagonal);}
    }

    T Psi() const
    {
        if(model==FIXED_COROTATED)
            return mu*(Se.To_Vector()-1).Magnitude_Squared()+lambda/2*sqr(Je-1);
        else if(model==HENCKY)
            return hencky_psi;
    }

    MATRIX<T,TV::m> P() const
    {
        if(model==FIXED_COROTATED)
            return 2*mu*(Fe-Re)+lambda*(Je-1)*He;
        else if(model==HENCKY)
            return hencky_P;
    }

    MATRIX<T,TV::m> Times_dP_dF(const MATRIX<T,TV::m>& dF) const
    {
        if(model==FIXED_COROTATED){
            MATRIX<T,TV::m> dHdF_times_dF=Times_Cofactor_Matrix_Derivative(Fe,dF);
            return 2*mu*(dF-dR(dF,Re,Fe))+lambda*(He.Times_Transpose(dF)).Trace()*He+lambda*(Je-1)*dHdF_times_dF;}
        else if(model==HENCKY){
            MATRIX<T,TV::m> UtdFV=Ue.Transposed()*dF*Ve;
            MATRIX<T,TV::m> action=hencky_dP_dF_diagonal.Differential(UtdFV);
            return (Ue*action).Times_Transpose(Ve);}
    }

    void Evaluate_Hencky_Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dPi_dF)
    {
        if(F(2)<=0)
            dPi_dF.x1111=dPi_dF.x2222=dPi_dF.x2211=dPi_dF.x2121=dPi_dF.x2112=(T)0;
        else{
            DIAGONAL_MATRIX<T,2> F_threshold=F.Clamp_Min(failure_threshold);
            // DIAGONAL_MATRIX<T,2> F_threshold=F;
            T lambda_plus_two_mu=lambda+2*mu,two_mu=2*mu,s0=F_threshold(1),s1=F_threshold(2);
            T l0=log(s0),l1=log(s1),sum_log=l0+l1;
            SYMMETRIC_MATRIX<T,2> F_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(VECTOR<T,2>(F_threshold(1),F_threshold(2)));
            T Psi_0_minus_Psi_1=(two_mu/s0)*(diff_log_over_diff(s0,s1)-l1/s1)-lambda*sum_log/F_outer.x21;
            T Psi_0_plus_Psi_1=((lambda*s0+lambda_plus_two_mu*s1)*l0+(lambda*s1+lambda_plus_two_mu*s0)*l1)/(F_outer.x21*(s0+s1));
            dPi_dF.x1111=(lambda_plus_two_mu*(1-l0)-lambda*l1)/F_outer.x11;
            dPi_dF.x2222=(lambda_plus_two_mu*(1-l1)-lambda*l0)/F_outer.x22;
            dPi_dF.x2211=lambda/F_outer.x21;
            dPi_dF.x2121=(T)0.5*(Psi_0_minus_Psi_1+Psi_0_plus_Psi_1);
            dPi_dF.x2112=(T)0.5*(Psi_0_minus_Psi_1-Psi_0_plus_Psi_1);}
    }

    void Evaluate_Hencky_Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF)
    {
        if(F(3)<=0){
            dPi_dF.x1111=dPi_dF.x2222=dPi_dF.x3333=dPi_dF.x2211=dPi_dF.x3311=dPi_dF.x3322=dPi_dF.x2121=dPi_dF.x3131=dPi_dF.x3232=dPi_dF.x2112=dPi_dF.x3113=dPi_dF.x3223=(T)0;}
        else{
            DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold);
            // DIAGONAL_MATRIX<T,3> F_threshold=F;
            T lambda_plus_two_mu=lambda+2*mu,two_mu=2*mu,s0=F_threshold(1),s1=F_threshold(2),s2=F_threshold(3);
            T l0=log(s0),l1=log(s1),l2=log(s2),sum_log=l0+l1+l2;
            SYMMETRIC_MATRIX<T,3> F_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR<T,3>(s0,s1,s2));
            MATRIX<T,3> sl_outer=MATRIX<T,3>::Outer_Product(VECTOR<T,3>(s0,s1,s2),VECTOR<T,3>(l0,l1,l2));
            //Hessian variables
            dPi_dF.x1111=(lambda_plus_two_mu-two_mu*l0-lambda*sum_log)/F_outer.x11;
            dPi_dF.x2222=(lambda_plus_two_mu-two_mu*l1-lambda*sum_log)/F_outer.x22;
            dPi_dF.x3333=(lambda_plus_two_mu-two_mu*l2-lambda*sum_log)/F_outer.x33;
            dPi_dF.x2211=lambda/F_outer.x21;
            dPi_dF.x3311=lambda/F_outer.x31;
            dPi_dF.x3322=lambda/F_outer.x32;
            //2x2 block-matrices
            T Psi_0_minus_Psi_1=(two_mu/s0)*(diff_log_over_diff(s0,s1)-l1/s1)-lambda*sum_log/F_outer.x21;
            T Psi_0_minus_Psi_2=(two_mu/s0)*(diff_log_over_diff(s0,s2)-l2/s2)-lambda*sum_log/F_outer.x31;
            T Psi_1_minus_Psi_2=(two_mu/s1)*(diff_log_over_diff(s1,s2)-l2/s2)-lambda*sum_log/F_outer.x32;
            T row_sum1=sl_outer(1,1)+sl_outer(1,2)+sl_outer(1,3);
            T row_sum2=sl_outer(2,1)+sl_outer(2,2)+sl_outer(2,3);
            T row_sum3=sl_outer(3,1)+sl_outer(3,2)+sl_outer(3,3);
            T Psi_0_plus_Psi_1=(two_mu*(s1*l0+s0*l1)+lambda*(row_sum1+row_sum2))/(F_outer.x21*(s0+s1));
            T Psi_0_plus_Psi_2=(two_mu*(s2*l0+s0*l2)+lambda*(row_sum1+row_sum3))/(F_outer.x31*(s0+s2));
            T Psi_1_plus_Psi_2=(two_mu*(s2*l1+s1*l2)+lambda*(row_sum2+row_sum3))/(F_outer.x32*(s1+s2));
            dPi_dF.x2121=(T)0.5*(Psi_0_minus_Psi_1+Psi_0_plus_Psi_1);
            dPi_dF.x2112=(T)0.5*(Psi_0_minus_Psi_1-Psi_0_plus_Psi_1); 
            dPi_dF.x3131=(T)0.5*(Psi_0_minus_Psi_2+Psi_0_plus_Psi_2);
            dPi_dF.x3113=(T)0.5*(Psi_0_minus_Psi_2-Psi_0_plus_Psi_2); 
            dPi_dF.x3232=(T)0.5*(Psi_1_minus_Psi_2+Psi_1_plus_Psi_2);
            dPi_dF.x3223=(T)0.5*(Psi_1_minus_Psi_2-Psi_1_plus_Psi_2);}
    }

    void Read(std::istream& input)
    {Read_Binary<T>(input,Fe,Fp,lambda0,mu0,plastic,stretching_yield,compression_yield,hardening_factor);}

    void Write(std::ostream& output) const
    {Write_Binary<T>(output,Fe,Fp,lambda0,mu0,plastic,stretching_yield,compression_yield,hardening_factor);}

    void Test()
    {
        RANDOM_NUMBERS<T> random;
        T e=(T)1e-6;
        for(int t=1;t<=100;t++){
            MATRIX<T,TV::m> F0,dF,oldF=Fe;
            random.Fill_Uniform(F0,-1,1);
            // if(F0.Determinant()<0.1) continue;
            random.Fill_Uniform(dF,-e,e);
            Fe=F0;
            Precompute();
            T psi1=Psi();
            MATRIX<T,TV::m> P1=P();
            MATRIX<T,TV::m> dP1=Times_dP_dF(dF);
            Fe=F0+dF;
            Precompute();
            T psi2=Psi();
            MATRIX<T,TV::m> P2=P();
            MATRIX<T,TV::m> dP2=Times_dP_dF(dF);
            Fe=oldF;
            Precompute();
            LOG::cout<<"E test: "<<((psi2-psi1)-(P2+P1).Times_Transpose(dF).Trace()/2)/e<<"    "<<maxabs(psi2,psi1)<<std::endl;
            LOG::cout<<"P test: "<<((P2-P1)-(dP2+dP1)/2).Frobenius_Norm()/e<<"    "<<max(P2.Frobenius_Norm(),P1.Frobenius_Norm())<<std::endl;}
    }
}; 
} // namespace Nova



#endif
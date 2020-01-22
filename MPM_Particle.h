//!#####################################################################
//! \file MPM_Particle.h
//!#####################################################################
// Class MPM_Particle
//######################################################################
#ifndef __MPM_Particle__
#define __MPM_Particle__

#include <nova/Tools/Vectors/Vector.h>
#include "MPM_Constitutive_Model.h"
#include <nova/Tools/Matrices/Matrix_3x3.h>

namespace Nova{
template<class T,int d>
class MPM_Particle
{
    using TV                        = Vector<T,d>;
    using Mat                       = Matrix<T,3>;
    using T_INDEX                   = Vector<int,d>;
public:
    bool valid;
    TV X,V;
    T_INDEX closest_node;
    T mass,volume;
    Mat weights;
    Mat dweights;
    T_INDEX base_node;              
    
    // Constitutive model
    MPM_Constitutive_Model<T,d> constitutive_model;
    Matrix<T,d> scp;

    MPM_Particle()
    {Initialize();}

    ~MPM_Particle() {}

    void Initialize()
    {
        valid=true;
        X=TV();
        V=TV();
        mass=(T)0.;
        volume=(T)0.;
        scp=Matrix<T,d>();

    }

    T Weight(T_INDEX index)
    {
        T weight=(T)1.;
        index+=1;
        for(int i=0;i<d;++i)
            weight*=weights(index(i),i);
        return weight;
    }

    TV Weight_Gradient(T_INDEX index)
    {
        TV weight_gradient=TV(1);
        index+=1;
        for(int i=0;i<d;++i) for(int j=0;j<d;++j)
        weight_gradient(i)*=(i==j)?dweights(index(j),j):weights(index(j),j);
        return weight_gradient;
    }

    void Update_Weights(const Grid<T,d>& grid)    
    {
        closest_node=grid.Closest_Node(X);
        TV dx=grid.dX;
        const TV X_eval=X-(grid.Node(closest_node)-dx);
        for(int i=0;i<d;++i){
            const T one_over_dx=grid.one_over_dX(i);
            const T x=X_eval(i);
            Compute_Weights(x*one_over_dx,one_over_dx,i);}
    }

    void Compute_Weights(T x,const T one_over_dx,const int k)
    {
        weights(0,k)=(T).5*x*x-(T)1.5*x+(T)1.125;
        dweights(0,k)=(x-(T)1.5)*one_over_dx;
        x-=(T)1;
        weights(1,k)=-x*x+(T).75;
        dweights(1,k)=(T)(-2)*x*one_over_dx;
        x-=(T)1;
        weights(2,k)=(T).5*x*x+(T)1.5*x+(T)1.125;
        dweights(2,k)=(x+(T)1.5)*one_over_dx;
    }
};
}
#include "Read_Write/Read_Write_MPM_Particles.h"
#endif

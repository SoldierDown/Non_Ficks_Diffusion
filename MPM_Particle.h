//!#####################################################################
//! \file MPM_Particle.h
//!#####################################################################
// Class MPM_Particle
//######################################################################
#ifndef __MPM_Particle__
#define __MPM_Particle__

#include <nova/Tools/Vectors/Vector.h>
#include "MPM_Constitutive_Model.h"

namespace Nova{
template<class T,int d>
class MPM_Particle
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
public:
    bool valid;
    TV X,V;
    T mass,volume;
    
    T_INDEX base_node;              

    // Hydrogel variables
    T mass_solid;
    T mass_fluid;
    T saturation;
    T volume_fraction_0;
    T div_Qc;
    
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

        // Hydrogel variables
        mass_solid=mass;
        mass_fluid=(T)0.;
        saturation=(T)0.;
        volume_fraction_0=(T)0.;
        div_Qc=(T)0.;
        scp=Matrix<T,d>();
    }
};
}
#include "Read_Write/Read_Write_MPM_Particles.h"
#endif

//!#####################################################################
//! \file MPM_Particle.h
//!#####################################################################
// Class MPM_Particle
//######################################################################
#ifndef __MPM_Particle__
#define __MPM_Particle__

#include <nova/Tools/Vectors/Vector.h>

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
    Vector<TV,d> weights;
    Vector<TV,d> dweight;

    // Constitutive model
    


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
    }
};
}
#include "Read_Write/Read_Write_MPM_Particles.h"
#endif

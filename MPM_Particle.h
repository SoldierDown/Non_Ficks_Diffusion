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
    TV X,V;
    T mass,volume;

    MPM_Particle()
    {Initialize();}

    ~MPM_Particle() {}

    void Initialize()
    {
        mass=(T)0.;
        volume=(T)0.;
    }
};
}
#include "Read_Write/Read_Write_MPM_Particles.h"
#endif

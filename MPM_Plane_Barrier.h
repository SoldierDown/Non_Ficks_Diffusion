//!#####################################################################
//! \file MPM_Plane_Barrier.h
//!#####################################################################
// Class MPM_Plane_Barrier
//###################################################################### 
#ifndef __MPM_Plane_Barrier__
#define __MPM_Plane_Barrier__

#include <nova/Tools/Vectors/Vector.h>

namespace Nova{

template<class T,int d>
class MPM_Plane_Barrier
{
    using TV                    = Vector<T,d>;
public:

    const T mu; // coefficient of friction
    Vector<T,d> axis_vector;
    Vector<T,d> normal;
    Vector<T,d> surface;
    const bool sticky;
    
    MPM_Plane_Barrier(const T mu_input,const TV normal_input,const TV surface_input,const bool sticky_input=false)
        :mu(mu_input),normal(normal_input),surface(surface_input),sticky(sticky_input)
    {}

    ~MPM_Plane_Barrier()
    {
    }
};

}
#endif
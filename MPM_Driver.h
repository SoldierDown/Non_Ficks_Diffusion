//!#####################################################################
//! \file MPM_Driver.h
//!#####################################################################
// Class MPM_Driver
//######################################################################
#ifndef __MPM_Driver__
#define __MPM_Driver__

#include <nova/Tools/Utilities/Driver.h>
#include "MPM_Example.h"

namespace Nova{
template<class T,int d>
class MPM_Driver: public Driver<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = Driver<T,d>;

    using Base::time;
    using Base::Compute_Dt;using Base::Write_Output_Files;

  public:
    MPM_Example<T,d>& example;

    MPM_Driver(MPM_Example<T,d>& example_input);
    ~MPM_Driver() {}

//######################################################################
    void Initialize();
    void Advance_To_Target_Time(const T target_time) override;
    void Simulate_To_Frame(const int frame) override;
    void Advance_Step(const T dt);
//######################################################################
};
}
#endif

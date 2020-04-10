//!#####################################################################
//! \file PF_Driver.h
//!#####################################################################
// Class PF_Driver
//######################################################################
#ifndef __PF_Driver__
#define __PF_Driver__

#include <nova/Tools/Utilities/Driver.h>
#include "PF_Example.h"

namespace Nova{
template<class T,int d>
class PF_Driver: public Driver<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = Driver<T,d>;

    using Base::time;
    using Base::Compute_Dt;using Base::Write_Output_Files;

  public:
    PF_Example<T,d>& example;

    PF_Driver(PF_Example<T,d>& example_input);
    ~PF_Driver() {}

//######################################################################
    void Initialize() override;
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time);
    void Advance_To_Target_Time(const T target_time) override;
    void Simulate_To_Frame(const int frame) override;
//######################################################################
};
}
#endif

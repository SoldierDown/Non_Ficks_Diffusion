//!#####################################################################
//! \file PC_Driver.h
//!#####################################################################
// Class PC_Driver
//######################################################################
#ifndef __PC_Driver__
#define __PC_Driver__

#include <nova/Tools/Utilities/Driver.h>
#include "PC_Example.h"

namespace Nova{
template<class T,int d>
class PC_Driver: public Driver<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = Driver<T,d>;

    using Base::time;
    using Base::Compute_Dt;using Base::Write_Output_Files;

  public:
    PC_Example<T,d>& example;

    PC_Driver(PC_Example<T,d>& example_input);
    ~PC_Driver() {}

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

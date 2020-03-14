//!#####################################################################
//! \file Smoke_Driver.h
//!#####################################################################
// Class Smoke_Driver
//######################################################################
#ifndef __Smoke_Driver__
#define __Smoke_Driver__

#include <nova/Tools/Utilities/Driver.h>
#include "Smoke_Example.h"

namespace Nova{
template<class T,int d>
class Smoke_Driver: public Driver<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = Driver<T,d>;

    using Base::time;
    using Base::Compute_Dt;using Base::Write_Output_Files;

  public:
    Smoke_Example<T,d>& example;

    Smoke_Driver(Smoke_Example<T,d>& example_input);
    ~Smoke_Driver() {}

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

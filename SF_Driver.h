//!#####################################################################
//! \file SF_Driver.h
//!#####################################################################
// Class SF_Driver
//######################################################################
#ifndef __SF_Driver__
#define __SF_Driver__

#include <nova/Tools/Utilities/Driver.h>
#include "SF_Example.h"

namespace Nova{
template<class T,int d>
class SF_Driver: public Driver<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = Driver<T,d>;

    using Base::time;
    using Base::Compute_Dt;using Base::Write_Output_Files;

  public:
    SF_Example<T,d>& example;

    SF_Driver(SF_Example<T,d>& example_input);
    ~SF_Driver() {}

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

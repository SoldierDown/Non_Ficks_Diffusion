//!#####################################################################
//! \file Multigrid_Data.h
//!#####################################################################
// Class Multigrid_Data
//######################################################################
#ifndef __Multigrid_Data__
#define __Multigrid_Data__

#include <stdint.h>

namespace Nova{
template<class T,class T_FLAGS=uint32_t>
struct Multigrid_Data
{
    typedef T_FLAGS Flags_type;

    T_FLAGS flags;  // flags
    T ch0;          // saturation_channel
    T ch1;          // rhs_channel
    T ch2;          // result_channel
    T ch3;          // x_gradient_channel
    T ch4;          // y_gradient_channel
    T ch5;          // z_gradient_channel
};
}
#endif

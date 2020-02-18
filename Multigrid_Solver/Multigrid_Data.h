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
    T ch0;          // x_channel
    T ch1;          // b_channel
    T ch2;          // r_channel
};
}
#endif

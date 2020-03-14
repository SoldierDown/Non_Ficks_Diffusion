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

    T_FLAGS flags;
    T ch0;
    T ch1;
    T ch2;
    T ch3;
    T ch4;
    T ch5;
};
}
#endif

//!#####################################################################
//! \file Viewer_Data.h
//!#####################################################################
// Class Viewer_Data
//######################################################################
#ifndef __Viewer_Data__
#define __Viewer_Data__

#include <stdint.h>

namespace Nova{
template<class T,class T_FLAGS=uint32_t>
struct Viewer_Data
{
    typedef T_FLAGS Flags_type;

    T_FLAGS flags;
    T ch0;          // density
    T ch1;          // velocity
    T ch2;
    T ch3;
    T ch4;
    T ch5;
    T ch6;
    T ch7;
};
}
#endif

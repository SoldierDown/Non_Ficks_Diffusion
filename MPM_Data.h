//!#####################################################################
//! \file MPM_Data.h
//!#####################################################################
// Class MPM_Data
//######################################################################
#ifndef __MPM_Data__
#define __MPM_Data__

#include <stdint.h>

namespace Nova{
template<class T,class T_FLAGS=uint32_t>
struct MPM_Data
{
    typedef T_FLAGS Flags_type;

    T_FLAGS         flags;
    T ch0;          // mass
    T ch1;          // X-velocity
    T ch2;          // Y-velocity
    T ch3;          // Z-velocity
    T ch4;          // X-velocity-star
    T ch5;          // Y-velocity-star
    T ch6;          // Z-velocity-star
    T ch7;          // X-force
    T ch8;          // Y-force
    T ch9;          // Z-force
    T ch10;         // valid nodes
    T ch11;         // collide nodes
    T ch12;         
    T ch13;
    T ch14;
};
}
#endif

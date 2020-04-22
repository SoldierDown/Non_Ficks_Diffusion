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

    T ch10;         // rhs-x
    T ch11;         // rhs-y
    T ch12;         // rhs-z
    T ch13;         // q-x
    T ch14;         // q-y
    T ch15;         // q-z
    T ch16;         // s-x
    T ch17;         // s-y
    T ch18;         // s-z
    T ch19;         // r-x
    T ch20;         // r-y
    T ch21;         // r-z
    T ch22;         // z-x
    T ch23;         // z-y
    T ch24;         // z-z
    T ch25;         // tmp;
};
}
#endif
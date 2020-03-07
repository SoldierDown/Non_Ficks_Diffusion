//!#####################################################################
//! \file Poisson_Data.h
//!#####################################################################
// Class Poisson_Data
//######################################################################
#ifndef __Poisson_Data__
#define __Poisson_Data__

#include <stdint.h>

namespace Nova{
template<class T,class T_FLAGS=uint32_t>
struct Poisson_Data
{
    typedef T_FLAGS Flags_type;

    T_FLAGS flags;
    T ch0;
    T ch1;
    T ch2;
    T ch3;
    T ch4;
    T ch5;
    T ch6;
    T ch7;
    T ch8;
    T ch9;
    T ch10;
    T ch11;
    T ch12;
    T ch13;
    T ch14;
};
}
#endif

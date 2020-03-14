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
    T ch0;      // X-face_velocity
    T ch1;      // Y-face_velocity
    T ch2;      // Z-face_velocity
    T ch3;      // pressure
    T ch4;      // density
    T ch5;      // X-node_velocity
    T ch6;      // Y-node_velocity
    T ch7;      // Z-node_velocity
    T ch8;      // node_density
    T ch9;      // temp
    T ch10;     
    T ch11;
    T ch12;
    T ch13;     
    T ch14;     // lap_density
};
}
#endif

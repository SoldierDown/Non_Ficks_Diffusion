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
    T ch10;         // ficks_rhs
    
    // Hydrogel channels
    T ch11;         // saturation
    T ch12;         // lap_saturation
    T ch13;         // void_mass_fluid
    T ch14;         // volume
    T ch15;         // div_Qc
    
    T ch16;         // q
    T ch17;         // s
    T ch18;         // r
    T ch19;         // k
    T ch20;         // z

    // Matrix components
    // T ch16;         // (0,0)
    // T ch17;         // (0,1)
    // T ch18;         // (0,2)
    // T ch19;         // (1,0)
    // T ch20;         // (1,1)
    // T ch21;         // (1,2)
    // T ch22;         // (2,0)
    // T ch23;         // (2,1)
    // T ch24;         // (2,2)
};
}
#endif

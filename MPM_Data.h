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


    // Hydrogel channels
    T ch10;         // saturation
    T ch11;         // lap_saturation
    T ch12;         // void_mass_fluid
    T ch13;         // volume
    T ch14;         // div_Qc

    T ch15;         // rhs-x
    T ch16;         // rhs-y
    T ch17;         // rhs-z

    T ch18;         // q-x
    T ch19;         // q-y
    T ch20;         // q-z
    T ch21;         // s-x
    T ch22;         // s-y
    T ch23;         // s-z
    T ch24;         // r-x
    T ch25;         // r-y
    T ch26;         // r-z
    T ch27;         // k-x
    T ch28;         // k-y
    T ch29;         // k-z
    T ch30;         // z-x
    T ch31;         // z-y
    T ch32;         // z-z

    // implicit update for velocity: 6*2 or 6*3 
    // matrix 4 or 9
    // 16 or 27 in total

    // matrix components
    T ch33;         // (0,0)
    T ch34;         // (0,1)
    T ch35;         // (0,2)
    T ch36;         // (1,0)
    T ch37;         // (1,1)
    T ch38;         // (1,2)
    T ch39;         // (2,0)
    T ch40;         // (2,1)
    T ch41;         // (2,2)

};
}
#endif
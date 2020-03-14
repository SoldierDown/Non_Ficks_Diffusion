//!#####################################################################
//! \file Diffusion_Data.h
//!#####################################################################
// Class Diffusion_Data
//######################################################################
#ifndef __Diffusion_Data__
#define __Diffusion_Data__

#include <stdint.h>

namespace Nova{
template<class T,class T_FLAGS=uint32_t>
struct Diffusion_Data
{
    typedef T_FLAGS Flags_type;

    T_FLAGS         flags;
    // Hydrogel channels
    T ch0;         // saturation
    T ch1;         // lap_saturation
    T ch2;         // void_mass_fluid
    T ch3;         // volume
    T ch4;         // div_Qc

    T ch5;         // rhs
    T ch6;         // q
    T ch7;         // s
    T ch8;         // r
    T ch9;         // z
};
}
#endif

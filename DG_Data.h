//!#####################################################################
//! \file DG_Data.h
//!#####################################################################
// Class DG_Data
//######################################################################
#ifndef __DG_Data__
#define __DG_Data__

#include <stdint.h>

namespace Nova{
template<class T,class T_FLAGS=uint32_t>
struct DG_Data
{
    typedef T_FLAGS Flags_type;

    T_FLAGS flags;
    T ch0;      // density
    T ch1;      // density_backup
    T ch2;      // T
    T ch3;      // T_backup
    T ch4;      // X_Qpc
    T ch5;      // Y_Qpc
    T ch6;      // Z_Qpc
    T ch7;      // X_Qpc backup
    T ch8;      // Y_Qpc backup
    T ch9;      // Z_Qpc backup
    T ch10;     // X_Qtc
    T ch11;     // Y_Qtc
    T ch12;     // Z_Qtc
    T ch13;     // X_Qtc_backup
    T ch14;     // Y_Qtc_backup
    T ch15;     // Z_Qtc_backup
    T ch16;     // X_velocity
    T ch17;     // Y_velocity
    T ch18;     // Z_velocity
    T ch19;     // epsilon
    T ch20;     // temp
    T ch21;     // temp
    T ch22;     // temp
    T ch23;     // temp
    T ch24;     // temp
    T ch25;     // temp
    T ch26;     // temp
    T ch27;     // temp
    T ch28;     // temp
    T ch29;     // temp
    T ch30;     // temp
};
}
#endif

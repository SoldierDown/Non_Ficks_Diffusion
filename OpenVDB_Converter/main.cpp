//!#####################################################################
//! \file main.cpp
//!#####################################################################

#include <fstream>
#include <iostream>
#include "OpenVDB_Converter.h"
#include <chrono>
#include <omp.h>

using namespace Nova;
using namespace std::chrono;

int main(int argc,char** argv)
{
    enum {d=3};
    typedef float T;
    std::string directory_name="/home/hertz/Nova/build/data_alienware/Implicit_Source_Smoke_3d_F_case_3_diff_0.010000_Fc_0.000000_tau_4.000000_bv_1.000000_sr_50.000000_Resolution_64x128x64";
    // File_Utilities::Create_Directory(directory_name+"/nodes_data");
    Nova::OpenVDB_Converter<T,d> converter(directory_name,1);
    for(int i=0;i<=500;++i) converter.Read_From_Frame(i);

    return 0;
}

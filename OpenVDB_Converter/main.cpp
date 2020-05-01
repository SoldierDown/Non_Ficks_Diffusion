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
    std::string directory_name="/home/hertz/Nova/build/data_alienware/test";
    Nova::OpenVDB_Converter<T,d> converter(directory_name,1);
    converter.Read_From_Frame(0);

    return 0;
}

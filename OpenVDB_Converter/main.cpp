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
    Nova::OpenVDB_Converter<T,d> converter;
    converter.Parse(argc,argv);
    boost::filesystem::path directory_path(converter.output_directory+"/converted_data");
    boost::filesystem::create_directory(directory_path);
    converter.Initialize();
    std::cout<<"first: "<<converter.first_frame<<", last: "<<converter.last_frame<<", step: "<<converter.step<<std::endl;
    for(int i=converter.first_frame;i<=converter.last_frame;i+=converter.step)
        converter.Convert_Frame(i);
    return 0;
}

//!#####################################################################
//! \file main.cpp
//!#####################################################################

#include <fstream>
#include <iostream>
#include "Snowflake_Thicken_Helper.h"
#include <chrono>
#include <omp.h>

using namespace Nova;
using namespace std::chrono;

int main(int argc,char** argv)
{
    enum {d=3};
    typedef float T;
    Nova::Snowflake_Thicken_Helper<T> thicken_helper;
    thicken_helper.Parse(argc,argv);
    boost::filesystem::path obj_path(thicken_helper.output_directory+"/obj_data");
    boost::filesystem::path vdb_path(thicken_helper.output_directory+"/vdb_data");
    boost::filesystem::create_directory(vdb_path);
    boost::filesystem::create_directory(obj_path);
    thicken_helper.Initialize();
    std::cout<<"first: "<<thicken_helper.first_frame<<", last: "<<thicken_helper.last_frame<<", step: "<<thicken_helper.step<<std::endl;
    for(int i=thicken_helper.first_frame;i<=thicken_helper.last_frame;i+=thicken_helper.step)
        thicken_helper.Convert_Frame(i);
    return 0;
}
//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <chrono>
#include "../MPM_Driver.h"
#include "Standard_Tests/Standard_Tests.h"
using namespace Nova;
using namespace std::chrono;

int main(int argc,char** argv)
{
    Log::cout.precision(20);
    high_resolution_clock::time_point tb = high_resolution_clock::now();
    enum {d=2};
    typedef float T;
    MPM_Example<T,d> *example=new Standard_Tests<T,d>();
    example->Parse(argc,argv);
    example->bbox=Range<T,d>(example->domain.max_corner,example->domain.min_corner);
    File_Utilities::Create_Directory(example->output_directory);
    File_Utilities::Create_Directory(example->output_directory+"/common");
    Log::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
        
    MPM_Driver<T,d> driver(*example);
    driver.Execute_Main_Program();
    
    Log::cout<<"ras: "<<example->ras_rt/example->ras_cnt<<std::endl;
    Log::cout<<"ras to voxel: "<<example->ras_vx_rt/example->ras_vx_cnt<<std::endl;
    Log::cout<<"explicit force: "<<example->explicit_force_rt/example->explicit_force_cnt<<std::endl;
    Log::cout<<"apply force: "<<example->apply_force_rt/example->apply_force_cnt<<std::endl;
    Log::cout<<"update x and v: "<<example->update_x_v_rt/example->update_x_v_cnt<<std::endl;

    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
	Log::cout<<"Total duration: "<<dur.count()<<std::endl;   
        
    delete example;

    return 0;
}

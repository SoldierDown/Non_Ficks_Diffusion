//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <nova/Tools/Utilities/Pthread_Queue.h>
#include "../Smoke_Driver.h"
#include "Standard_Tests/Standard_Tests.h"
#include <chrono>
using namespace std::chrono;
using namespace Nova;

namespace Nova{
int number_of_threads=0;
}

extern Pthread_Queue* pthread_queue;

int main(int argc,char** argv)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    enum {d=3};
    typedef float T;typedef Vector<T,d> TV;
    typedef Vector<int,d> T_INDEX;

    Smoke_Example<T,d> *example=new Standard_Tests<T,d>();
    example->Parse(argc,argv);

    if(number_of_threads) pthread_queue=new Pthread_Queue(number_of_threads);

    File_Utilities::Create_Directory(example->output_directory);
    File_Utilities::Create_Directory(example->output_directory+"/common");
    File_Utilities::Create_Directory(example->output_directory+"/density_data");
    File_Utilities::Create_Directory(example->output_directory+"/converted_data");
    Log::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    Smoke_Driver<T,d> driver(*example);
    driver.Execute_Main_Program();
    high_resolution_clock::time_point te=high_resolution_clock::now();
    T total_rt=duration_cast<duration<double>>(te-tb).count();
    int substeps=driver.substep_counter;
    Log::cout<<"Average: "<<std::endl;
    Log::cout<<"Total substeps: "<<substeps<<std::endl;
    Log::cout<<"full timestep: "<<total_rt/substeps<<std::endl;
    Log::cout<<"diffusion: "<<driver.diffision_rt/substeps<<std::endl;
    Log::cout<<"qc advection: "<<driver.qc_advection_rt/substeps<<std::endl;
    Log::cout<<"density advection: "<<driver.density_advection_rt/substeps<<std::endl;
    Log::cout<<"velocity advection: "<<driver.velocity_advection_rt/substeps<<std::endl;
    Log::cout<<"source modification: "<<driver.source_modification_rf/substeps<<std::endl;
    Log::cout<<"projection: "<<driver.projection_rt/substeps<<std::endl;
    delete example;

    return 0;
}

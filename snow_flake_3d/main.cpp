//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <nova/Tools/Utilities/Pthread_Queue.h>
#include "../SF_Driver.h"
#include "Standard_Tests/Standard_Tests.h"
using namespace Nova;

namespace Nova{
int number_of_threads=0;
}

extern Pthread_Queue* pthread_queue;

int main(int argc,char** argv)
{
    enum {d=3};
    typedef float T;typedef Vector<T,d> TV;
    typedef Vector<int,d> T_INDEX;

    SF_Example<T,d> *example=new Standard_Tests<T,d>();
    example->Parse(argc,argv);

    if(number_of_threads) pthread_queue=new Pthread_Queue(number_of_threads);

    File_Utilities::Create_Directory(example->output_directory);
    File_Utilities::Create_Directory(example->output_directory+"/common");
    File_Utilities::Create_Directory(example->output_directory+"/density_data");
    Log::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    SF_Driver<T,d> driver(*example);
    driver.Execute_Main_Program();
    int substeps=example->substep_counter;
    Log::cout<<"Average: "<<std::endl;
    Log::cout<<"Total substeps: "<<substeps<<std::endl;
    Log::cout<<"full timestep: "<<example->total_rt/substeps<<std::endl;
    Log::cout<<"Advect scalar: "<<example->advect_scalar_rt/substeps<<std::endl;
    Log::cout<<"Advect Q: "<<example->advect_Q_rt/substeps<<std::endl;
    Log::cout<<"Update density: "<<example->update_s_rt/substeps<<std::endl;
    Log::cout<<"Diffuse: "<<example->diffusion_rt/substeps<<std::endl;
    Log::cout<<"Update T: "<<example->update_t_rt/substeps<<std::endl;
    Log::cout<<"Update Qs: "<<example->update_qs_rt/substeps<<std::endl;
    Log::cout<<"Update Qt: "<<example->update_qt_rt/substeps<<std::endl;
    delete example;

    return 0;
}

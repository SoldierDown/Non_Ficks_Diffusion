//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <chrono>
#include "../MPM_Driver.h"
#include "Standard_Tests/Standard_Tests.h"
using namespace Nova;
using namespace std::chrono;

namespace Nova{
    int number_of_threads=0;
}
int main(int argc,char** argv)
{
    enum {d=3};
    typedef float T;
    MPM_Example<T,d> *example=new Standard_Tests<T,d>();
    example->Parse(argc,argv);
    example->bbox=Range<T,d>(example->domain.max_corner,example->domain.min_corner);
    File_Utilities::Create_Directory(example->output_directory);
    File_Utilities::Create_Directory(example->output_directory+"/common");
    File_Utilities::Create_Directory(example->output_directory+"/particle_indicator");
    Log::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    MPM_Driver<T,d> driver(*example);
    driver.Execute_Main_Program();
 
    int substeps=driver.substep_counter;
    Log::cout<<"Average: "<<std::endl;
    Log::cout<<"Total substeps: "<<substeps<<std::endl;
    Log::cout<<"Initialize SPGrid: "<<driver.init_spgrid_rt/substeps<<std::endl;
    Log::cout<<"Reset grid variables: "<<driver.reset_grid_var_rt/substeps<<std::endl;
    Log::cout<<"Update particle weights: "<<driver.update_p_weights_rt/substeps<<std::endl;
    Log::cout<<"Group Particles: "<<driver.group_p_rt/substeps<<std::endl;
    Log::cout<<"Rasterize: "<<driver.rasterize_rt/substeps<<std::endl;
    Log::cout<<"Process waiting particles: "<<driver.procee_waiting_p_rt/substeps<<std::endl;
    Log::cout<<"Diffusion: "<<driver.diffusion_rt/substeps<<std::endl;
    Log::cout<<"Update constitutive model state: "<<driver.update_constitutive_model_state_rt/substeps<<std::endl;
    Log::cout<<"Update particle velocity and position: "<<driver.update_p_v_x_rt/substeps<<std::endl;
    Log::cout<<"Populate simulated particles: "<<driver.pop_sim_p_rt/substeps<<std::endl;
    Log::cout<<"Total iterations: "<<example->iteration_counter<<std::endl;
    Log::cout<<"average iterations: "<<(T)example->iteration_counter/(T)substeps<<std::endl;
    Log::cout<<"Total: "<<driver.total_rt/substeps<<std::endl;

    delete example;

    return 0;
}

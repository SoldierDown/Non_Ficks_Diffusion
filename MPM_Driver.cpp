//!#####################################################################
//! \file MPM_Driver.cpp
//!#####################################################################
#include <chrono>
#include "MPM_Driver.h"
using namespace std::chrono;
using namespace Nova;
//######################################################################
// Constructor
//######################################################################
template<class T,int d> MPM_Driver<T,d>::
MPM_Driver(MPM_Example<T,d>& example_input)
    :Base(example_input),example(example_input)
{}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Initialize()
{
    if(!example.restart){
        example.Initialize();
        example.Reset_Grid_Based_Variables();
        example.Update_Particle_Weights();
        example.Group_Particles();
        example.Rasterize_Voxels();
        example.Rasterize();
        example.Process_Waiting_Particles();
    }
    else example.Read_Output_Files(example.restart_frame);
}
//######################################################################
// Test
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Test()
{
    example.Initialize();
    example.Test();
}
//######################################################################
// Execute_Main_Program
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Execute_Main_Program() 
{
    // Test();
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//######################################################################
// Advance_To_Target_Time
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Advance_To_Target_Time(const T target_time)
{
    T min_dt=(T)1e-6;
    T max_dt=(T).005;
    T cfl=example.cfl;
    T dx_min=example.hierarchy->Lattice(0).dX(0);
    bool done=false;
    example.Populate_Simulated_Particles(); // add only valid particles to array
    for(int substep=1;!done;substep++){
        T max_v=example.Max_Particle_Velocity();
        Log::Scope scope("SUBSTEP","substep "+std::to_string(substep));
        T dt=std::max(min_dt,std::min(max_dt,cfl*dx_min/std::max(max_v,(T)1e-2)));
        if(target_time-time<dt*1.001){
            dt=target_time-time;
            done=true;
        }
        else if(target_time-time<(T)2.*dt){
            dt=(target_time-time)*(T).5;
        }
        Advance_Step(dt);
        Log::cout<<"dt: "<<dt<<std::endl;
        // if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;
        }
}
//######################################################################
// Simulate_To_Frame
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Simulate_To_Frame(const int target_frame)
{
    example.frame_title="Frame "+std::to_string(example.current_frame);
    if(!example.restart) Write_Output_Files(example.current_frame);

    while(example.current_frame<target_frame){
        Log::Scope scope("FRAME","Frame "+std::to_string(++example.current_frame));

        Advance_To_Target_Time(example.Time_At_Frame(example.current_frame));

        example.frame_title="Frame "+std::to_string(example.current_frame);
        Write_Output_Files(++example.output_number);

        *(example.output)<<"TIME = "<<time<<std::endl;}
}
//######################################################################
// Advance_Step
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Advance_Step(const T dt)
{
    example.Initialize_SPGrid();
    example.Reset_Grid_Based_Variables();
    example.Add_Fluid_Source();
    example.Update_Particle_Weights();
    example.Group_Particles();
    example.Rasterize_Voxels();
    example.Rasterize();
    example.Process_Waiting_Particles();
    example.Update_Constitutive_Model_State();
    example.Update_Particle_Velocities_And_Positions(dt);
}
//######################################################################
template class Nova::MPM_Driver<float,2>;
template class Nova::MPM_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::MPM_Driver<double,2>;
template class Nova::MPM_Driver<double,3>;
#endif

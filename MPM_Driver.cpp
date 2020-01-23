//!#####################################################################
//! \file MPM_Driver.cpp
//!#####################################################################
#include <chrono>
#include "MPM_Driver.h"
using namespace Nova;
using namespace std::chrono;
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
    const T dt=1e-5;
    if(!example.restart){
        example.Initialize();
        example.Reset_Grid_Based_Variables();
        example.Update_Particle_Weights();
        example.Group_Particles();
        example.Rasterize();
        if(example.FICKS) example.Ficks_Diffusion(dt);
        else example.Non_Ficks_Diffusion(dt);
        example.Estimate_Particle_Volumes();
        example.Update_Constitutive_Model_State();
        example.Update_Particle_Velocities_And_Positions(dt);
        time=dt;
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
    example.Reset_Grid_Based_Variables();
    // example.Test();
}
//######################################################################
// Execute_Main_Program
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Execute_Main_Program() 
{
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
    T max_dt=(T).001;
    T cfl=example.cfl;
    T dx_min=example.hierarchy->Lattice(0).dX(0);
    bool done=false;
    for(int substep=1;!done;substep++){
        example.Populate_Simulated_Particles(); // add only valid particles to array
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
        Log::cout<<"dt: "<<dt<<std::endl;
        Advance_Step(dt);
        if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;
        Log::cout<<"time: "<<time<<std::endl;}
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
    bool SHOW_RUNNING_TIME=false;
    high_resolution_clock::time_point tb1 = high_resolution_clock::now();
    example.Initialize_SPGrid();
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te1 = high_resolution_clock::now();
	duration<double> d1 = duration_cast<duration<double>>(te1 - tb1);
	std::printf("Initialize SPGrid duration: %f\n", d1.count());}

    high_resolution_clock::time_point tb2 = high_resolution_clock::now();
    example.Reset_Grid_Based_Variables();
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te2 = high_resolution_clock::now();
	duration<double> d2 = duration_cast<duration<double>>(te2 - tb2);
	std::printf("Reset Variables duration: %f\n", d2.count());}

    high_resolution_clock::time_point tb8 = high_resolution_clock::now();
    example.Update_Particle_Weights();
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te8 = high_resolution_clock::now();
	duration<double> d8 = duration_cast<duration<double>>(te8 - tb8);
	std::printf("Update Particle Weights: %f\n", d8.count());}

    high_resolution_clock::time_point tb3 = high_resolution_clock::now();
    example.Group_Particles();
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te3 = high_resolution_clock::now();
	duration<double> d3 = duration_cast<duration<double>>(te3 - tb3);
	std::printf("Group Particles duration: %f\n", d3.count());}

    high_resolution_clock::time_point tb4 = high_resolution_clock::now();
    example.Rasterize();
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te4 = high_resolution_clock::now();
	duration<double> d4 = duration_cast<duration<double>>(te4 - tb4);
	std::printf("Rasterize duration: %f\n", d4.count());}

    high_resolution_clock::time_point tb9 = high_resolution_clock::now();
    if(example.FICKS) example.Ficks_Diffusion(dt);
    else example.Non_Ficks_Diffusion(dt);
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te9 = high_resolution_clock::now();
	duration<double> d9 = duration_cast<duration<double>>(te9 - tb9);
	std::printf("Diffusion duration: %f\n", d9.count());}

    high_resolution_clock::time_point tb5 = high_resolution_clock::now();
    example.Update_Constitutive_Model_State();
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te5 = high_resolution_clock::now();
	duration<double> d5 = duration_cast<duration<double>>(te5 - tb5);
	std::printf("Update Constitutive Model State duration: %f\n", d5.count());}

    high_resolution_clock::time_point tb6 = high_resolution_clock::now();
    example.Update_Particle_Velocities_And_Positions(dt);
    if (SHOW_RUNNING_TIME){high_resolution_clock::time_point te6 = high_resolution_clock::now();
	duration<double> d6 = duration_cast<duration<double>>(te6 - tb6);
	std::printf("Update Particle Velocities And Positions duration: %f\n", d6.count());}
}
//######################################################################
template class Nova::MPM_Driver<float,2>;
template class Nova::MPM_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::MPM_Driver<double,2>;
template class Nova::MPM_Driver<double,3>;
#endif

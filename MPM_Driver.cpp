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
    substep_counter=0;
    total_rt=(T)0.; init_spgrid_rt=(T)0.; reset_grid_var_rt=(T)0.; update_p_weights_rt=(T)0.;
    group_p_rt=(T)0.; rasterize_rt=(T)0.; procee_waiting_p_rt=(T)0.; diffusion_rt=(T)0.;
    update_constitutive_model_state_rt=(T)0.; update_p_v_x_rt=(T)0.; pop_sim_p_rt=(T)0.;
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
    example.Reset_Grid_Based_Variables();
}
//######################################################################
// Execute_Main_Program
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Execute_Main_Program() 
{
    // example.Test();
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//######################################################################
// Advance_To_Target_Time
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Advance_To_Target_Time(const T target_time)
{
    Log::cout<<"min dt: "<<example.min_dt<<", max_dt: "<<example.max_dt<<std::endl;
    T min_dt=example.min_dt;
    T max_dt=example.max_dt;
    T cfl=example.cfl;
    T dx_min=example.mpm_hierarchy->Lattice(0).dX(0);
    bool done=false;
    for(int substep=1;!done;substep++){
        substep_counter++;
        high_resolution_clock::time_point tb=high_resolution_clock::now();
        example.Populate_Simulated_Particles(); // add only valid particles to array
        high_resolution_clock::time_point te=high_resolution_clock::now();
        pop_sim_p_rt+=duration_cast<duration<T>>(te-tb).count();

        Log::cout<<"number of particles: "<<example.simulated_particles.size()<<std::endl;
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
        Log::cout<<"time: "<<time<<std::endl;
        te=high_resolution_clock::now();
        total_rt+=duration_cast<duration<T>>(te-tb).count();
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
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    example.Initialize_SPGrid();
    high_resolution_clock::time_point te=high_resolution_clock::now();
    init_spgrid_rt+=duration_cast<duration<T>>(te-tb).count();


    tb=high_resolution_clock::now();
    example.Reset_Grid_Based_Variables();
    te=high_resolution_clock::now();
    reset_grid_var_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Update_Particle_Weights();
    te=high_resolution_clock::now();
    update_p_weights_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Group_Particles();
    te=high_resolution_clock::now();
    group_p_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Rasterize_Voxels();
    te=high_resolution_clock::now();
    rasterize_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Rasterize();
    te=high_resolution_clock::now();
    rasterize_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Process_Waiting_Particles();
    te=high_resolution_clock::now();
    procee_waiting_p_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    if(example.FICKS) example.Ficks_Diffusion(dt);
    else example.Non_Ficks_Diffusion(dt);
    te=high_resolution_clock::now();
    diffusion_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Update_Constitutive_Model_State();
    te=high_resolution_clock::now();
    update_constitutive_model_state_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Update_Particle_Velocities_And_Positions(dt);
    te=high_resolution_clock::now();
    update_p_v_x_rt+=duration_cast<duration<T>>(te-tb).count();
}
//######################################################################
template class Nova::MPM_Driver<float,2>;
template class Nova::MPM_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::MPM_Driver<double,2>;
template class Nova::MPM_Driver<double,3>;
#endif

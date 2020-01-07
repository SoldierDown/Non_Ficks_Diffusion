//!#####################################################################
//! \file MPM_Driver.cpp
//!#####################################################################
#include "MPM_Driver.h"
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
        example.Rasterize();
        example.Estimate_Particle_Volumes();
        example.Update_Constitutive_Model_State();
        example.Update_Particle_Velocities_And_Positions(1e-5);
    }
    else example.Read_Output_Files(example.restart_frame);
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
    bool done=false;
    T min_dt=(T)1e-6;
    T max_dt=(T).005;
    T cfl=(T).1;
    T max_v;
    T dx_min=example.hierarchy->Lattice(0).dX(0);
    for(int substep=1;!done;substep++){
        example.Populate_Simulated_Particles(); // add only valid particles to array
        Log::Scope scope("SUBSTEP","substep "+std::to_string(substep));
        //T dt=Compute_Dt(time,target_time);
        max_v=example.Max_Particle_Velocity();
        T dt=std::max(min_dt,std::min(max_dt,cfl*dx_min/std::max(max_v,(T)1e-2)));
        dt=(T)1e-3;
        Example<T,d>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
        Log::cout<<"dt: "<<dt<<std::endl;
        Advance_Step(dt);
        if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;}
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
    // std::cout<<"1"<<std::endl;
    example.Initialize_SPGrid();
    // std::cout<<"2"<<std::endl;
    example.Reset_Grid_Based_Variables();
    // std::cout<<"3"<<std::endl;
    example.Rasterize();
    // std::cout<<"4"<<std::endl;
    example.Update_Constitutive_Model_State();
    // std::cout<<"5"<<std::endl;
    example.Update_Particle_Velocities_And_Positions(dt);
    // std::cout<<"6"<<std::endl;
}
//######################################################################
template class Nova::MPM_Driver<float,2>;
template class Nova::MPM_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::MPM_Driver<double,2>;
template class Nova::MPM_Driver<double,3>;
#endif

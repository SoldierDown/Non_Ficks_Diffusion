//!#####################################################################
//! \file Smoke_Driver.cpp
//!#####################################################################
#include <chrono>
#include "Smoke_Driver.h"
using namespace std::chrono;
using namespace Nova;
//######################################################################
// Constructor
//######################################################################
template<class T,int d> Smoke_Driver<T,d>::
Smoke_Driver(Smoke_Example<T,d>& example_input)
    :Base(example_input),example(example_input)
{}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Smoke_Driver<T,d>::
Initialize()
{
    substep_counter=0;
    diffision_rt=(T)0.; qc_advection_rt=(T)0.; density_advection_rt=(T)0.;
    velocity_advection_rt=(T)0.; source_modification_rf=(T)0.; projection_rt=(T)0.;
    Base::Initialize();

    example.Log_Parameters();
    example.Initialize_Sources(example.test_number);

    if(!example.restart) example.Initialize();
    else example.Read_Output_Files(example.restart_frame);

    example.Initialize_Velocity_Field();
    // divergence free
    if(!example.uvf) example.Project();
}
//######################################################################
// Advance_One_Time_Step_Explicit_Part
//######################################################################
template<class T,int d> void Smoke_Driver<T,d>::
Advance_One_Time_Step_Explicit_Part(const T dt,const T time)
{
    if(!example.nd) {
        high_resolution_clock::time_point tb=high_resolution_clock::now();
        if(example.FICKS) example.Ficks_Diffusion(dt);
        else example.Non_Ficks_Diffusion(dt);
        high_resolution_clock::time_point te=high_resolution_clock::now();
	    diffision_rt+=duration_cast<duration<T>>(te-tb).count();
        tb=high_resolution_clock::now();
        if(!example.FICKS) example.Advect_Face_Vector(dt);
        te=high_resolution_clock::now();
	    qc_advection_rt+=duration_cast<duration<T>>(te-tb).count();
    example.Backup_Density();}
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    example.Advect_Density(dt);
    high_resolution_clock::time_point te=high_resolution_clock::now();
    density_advection_rt+=duration_cast<duration<T>>(te-tb).count();
    tb=high_resolution_clock::now();
    if(example.const_density_source) example.Modify_Density_With_Sources();
    else example.Add_Source(dt);
    te=high_resolution_clock::now();
    source_modification_rf+=duration_cast<duration<T>>(te-tb).count();
    // convect
    tb=high_resolution_clock::now();
    if(!example.uvf) example.Advect_Face_Velocities(dt);
    te=high_resolution_clock::now();
    velocity_advection_rt+=duration_cast<duration<T>>(te-tb).count();

}
//######################################################################
// Advance_One_Time_Step_Implicit_Part
//######################################################################
template<class T,int d> void Smoke_Driver<T,d>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    if(!example.uvf) example.Project();
    high_resolution_clock::time_point te=high_resolution_clock::now();
    projection_rt+=duration_cast<duration<T>>(te-tb).count();
}
//######################################################################
// Advance_To_Target_Time
//######################################################################
template<class T,int d> void Smoke_Driver<T,d>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        Log::Scope scope("SUBSTEP","substep "+std::to_string(substep));
        substep_counter++;
        T dt=Compute_Dt(time,target_time);
        if(example.explicit_diffusion) dt/=(T)100.;
        Example<T,d>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
        Advance_One_Time_Step_Explicit_Part(dt,time);
        Advance_One_Time_Step_Implicit_Part(dt,time);
        Log::cout<<"dt: "<<dt<<std::endl;
        if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;}
}
//######################################################################
// Simulate_To_Frame
//######################################################################
template<class T,int d> void Smoke_Driver<T,d>::
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
template class Nova::Smoke_Driver<float,2>;
template class Nova::Smoke_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Smoke_Driver<double,2>;
template class Nova::Smoke_Driver<double,3>;
#endif

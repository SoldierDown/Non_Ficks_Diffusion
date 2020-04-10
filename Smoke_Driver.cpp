//!#####################################################################
//! \file Smoke_Driver.cpp
//!#####################################################################
#include "Smoke_Driver.h"
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
    Base::Initialize();

    example.Log_Parameters();
    example.Initialize_Sources();

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
    // scalar advance
    example.Add_Source(dt);
    example.Backup_Density();
    example.Advect_Density(dt);
    // example.Modify_Density_With_Sources();
    example.Diffuse_Density(dt);
    example.Backup_Density();
    // convect
    if(!example.uvf) example.Advect_Face_Velocities(dt);
}
//######################################################################
// Advance_One_Time_Step_Implicit_Part
//######################################################################
template<class T,int d> void Smoke_Driver<T,d>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time)
{
    if(!example.uvf) example.Project();
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
        T dt=Compute_Dt(time,target_time);
        if(example.explicit_diffusion) dt/=(T)100.;
        Example<T,d>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
        Advance_One_Time_Step_Explicit_Part(dt,time);
        Advance_One_Time_Step_Implicit_Part(dt,time);
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

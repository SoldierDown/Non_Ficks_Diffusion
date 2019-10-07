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
    Base::Initialize();

    if(!example.restart){example.Initialize();
        example.Rasterize();
    }
    else example.Read_Output_Files(example.restart_frame);
}
//######################################################################
// Advance_To_Target_Time
//######################################################################
template<class T,int d> void MPM_Driver<T,d>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        Log::Scope scope("SUBSTEP","substep "+std::to_string(substep));
        T dt=Compute_Dt(time,target_time);
        Example<T,d>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
        Log::cout<<"Dt: "<<dt<<std::endl;

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
template class Nova::MPM_Driver<float,2>;
template class Nova::MPM_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::MPM_Driver<double,2>;
template class Nova::MPM_Driver<double,3>;
#endif

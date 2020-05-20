//!#####################################################################
//! \file SF_Driver.cpp
//!#####################################################################
#include "SF_Driver.h"
#include <chrono>
using namespace std::chrono;
using namespace Nova;
//######################################################################
// Constructor
//######################################################################
template<class T,int d> SF_Driver<T,d>::
SF_Driver(SF_Example<T,d>& example_input)
    :Base(example_input),example(example_input)
{
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void SF_Driver<T,d>::
Initialize()
{
    Base::Initialize();

    example.Log_Parameters();
    example.Initialize_Sources();

    if(!example.restart) example.Initialize();
    else example.Read_Output_Files(example.restart_frame);

    example.Initialize_Velocity_Field();
    example.Backup();
}
//######################################################################
// Advance_One_Time_Step_Explicit_Part
//######################################################################
template<class T,int d> void SF_Driver<T,d>::
Advance_One_Time_Step_Explicit_Part(const T dt,const T time)
{
    // scalar advance
    // example.Add_Source(dt);
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    example.Advect_Scalar_Field(dt);
    high_resolution_clock::time_point te=high_resolution_clock::now();
    example.advect_scalar_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    if(!example.FICKS) example.Advect_Face_Vector_Field(dt);
    te=high_resolution_clock::now();
    example.advect_Q_rt+=duration_cast<duration<T>>(te-tb).count();

    if(example.FICKS){Log::cout<<"Fick's"<<std::endl;
    tb=high_resolution_clock::now();
    example.Update_Density(dt); 
    te=high_resolution_clock::now();
    example.update_s_rt+=duration_cast<duration<T>>(te-tb).count();
    
    tb=high_resolution_clock::now();
    example.Update_Temperature(dt);
    te=high_resolution_clock::now();
    example.update_t_rt+=duration_cast<duration<T>>(te-tb).count();}
    else{Log::cout<<"Non-Fick's"<<std::endl;
    tb=high_resolution_clock::now();
    example.Update_Density(dt);
    te=high_resolution_clock::now();
    example.update_s_rt+=duration_cast<duration<T>>(te-tb).count();
    
    tb=high_resolution_clock::now();
    example.Update_Face_Qsc(dt);
    te=high_resolution_clock::now();
    example.update_qs_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Update_Temperature(dt);
    te=high_resolution_clock::now();
    example.update_t_rt+=duration_cast<duration<T>>(te-tb).count();

    tb=high_resolution_clock::now();
    example.Update_Face_Qtc(dt);
    te=high_resolution_clock::now();
    example.update_qt_rt+=duration_cast<duration<T>>(te-tb).count();
    }
    example.Backup();

    // convect
    // example.Advect_Face_Velocities(dt);
}
//######################################################################
// Advance_One_Time_Step_Implicit_Part
//######################################################################
template<class T,int d> void SF_Driver<T,d>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time)
{
    // example.Project(dt);
}
//######################################################################
// Advance_To_Target_Time
//######################################################################
template<class T,int d> void SF_Driver<T,d>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        example.substep_counter++;
        high_resolution_clock::time_point tb=high_resolution_clock::now();
        Log::Scope scope("SUBSTEP","substep "+std::to_string(substep));
        T dt=Compute_Dt(time,target_time);
        // dt=std::max(min_dt,std::min(max_dt,dt));
        // if(example.explicit_diffusion) dt/=(T)10.;
        dt=example.const_time_step;
        Example<T,d>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
        
        Advance_One_Time_Step_Explicit_Part(dt,time);
        Advance_One_Time_Step_Implicit_Part(dt,time);
        Log::cout<<"dt: "<<dt<<std::endl;
        // done=true;
        if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;
        high_resolution_clock::time_point te=high_resolution_clock::now();
        example.total_rt+=duration_cast<duration<T>>(te-tb).count();
    }
}
//######################################################################
// Simulate_To_Frame
//######################################################################
template<class T,int d> void SF_Driver<T,d>::
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
template class Nova::SF_Driver<float,2>;
template class Nova::SF_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::SF_Driver<double,2>;
template class Nova::SF_Driver<double,3>;
#endif

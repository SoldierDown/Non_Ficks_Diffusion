//!#####################################################################
//! \file Smoke_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Advection/Grid_Hierarchy_Advection.h>
#include <nova/Tools/Grids/Grid_Iterator_Face.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include "Apply_Pressure.h"
#include "Boundary_Value_Initializer.h"
#include "Uniform_Velocity_Field_Initializer.h"
#include "Velocity_Field_Traverser.h"
#include "Boundary_Condition_Helper.h"
#include "Poisson_Solver/Poisson_CG_System.h"
#include "Compute_Time_Step.h"
#include "Density_Modifier.h"
#include "Source_Velocity_Setup.h"
#include "Source_Adder.h"
#include "Initialize_Dirichlet_Cells.h"
#include "Poisson_Solver/Multigrid_Data.h"
#include "Smoke_Example.h"
#include "Ficks_Solver/Ficks_CG_System.h"
#include "Non_Ficks_Solver/Non_Ficks_CG_System.h"
#include "Ficks_RHS_Helper.h"
#include "Non_Ficks_RHS_Helper.h"
#include "Non_Ficks_Smoke_Density_Explicit_Update_Helper.h"
#include "Clamp_Helper.h"
#include "Lap_Calculator.h"
#include "Density_Backup_Helper.h"
#include "Explicit_Face_Qc_Updater.h"
#include "Implicit_Face_Qc_Updater.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Advection_Helper.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Averaging_Helper.h"
#include "Density_Clamp_Helper.h"
#include "Flip_Helper.h"
#include "Write_To_File_Helper.h"
#include <omp.h>
using namespace Nova;
namespace Nova{
extern int number_of_threads;
}
//######################################################################
// Constructor
//######################################################################
template<class T,int d> Smoke_Example<T,d>::
Smoke_Example()
    :Base(),hierarchy(nullptr),rasterizer(nullptr)
{
    face_velocity_channels(0)           = &Struct_type::ch0;
    face_velocity_channels(1)           = &Struct_type::ch1;
    if(d==3) face_velocity_channels(2)  = &Struct_type::ch2;
    face_qc_channels(0)                 = &Struct_type::ch3;
    face_qc_channels(1)                 = &Struct_type::ch4;
    if(d==3) face_qc_channels(2)        = &Struct_type::ch5;
    density_channel                     = &Struct_type::ch6;
    density_backup_channel              = &Struct_type::ch7;
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Initialize()
{
    Initialize_SPGrid();
    Initialize_Fluid_State(test_number);
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Initialize_SPGrid()
{
    Log::Scope scope("Initialize_SPGrid");
    Initialize_Rasterizer(test_number);
    for(Grid_Hierarchy_Iterator<d,Hierarchy_Rasterizer> iterator(hierarchy->Lattice(levels-1).Cell_Indices(),levels-1,*rasterizer);iterator.Valid();iterator.Next());
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Cells(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Valid_Faces(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Active_Faces(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Active_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Shared_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_T_Junction_Nodes(*hierarchy);
    Initialize_Dirichlet_Cells<Struct_type,T,d>(*hierarchy,domain_walls);
    //Set_Neumann_Faces_Inside_Sources();
    hierarchy->Update_Block_Offsets();
    hierarchy->Initialize_Red_Black_Partition(2*number_of_threads);
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Limit_Dt(T& dt,const T time)
{
    T dt_convection=(T)0.;

    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);

    for(int level=0;level<levels;++level)
        Compute_Time_Step<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,
                                           other_face_offsets,level,dt_convection);

    if(dt_convection>(T)1e-5) dt=cfl/dt_convection;
    Log::cout<<"Time Step: "<<dt<<std::endl;
}
//######################################################################
// Advect_Density
//######################################################################
template    <class T,int d> void Smoke_Example<T,d>::
Advect_Density(const T dt)
{
    Channel_Vector cell_velocity_channels;
    cell_velocity_channels(0)               = &Struct_type::ch8;
    cell_velocity_channels(1)               = &Struct_type::ch9;
    if(d==3) cell_velocity_channels(2)      = &Struct_type::ch10;
    T Struct_type::* temp_channel           = &Struct_type::ch11;
    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
    Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Cells(*hierarchy,hierarchy->Allocator(0),hierarchy->Blocks(0),face_velocity_channels,cell_velocity_channels,other_face_offsets);
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Density(*hierarchy,cell_velocity_channels,density_channel,temp_channel,dt);
}
//######################################################################
// Diffuse_Density
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Diffuse_Density(const T dt)
{
    if(FICKS) Ficks_Diffusion(dt);
    else Non_Ficks_Diffusion(dt); 
}
//######################################################################
// Backup_Density
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Backup_Density()
{
    for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel);
    for(int level=0;level<levels;++level) SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,
                                                                                density_backup_channel,(unsigned)Cell_Type_Interior);
}
//######################################################################
// Ficks_Diffusion
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Ficks_Diffusion(const T dt)
{
    using Multigrid_struct_type 	= Multigrid_Data<T>;
    Log::cout<<"Fick's Diffusion"<<std::endl;
	const Grid<T,d>& grid=hierarchy->Lattice(0);
    const T one_over_dx2=Nova_Utilities::Sqr(grid.one_over_dX(0));
    const T a=diff_coeff*dt*one_over_dx2; const T two_d_a_plus_one=(T)2*d*a+(T)1.;
	if(!explicit_diffusion){
	    Ficks_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,diff_coeff*dt,3,1,200);
        T Struct_type::* q_channel              = &Struct_type::ch8;
        T Struct_type::* r_channel              = &Struct_type::ch9;
        T Struct_type::* s_channel              = &Struct_type::ch10;
        T Struct_type::* k_channel              = &Struct_type::ch11;
        T Struct_type::* z_channel              = &Struct_type::ch11;
        T Struct_type::* b_channel              = &Struct_type::ch12;
    
    // clear all channels
    for(int level=0;level<levels;++level){
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),b_channel);}

	    for(int level=0;level<levels;++level) Ficks_RHS_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,b_channel,a);
	    CG_Vector<Struct_type,T,d> x_V(*hierarchy,density_channel),b_V(*hierarchy,b_channel),q_V(*hierarchy,q_channel),
                                    s_V(*hierarchy,s_channel),r_V(*hierarchy,r_channel),k_V(*hierarchy,k_channel),z_V(*hierarchy,z_channel);   
	    Conjugate_Gradient<T> cg;
        cg_system.Multiply(x_V,r_V);
        r_V-=b_V;
        const T b_norm=cg_system.Convergence_Norm(r_V);
        Log::cout<<"Norm: "<<b_norm<<std::endl;
        // cg.print_residuals=true;
        // cg.print_diagnostics=true;
        cg.restart_iterations=cg_restart_iterations;
        const T tolerance=std::max((T)1e-6*b_norm,(T)1e-6);
        cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);
        for(int level=0;level<levels;++level) Clamp_Heler<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);}
    else{
        T Struct_type::* lap_density_channel    = &Struct_type::ch8;
        for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),lap_density_channel);
        for(int level=0;level<levels;++level) Lap_Calculator<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel,lap_density_channel,one_over_dx2);        
        for(int level=0;level<levels;++level) SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dt,lap_density_channel,density_backup_channel,density_channel,Cell_Type_Interior);
        for(int level=0;level<levels;++level) Clamp_Heler<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);}
        Log::cout<<"Fick's Diffusion finished"<<std::endl;
}
//######################################################################
// Non_Ficks_Diffusion
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Non_Ficks_Diffusion(const T dt)
{
    using Multigrid_struct_type             = Multigrid_Data<T>;
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;
    Log::cout<<"Non-Fick's Diffusion"<<std::endl;
    const Grid<T,d>& grid=hierarchy->Lattice(0);
	const T dx2=Nova_Utilities::Sqr(grid.dX(0));        const T one_over_dx2=(T)1./dx2;
    Log::cout<<"diff_coeff: "<<diff_coeff<<", fc: "<<Fc<<", tau: "<<tau<<", bv: "<<bv<<", sr: "<<source_rate<<std::endl;
    if(explicit_diffusion){
        const T coeff3=dt*diff_coeff*Fc;                    const T coeff4=-dt;
        const T coeff5=-dt/tau;                             const T coeff6=-diff_coeff*dt*(1-Fc)/tau;
        T Struct_type::* div_qc_channel                     = &Struct_type::ch8;
        for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_qc_channel);
        // compute div(Qc) at time n
        Hierarchy_Projection::Compute_Divergence(*hierarchy,face_qc_channels,div_qc_channel); 
	    for(int level=0;level<levels;++level) Flip_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_qc_channel);
        
        T Struct_type::* lap_density_channel                = &Struct_type::ch9;
        for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),lap_density_channel);
        for(int level=0;level<levels;++level) Lap_Calculator<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel,lap_density_channel,one_over_dx2);        
        for(int level=0;level<levels;++level) 
            Non_Ficks_Smoke_Density_Explicit_Update_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,lap_density_channel,div_qc_channel,coeff3,coeff4);
        for(int level=0;level<levels;++level) Density_Clamp_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);
        
        // Advect face vector by axis
        Channel_Vector interpolated_face_velocity_channels;
        interpolated_face_velocity_channels(0)              = &Struct_type::ch8;
        interpolated_face_velocity_channels(1)              = &Struct_type::ch9;
        if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch10; 
        Channel_Vector face_qc_backup_channels;
        face_qc_backup_channels(0)                          = &Struct_type::ch11;
        face_qc_backup_channels(1)                          = &Struct_type::ch12;
        if(d==3) face_qc_backup_channels(2)                 = &Struct_type::ch13; 
        T Struct_type::* temp_channel                       = &Struct_type::ch14;  
        // Clear
        for(int level=0;level<levels;++level) {for(int v=0;v<d;++v) {
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),interpolated_face_velocity_channels(v));
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qc_backup_channels(v));}
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),temp_channel);} 
        
        for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis) 
            SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qc_channels(axis),
                                                face_qc_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));
        Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Vector(*hierarchy,face_qc_channels,face_velocity_channels,interpolated_face_velocity_channels,temp_channel,dt);
        // update Qc
        for(int level=0;level<levels;++level)
            Explicit_Face_Qc_Updater<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_qc_channels,face_qc_channels,density_backup_channel,coeff5,coeff6,level);}
    else{
        const T coeff1=dt*diff_coeff*(Fc*tau+dt)/(dt+tau);  const T coeff2=dt*tau/(dt+tau);                    
        const T coeff5=tau/(tau+dt);                        const T coeff6=-diff_coeff*dt*(1-Fc)/(tau+dt);        
        // compute div(Qc^n)
        T Struct_type::* div_qc_channel                 = &Struct_type::ch8;
        for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_qc_channel);
        Hierarchy_Projection::Compute_Divergence(*hierarchy,face_qc_channels,div_qc_channel);
        for(int level=0;level<levels;++level) Flip_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_qc_channel);
        
        T Struct_type::* q_channel                      = &Struct_type::ch9;
        T Struct_type::* r_channel                      = &Struct_type::ch10;
        T Struct_type::* s_channel                      = &Struct_type::ch11;
        T Struct_type::* k_channel                      = &Struct_type::ch11;
        T Struct_type::* z_channel                      = &Struct_type::ch12;
        T Struct_type::* b_channel                      = &Struct_type::ch13;
        // clear all channels
        for(int level=0;level<levels;++level){
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),b_channel);}
        Non_Ficks_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,coeff1,3,1,200);
        for(int level=0;level<levels;++level) Non_Ficks_RHS_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,div_qc_channel,b_channel,
                                                                                    coeff1*one_over_dx2,coeff2);
        for(int level=0;level<levels;++level) SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel,
                                                                                density_channel,(unsigned)Cell_Type_Interior);
        CG_Vector<Struct_type,T,d> x_V(*hierarchy,density_channel),b_V(*hierarchy,b_channel),q_V(*hierarchy,q_channel),
                                                    s_V(*hierarchy,s_channel),r_V(*hierarchy,r_channel),k_V(*hierarchy,z_channel),z_V(*hierarchy,z_channel);   
        Conjugate_Gradient<T> cg;
        cg_system.Multiply(x_V,r_V);
        r_V-=b_V;
        const T b_norm=cg_system.Convergence_Norm(r_V);
        Log::cout<<"Norm: "<<b_norm<<std::endl;
        cg.print_residuals=true;
        cg.print_diagnostics=true;
        cg.restart_iterations=cg_restart_iterations;
        const T tolerance=std::max((T)1e-6*b_norm,(T)1e-6);
        cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);
        for(int level=0;level<levels;++level) Density_Clamp_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);

        // Advect face vector by axis
        Channel_Vector interpolated_face_velocity_channels;
        interpolated_face_velocity_channels(0)              = &Struct_type::ch8;
        interpolated_face_velocity_channels(1)              = &Struct_type::ch9;
        if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch10; 
        Channel_Vector face_qc_backup_channels;
        face_qc_backup_channels(0)                          = &Struct_type::ch11;
        face_qc_backup_channels(1)                          = &Struct_type::ch12;
        if(d==3) face_qc_backup_channels(2)                 = &Struct_type::ch13; 
        T Struct_type::* temp_channel                       = &Struct_type::ch14;  

        // Clear
        for(int level=0;level<levels;++level) {for(int v=0;v<d;++v) {
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),interpolated_face_velocity_channels(v));
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qc_backup_channels(v));}
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),temp_channel);} 
            
        Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Vector(*hierarchy,face_qc_channels,face_velocity_channels,interpolated_face_velocity_channels,temp_channel,dt);
        // Back up face qc
        for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis) 
            SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qc_channels(axis),
                                                face_qc_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));
        // update Qc
        for(int level=0;level<levels;++level)
            Implicit_Face_Qc_Updater<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_qc_channels,face_qc_backup_channels,density_channel,coeff5,coeff6,level);
    }
}

//######################################################################
// Modify_Density_With_Sources
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Modify_Density_With_Sources()
{
    for(int level=0;level<levels;++level)
        Density_Modifier<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,density_sources,level);
}
//######################################################################
// Add_Source
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Add_Source(const T dt)
{
    for(int level=0;level<levels;++level)
        Source_Adder<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,density_sources,source_rate,dt,level);
}
//######################################################################
// Advect_Face_Velocities
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Advect_Face_Velocities(const T dt)
{
    Channel_Vector interpolated_face_velocity_channels;
    interpolated_face_velocity_channels(0)              = &Struct_type::ch8;
    interpolated_face_velocity_channels(1)              = &Struct_type::ch9;
    if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch10;
    Channel_Vector face_velocity_backup_channels;
    face_velocity_backup_channels(0)                    = &Struct_type::ch11;
    face_velocity_backup_channels(1)                    = &Struct_type::ch12;
    if(d==3) face_velocity_backup_channels(2)           = &Struct_type::ch13;
    T Struct_type::* temp_channel                       = &Struct_type::ch14;
    // backup face velocity
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis) {
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity_channels(axis),
                                                face_velocity_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));}  
        Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Velocities(*hierarchy,face_velocity_channels,face_velocity_backup_channels,interpolated_face_velocity_channels,temp_channel,dt);
}
//######################################################################
// Set_Neumann_Faces_Inside_Sources
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Set_Neumann_Faces_Inside_Sources()
{
    for(int level=0;level<levels;++level){
        auto flags=hierarchy->Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
        for(Grid_Iterator_Face<T,d> iterator(hierarchy->Lattice(level));iterator.Valid();iterator.Next()){
            const int axis=iterator.Axis();const T_INDEX& face_index=iterator.Face_Index();
            uint64_t face_offset=Flag_array_mask::Linear_Offset(face_index._data);
            const uint64_t face_active_mask=Topology_Helper::Face_Active_Mask(axis);
            const uint64_t face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
            const TV X=hierarchy->Lattice(level).Face(axis,face_index);
            for(size_t i=0;i<velocity_sources.size();++i) if(velocity_sources(i)->Inside(X)){
                if(hierarchy->template Set<unsigned>(level,&Struct_type::flags).Is_Set(face_offset,face_valid_mask))
                    flags(face_offset)&=~face_active_mask;
                break;}}}
}
//######################################################################
// Initialize_Velocity_Field
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Initialize_Velocity_Field()
{
    Array<TV> source_velocity;
    source_velocity.Append(TV::Axis_Vector(1)*bv);
    for(int level=0;level<levels;++level)
        if(uvf)Uniform_Velocity_Field_Initializer<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,bv,level);
        else Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity_sources,source_velocity,face_velocity_channels,level);
}
//######################################################################
// Project
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Project()
{
    using Multigrid_struct_type             = Multigrid_Data<T>;
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;

    // set up divergence channel
    T Struct_type::* divergence_channel     = &Struct_type::ch8;
    T Struct_type::* pressure_channel       = &Struct_type::ch9;

    // clear
    for(int level=0;level<levels;++level){
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),pressure_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),divergence_channel);}

    Array<TV> source_velocity;
    source_velocity.Append(TV::Axis_Vector(1)*bv);
    // set boundary conditions
    for(int level=0;level<levels;++level){Boundary_Condition_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,pressure_channel,level);
       if(!uvf) Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity_sources,source_velocity,face_velocity_channels,level);}
        
    // compute divergence
    Hierarchy_Projection::Compute_Divergence(*hierarchy,face_velocity_channels,divergence_channel);

    Poisson_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,3,1,200);

    T Struct_type::* q_channel              = &Struct_type::ch10;
    T Struct_type::* r_channel              = &Struct_type::ch11;
    T Struct_type::* s_channel              = &Struct_type::ch12;
    T Struct_type::* k_channel              = &Struct_type::ch12;
    T Struct_type::* z_channel              = &Struct_type::ch13;

    // clear all channels
    for(int level=0;level<levels;++level){
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);}

    CG_Vector<Struct_type,T,d> x_V(*hierarchy,pressure_channel);
    CG_Vector<Struct_type,T,d> b_V(*hierarchy,divergence_channel);
    CG_Vector<Struct_type,T,d> q_V(*hierarchy,q_channel);
    CG_Vector<Struct_type,T,d> r_V(*hierarchy,r_channel);
    CG_Vector<Struct_type,T,d> s_V(*hierarchy,s_channel);
    CG_Vector<Struct_type,T,d> k_V(*hierarchy,k_channel);
    CG_Vector<Struct_type,T,d> z_V(*hierarchy,z_channel);

    Conjugate_Gradient<T> cg;
    cg_system.Multiply(x_V,r_V);
    r_V-=b_V;
    const T b_norm=cg_system.Convergence_Norm(r_V);
    Log::cout<<"Norm: "<<b_norm<<std::endl;
    // cg.print_residuals=true;
    // cg.print_diagnostics=true;
    cg.restart_iterations=cg_restart_iterations;
    const T tolerance=std::max(cg_tolerance*b_norm,(T)1e-6);
    cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);

    for(int level=0;level<levels;++level)
        Apply_Pressure<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,pressure_channel,level);
}
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Register_Options()
{
    Base::Register_Options();

    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL number.");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Integer_Argument("-test_number",1,"Test number.");
    parse_args->Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(64.),"n","Grid resolution");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(64.),"n","Grid resolution");
    parse_args->Add_Double_Argument("-diff_coeff",(T)1e-3,"diffusion coefficient.");
    parse_args->Add_Double_Argument("-fc",(T)0.,"fc.");
    parse_args->Add_Double_Argument("-bv",(T)1.,"Background velocity(along y axis).");
    parse_args->Add_Double_Argument("-sr",(T)1.,"Source rate");
    parse_args->Add_Double_Argument("-tau",(T)1.,"tau.");
    parse_args->Add_Option_Argument("-ficks","Fick's diffusion.");
    parse_args->Add_Option_Argument("-nd","Turn off diffusion");
    parse_args->Add_Option_Argument("-ed","Explicit diffusion");
    parse_args->Add_Option_Argument("-uvf","Uniform velocity field");
    // for CG
    parse_args->Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args->Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args->Add_Double_Argument("-cg_tolerance",1e-4,"CG tolerance");
}
//######################################################################
// Parse_Options
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Parse_Options()
{
    Base::Parse_Options();

    cfl=(T)parse_args->Get_Double_Value("-cfl");
    levels=parse_args->Get_Integer_Value("-levels");
    test_number=parse_args->Get_Integer_Value("-test_number");
    mg_levels=parse_args->Get_Integer_Value("-mg_levels");
    number_of_threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(number_of_threads);
    if(d==2){auto cell_counts_2d=parse_args->Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_2d(v);}
    else{auto cell_counts_3d=parse_args->Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_3d(v);}
    diff_coeff=parse_args->Get_Double_Value("-diff_coeff");
    bv=parse_args->Get_Double_Value("-bv");
    source_rate=parse_args->Get_Double_Value("-sr");
    Fc=parse_args->Get_Double_Value("-fc");
    tau=parse_args->Get_Double_Value("-tau");
    FICKS=parse_args->Get_Option_Value("-ficks");
    uvf=parse_args->Get_Option_Value("-uvf");
    nd=parse_args->Get_Option_Value("-nd");
    explicit_diffusion=parse_args->Get_Option_Value("-ed");
    cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    cg_restart_iterations=parse_args->Get_Integer_Value("-cg_restart_iterations");
    cg_tolerance=(T)parse_args->Get_Double_Value("-cg_tolerance");
    if(nd){diff_coeff=(T)0.;tau=(T)1.;Fc=(T)0.;}
    switch (test_number){
    case 1:{const_density_source=false;const_density_value=(T)0.;uvf=false;}break;
    case 2:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;
    case 3:{const_density_source=false;const_density_value=(T)0.;uvf=false;}break;
    case 4:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;
    case 5:{const_density_source=false;const_density_value=(T)0.;uvf=true;}break;
    case 6:{const_density_source=true;const_density_value=(T)1.;uvf=true;}break;
    case 7:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;
    case 8:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;}
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Write_Output_Files(const int frame) const
{
    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/frame_title",frame_title);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_density",density_channel);
    // hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_u",face_velocity_channels(0));
    // hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_v",face_velocity_channels(1));
    // if(d==3) hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_w",face_velocity_channels(2));
    // Write_To_File_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,output_directory+"/density_data/"+std::to_string(frame)+".txt");
}
//######################################################################
// Read_Output_Files
//######################################################################
template<class T,int d> void Smoke_Example<T,d>::
Read_Output_Files(const int frame)
{
}
//######################################################################
template class Nova::Smoke_Example<float,2>;
template class Nova::Smoke_Example<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Smoke_Example<double,2>;
template class Nova::Smoke_Example<double,3>;
#endif

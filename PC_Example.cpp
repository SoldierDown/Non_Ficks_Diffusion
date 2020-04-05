//!#####################################################################
//! \file PC_Example.cpp
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
#include "Source_Adder.h"
#include "Initialize_Dirichlet_Cells.h"
#include "Poisson_Solver/Multigrid_Data.h"
#include "PC_Example.h"
#include "Clamp_Helper.h"
#include "Lap_Calculator.h"
#include "Density_Backup_Helper.h"
#include "Face_Qc_Updater.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Advection_Helper.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Averaging_Helper.h"
#include "Density_Clamp_Helper.h"
#include "Flip_Helper.h"
#include "Write_To_File_Helper.h"
#include "M_Calculator.h"
#include "Face_Vector_Copy.h"
#include "Density_Plus_Poly_Term_Helper.h"
#include "Axis_Finite_Differential_Helper.h"
#include "Psi_Evaluation_Helper.h"
#include "Epsilon_Epsilon_Prime_dSdX_Helper.h"
#include "Epsilon_Squared_Grad_S_Helper.h"
#include "Add_Constant.h"
#include "K_Gradient_T_Helper.h"
#include "Add_Random_Term.h"
#include <omp.h>
using namespace Nova;
namespace Nova{
extern int number_of_threads;
}
//######################################################################
// Constructor
//######################################################################
template<class T,int d> PC_Example<T,d>::
PC_Example()
    :Base(),hierarchy(nullptr),rasterizer(nullptr)
{
    explicit_diffusion=true;
    // for updating density
    tau_p=(T)3.e-4;
    // for updating face_qs_channels
    tau_1=(T)0.; SR=(T)0.; K=(T)1.2;
    // for updating face_qt_channels
    tau_2=(T)0.; k=(T)1.;
    // for updating m_channel
    m_alpha=(T).9; gamma=(T)10.; Teq=(T)1.;
    // for updating epsilon_channel
    omega=6; delta=(T).05;
    bv=(T)1.;

    density_channel                     = &Struct_type::ch0;            // intermedia 
    density_backup_channel              = &Struct_type::ch1;            // S^n
    T_channel                           = &Struct_type::ch2;
    T_backup_channel                    = &Struct_type::ch3;

    face_qs_channels(0)                 = &Struct_type::ch4;
    face_qs_channels(1)                 = &Struct_type::ch5;
    if(d==3) face_qs_channels(2)        = &Struct_type::ch6;
    face_qs_backup_channels(0)          = &Struct_type::ch7;
    face_qs_backup_channels(1)          = &Struct_type::ch8;
    if(d==3) face_qs_backup_channels(2) = &Struct_type::ch9;
    
    face_qt_channels(0)                 = &Struct_type::ch10;
    face_qt_channels(1)                 = &Struct_type::ch11;
    if(d==3) face_qt_channels(2)        = &Struct_type::ch12;
    face_qt_backup_channels(0)          = &Struct_type::ch13;
    face_qt_backup_channels(1)          = &Struct_type::ch14;
    if(d==3) face_qt_backup_channels(2) = &Struct_type::ch15;
    
    face_velocity_channels(0)           = &Struct_type::ch16;
    face_velocity_channels(1)           = &Struct_type::ch17;
    if(d==3) face_velocity_channels(2)  = &Struct_type::ch18;
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Initialize()
{
    Initialize_SPGrid();
    Initialize_State();
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Initialize_SPGrid()
{
    Log::Scope scope("Initialize_SPGrid");
    Initialize_Rasterizer();
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
template<class T,int d> void PC_Example<T,d>::
Limit_Dt(T& dt,const T time)
{
    T dt_convection=(T)0.;

    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);

    for(int level=0;level<levels;++level)
        Compute_Time_Step<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,
                                           other_face_offsets,level,dt_convection);
    Log::cout<<"dt_convection: "<<dt_convection<<std::endl;
    if(dt_convection>(T)1e-5) dt=cfl/dt_convection;
    Log::cout<<"Time Step: "<<dt<<std::endl;
}
//######################################################################
// Advect_Scalar_Field
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Advect_Scalar_Field(const T dt)
{
    Advect_Density(dt);
    Advect_Temperature(dt);
}
//######################################################################
// Advect_Temperature
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Advect_Density(const T dt)
{
    Log::Scope scope("Advect Density");
    Channel_Vector cell_velocity_channels;
    cell_velocity_channels(0)               = &Struct_type::ch20;
    cell_velocity_channels(1)               = &Struct_type::ch21;
    if(d==3) cell_velocity_channels(2)      = &Struct_type::ch22;
    T Struct_type::* temp_channel           = &Struct_type::ch23;
    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
    Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Cells(*hierarchy,hierarchy->Allocator(0),hierarchy->Blocks(0),face_velocity_channels,cell_velocity_channels,other_face_offsets);
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Density(*hierarchy,cell_velocity_channels,density_channel,temp_channel,dt);
}
//######################################################################
// Advect_Density
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Advect_Temperature(const T dt)
{
    Log::Scope scope("Advect Temperature");
    Channel_Vector cell_velocity_channels;
    cell_velocity_channels(0)               = &Struct_type::ch20;
    cell_velocity_channels(1)               = &Struct_type::ch21;
    if(d==3) cell_velocity_channels(2)      = &Struct_type::ch22;
    T Struct_type::* temp_channel           = &Struct_type::ch23;
    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
    Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Cells(*hierarchy,hierarchy->Allocator(0),hierarchy->Blocks(0),face_velocity_channels,cell_velocity_channels,other_face_offsets);
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Density(*hierarchy,cell_velocity_channels,T_channel,temp_channel,dt);
}
//######################################################################
// Advect_Face_Vector_Field
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Advect_Face_Vector_Field(const T dt)
{
    Advect_Face_Qs(dt);
    Advect_Face_Qt(dt);
}
//######################################################################
// Advect_Face_Qs
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Advect_Face_Qs(const T dt)
{
    Log::Scope scope("Advect Face Qs");
    // Advect face vector by axis
    Channel_Vector interpolated_face_velocity_channels;
    interpolated_face_velocity_channels(0)              = &Struct_type::ch20;
    interpolated_face_velocity_channels(1)              = &Struct_type::ch21;
    if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch22; 
    T Struct_type::* temp_channel                       = &Struct_type::ch23;  
    // Clear
    for(int level=0;level<levels;++level) {for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),interpolated_face_velocity_channels(v));
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),temp_channel);} 
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Vector(*hierarchy,face_qs_channels,face_velocity_channels,interpolated_face_velocity_channels,temp_channel,dt);
}
//######################################################################
// Advect_Face_Qt
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Advect_Face_Qt(const T dt)
{
    Log::Scope scope("Advect Face Qt");
    // Advect face vector by axis
    Channel_Vector interpolated_face_velocity_channels;
    interpolated_face_velocity_channels(0)              = &Struct_type::ch20;
    interpolated_face_velocity_channels(1)              = &Struct_type::ch21;
    if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch22; 
    T Struct_type::* temp_channel                       = &Struct_type::ch23;  
    // Clear
    for(int level=0;level<levels;++level) {for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),interpolated_face_velocity_channels(v));
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),temp_channel);} 
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Vector(*hierarchy,face_qt_channels,face_velocity_channels,interpolated_face_velocity_channels,temp_channel,dt);
}
//######################################################################
// Update_Density
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Update_Density(const T dt)
{
    Log::Scope scope("Update Density");
    if(explicit_diffusion) Explicitly_Update_Density(dt);
    else Implicitly_Update_Density(dt);    
}
//######################################################################
// Explicitly_Update_Density
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Explicitly_Update_Density(const T dt)
{
    Add_Divergence_Term_To_Density(dt);
    Add_Differential_Term_To_Density(dt);
    Add_Poly_Term_To_Density(dt);
    Add_Random_Term_To_Density(dt);
}
//######################################################################
// Add_Divergence_Term_To_Density
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Divergence_Term_To_Density(const T dt)
{
    Log::Scope scope("Add_Divergence_Term_To_Density");
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;
    T Struct_type::* div_qs_channel         = &Struct_type::ch20;
    for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_qs_channel);
    // compute -div(Qc) at time n
    Hierarchy_Projection::Compute_Divergence(*hierarchy,face_qs_backup_channels,div_qs_channel); 
    const T dt_over_tau_p=(T)dt/tau_p;
    for(int level=0;level<levels;++level) SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dt_over_tau_p,div_qs_channel,density_channel,density_channel,Cell_Type_Interior);
}
//######################################################################
// Add_Differential_Term_To_Density
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Differential_Term_To_Density(const T dt)
{
    Log::Scope scope("Add_Differential_Term_To_Density");
    epsilon_channel                 = &Struct_type::ch19;
    Channel_Vector dSdX_channels;
    dSdX_channels(0)                = &Struct_type::ch20;
    dSdX_channels(1)                = &Struct_type::ch21;
    if(d==3) dSdX_channels(2)       = &Struct_type::ch22;

    T Struct_type::* psi_channel    = &Struct_type::ch23;
    for(int level=0;level<levels;++level){ for(int axis=0;axis<d;++axis) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dSdX_channels(axis));
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),psi_channel);
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),epsilon_channel);}
    // compute dS/dx and dS/dy (dS/dz)
    const T one_over_2dx=(T).5*hierarchy->Lattice(0).one_over_dX(0);
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis)
        Axis_Finite_Differential_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel,dSdX_channels(axis),one_over_2dx,axis);
    // evaluate psi
    for(int level=0;level<levels;++level) 
        Psi_Evaluation_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dSdX_channels,psi_channel);
    // eps * eps' 
    for(int level=0;level<levels;++level) 
        Epsilon_Epsilon_Prime_dSdX_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dSdX_channels,psi_channel,epsilon_channel,
                                                            omega,delta);
    const T dt_over_tau_p=dt/tau_p;
    T Struct_type::* temp_channel = &Struct_type::ch24;    
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis){
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),temp_channel);
        Axis_Finite_Differential_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dSdX_channels((axis+1)%d),temp_channel,one_over_2dx,axis);
        SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dt_over_tau_p,temp_channel,density_channel,density_channel,Cell_Type_Interior);}

}
//######################################################################
// Add_Poly_Term_To_Density
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Poly_Term_To_Density(const T dt)
{
    Log::Scope scope("Add_Poly_Term_To_Density");
    T Struct_type::* m_channel                     = &Struct_type::ch20;
    for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),m_channel);
    for(int level=0;level<levels;++level) M_Calculator<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),m_channel,T_backup_channel);
    const T dt_over_tau_p=dt/tau_p;
    for(int level=0;level<levels;++level) Density_Plus_Poly_Term_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,density_backup_channel,m_channel,dt_over_tau_p);
}
//######################################################################
// Add_Random_Term_To_Density
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Random_Term_To_Density(const T dt)
{
    const T alpha=(T).01;
    for(int level=0;level<levels;++level) Add_Random_Term<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,density_backup_channel,alpha);
}
//######################################################################
// Implicitly_Update_Density
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Implicitly_Update_Density(const T dt)
{
}
//######################################################################
// Update_Face_Qs
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Update_Face_Qs(const T dt)
{
    if(explicit_diffusion) Explicitly_Update_Face_Qs(dt);
    else Implicitly_Update_Face_Qs(dt);    
}
//######################################################################
// Explicitly_Update_Face_Qs
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Explicitly_Update_Face_Qs(const T dt)
{
    Add_Linear_Term_To_Face_Qs(dt);
    Add_Gradient_Term_To_Face_Qs(dt);
}
//######################################################################
// Add_Linear_Term_To_Face_Qs
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Linear_Term_To_Face_Qs(const T dt)
{
    const T minus_dt_over_tau_1=-dt/tau_1;
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis)  
        SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),
            minus_dt_over_tau_1,face_qs_backup_channels(axis),face_qs_channels(axis),face_qs_channels(axis),Topology_Helper::Face_Active_Mask(axis));
}
//######################################################################
// Add_Gradient_Term_To_Face_Qs
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Gradient_Term_To_Face_Qs(const T dt)
{
    const T minus_dt_over_tau_1=-dt/tau_1; const T one_over_dx=hierarchy->Lattice(0).one_over_dX(0);
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis)  
    Epsilon_Squared_Grad_S_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),
            face_qs_channels,epsilon_channel,density_backup_channel,minus_dt_over_tau_1,one_over_dx,level);
}
//######################################################################
// Implicitly_Update_Face_Qs
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Implicitly_Update_Face_Qs(const T dt)
{
}
//######################################################################
// Update_Temperature
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Update_Temperature(const T dt)
{
    if(explicit_diffusion) Explicitly_Update_Temperature(dt);
    else Implicitly_Update_Temperature(dt);    
}
//######################################################################
// Explicitly_Update_Temperature
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Explicitly_Update_Temperature(const T dt)
{
    Add_Constant_Term_To_Temperature(dt);
    Add_dSdt_Term_To_Temperature(dt);
    Add_Divergence_Term_To_Temperature(dt);
}
//######################################################################
// Add_Constant_Term_To_Temperature
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Add_Constant_Term_To_Temperature(const T dt)
{
    for(int level=0;level<levels;++level) Add_Constant<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),T_channel,SR*dt);
}
//######################################################################
// Add_dSdt_Term_To_Temperature
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Add_dSdt_Term_To_Temperature(const T dt)
{
    T Struct_type::* temp_channel         = &Struct_type::ch20;
    for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),temp_channel);
    // compute S^(n+1) - S^(n)
    for(int level=0;level<levels;++level) SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,density_backup_channel,temp_channel,Cell_Type_Interior);
    // add K*dS/dt to T
    for(int level=0;level<levels;++level) SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),K,temp_channel,T_channel,T_channel,Cell_Type_Interior);
}
//######################################################################
// Add_Divergence_Term_To_Temperature
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Add_Divergence_Term_To_Temperature(const T dt)
{
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;
    T Struct_type::* div_qt_channel         = &Struct_type::ch20;
    for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_qt_channel);
    // compute -div(Qc) at time n
    Hierarchy_Projection::Compute_Divergence(*hierarchy,face_qs_channels,div_qt_channel); 
    for(int level=0;level<levels;++level) SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dt,div_qt_channel,T_channel,T_channel,Cell_Type_Interior);
}
//######################################################################
// Implicitly_Update_Temperature
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Implicitly_Update_Temperature(const T dt)
{
}
//######################################################################
// Update_Face_Qt
//######################################################################
template    <class T,int d> void PC_Example<T,d>::
Update_Face_Qt(const T dt)
{
    if(explicit_diffusion) Explicitly_Update_Face_Qt(dt);
    else Implicitly_Update_Face_Qt(dt);    
}
//######################################################################
// Explicitly_Update_Face_Qt
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Explicitly_Update_Face_Qt(const T dt)
{
    Add_Linear_Term_To_Face_Qt(dt);
    Add_Gradient_Term_To_Face_Qt(dt);
}
//######################################################################
// Add_Linear_Term_To_Face_Qt
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Linear_Term_To_Face_Qt(const T dt)
{
    const T minus_dt_over_tau_2=-dt/tau_2;
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis)  
        SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),
            minus_dt_over_tau_2,face_qt_backup_channels(axis),face_qt_channels(axis),face_qt_channels(axis),Topology_Helper::Face_Active_Mask(axis));
}
//######################################################################
// Add_Gradient_Term_To_Face_Qt
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Gradient_Term_To_Face_Qt(const T dt)
{
    const T minus_dt_over_tau_2=-dt/tau_1; const T one_over_dx=hierarchy->Lattice(0).one_over_dX(0);
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis)  
    K_Gradient_T_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),
            face_qt_channels,T_backup_channel,k,minus_dt_over_tau_2,one_over_dx,level);
}
//######################################################################
// Implicitly_Update_Face_Qt
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Implicitly_Update_Face_Qt(const T dt)
{
}
//######################################################################
// Backup
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Backup()
{
    Backup_Density(); Backup_Temperature(); Backup_Qs(); Backup_Qt();
}
//######################################################################
// Backup_Density
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Backup_Density()
{
    for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel);
    for(int level=0;level<levels;++level) SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,
                                                                                density_backup_channel,(unsigned)Cell_Type_Interior);
}
//######################################################################
// Backup_Temperature
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Backup_Temperature()
{
    for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),T_backup_channel);
    for(int level=0;level<levels;++level) SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),T_channel,
                                                                                T_backup_channel,(unsigned)Cell_Type_Interior);
}
//######################################################################
// Backup_Qs
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Backup_Qs()
{
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qs_backup_channels(axis));
    for(int level=0;level<levels;++level) Face_Vector_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qs_channels,face_qs_backup_channels);
}
//######################################################################
// Backup_Qt
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Backup_Qt()
{
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qt_backup_channels(axis));
    for(int level=0;level<levels;++level) Face_Vector_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_qt_channels,face_qt_backup_channels);
}
//######################################################################
// Modify_Density_With_Sources
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Modify_Density_With_Sources()
{
    for(int level=0;level<levels;++level)
        Density_Modifier<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,sources,level);
}
//######################################################################
// Add_Source
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Add_Source(const T dt)
{
    // for(int level=0;level<levels;++level)
    //     Source_Adder<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,sources,source_rate,dt,level);
}
//######################################################################
// Advect_Face_Velocities
//######################################################################
// template<class T,int d> void PC_Example<T,d>::
// Advect_Face_Velocities(const T dt)
// {
//     Channel_Vector interpolated_face_velocity_channels;
//     interpolated_face_velocity_channels(0)              = &Struct_type::ch8;
//     interpolated_face_velocity_channels(1)              = &Struct_type::ch9;
//     if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch10;
//     T Struct_type::* temp_channel                       = &Struct_type::ch11;
//     Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Velocities(*hierarchy,face_velocity_channels,interpolated_face_velocity_channels,temp_channel,dt);
// }
//######################################################################
// Set_Neumann_Faces_Inside_Sources
//######################################################################
template<class T,int d> void PC_Example<T,d>::
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
            for(size_t i=0;i<sources.size();++i) if(sources(i)->Inside(X)){
                if(hierarchy->template Set<unsigned>(level,&Struct_type::flags).Is_Set(face_offset,face_valid_mask))
                    flags(face_offset)&=~face_active_mask;
                break;}}}
}
//######################################################################
// Initialize_Velocity_Field
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Initialize_Velocity_Field()
{
    // for(int level=0;level<levels;++level)
    //     Uniform_Velocity_Field_Initializer<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,bv,level);
}
//######################################################################
// Project
//######################################################################
// template<class T,int d> void PC_Example<T,d>::
// Project(const T dt)
// {
//     using Multigrid_struct_type             = Multigrid_Data<T>;
//     using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;

//     // set up divergence channel
//     T Struct_type::* divergence_channel     = &Struct_type::ch8;

//     // clear
//     for(int level=0;level<levels;++level){
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),pressure_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),divergence_channel);}

//     // set boundary conditions
//     for(int level=0;level<levels;++level)
//         Boundary_Condition_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,pressure_channel,bv,level);

//     // compute divergence
//     Hierarchy_Projection::Compute_Divergence(*hierarchy,face_velocity_channels,divergence_channel);

//     Poisson_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,3,1,200);

//     T Struct_type::* q_channel              = &Struct_type::ch9;
//     T Struct_type::* r_channel              = &Struct_type::ch10;
//     T Struct_type::* s_channel              = &Struct_type::ch11;
//     T Struct_type::* k_channel              = &Struct_type::ch11;
//     T Struct_type::* z_channel              = &Struct_type::ch12;

//     // clear all channels
//     for(int level=0;level<levels;++level){
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);}

//     CG_Vector<Struct_type,T,d> x_V(*hierarchy,pressure_channel);
//     CG_Vector<Struct_type,T,d> b_V(*hierarchy,divergence_channel);
//     CG_Vector<Struct_type,T,d> q_V(*hierarchy,q_channel);
//     CG_Vector<Struct_type,T,d> r_V(*hierarchy,r_channel);
//     CG_Vector<Struct_type,T,d> s_V(*hierarchy,s_channel);
//     CG_Vector<Struct_type,T,d> k_V(*hierarchy,k_channel);
//     CG_Vector<Struct_type,T,d> z_V(*hierarchy,z_channel);

//     Conjugate_Gradient<T> cg;
//     cg_system.Multiply(x_V,r_V);
//     r_V-=b_V;
//     const T b_norm=cg_system.Convergence_Norm(r_V);
//     Log::cout<<"Norm: "<<b_norm<<std::endl;
//     // cg.print_residuals=true;
//     // cg.print_diagnostics=true;
//     cg.restart_iterations=cg_restart_iterations;
//     const T tolerance=std::max(cg_tolerance*b_norm,(T)1e-6);
//     cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);

//     for(int level=0;level<levels;++level)
//         Apply_Pressure<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,pressure_channel,level);
// }
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Register_Options()
{
    Base::Register_Options();

    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL number.");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(64.),"n","Grid resolution");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(64.),"n","Grid resolution");
    parse_args->Add_Option_Argument("-ed","Explicit diffusion");
    parse_args->Add_Integer_Argument("-omega",4,"Number of branches.");
    // for CG
    parse_args->Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args->Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args->Add_Double_Argument("-cg_tolerance",1e-4,"CG tolerance");
}
//######################################################################
// Parse_Options
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Parse_Options()
{
    Base::Parse_Options();
    cfl=(T)parse_args->Get_Double_Value("-cfl");
    levels=parse_args->Get_Integer_Value("-levels");
    mg_levels=parse_args->Get_Integer_Value("-mg_levels");
    number_of_threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(number_of_threads);
    if(d==2){auto cell_counts_2d=parse_args->Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_2d(v);}
    else{auto cell_counts_3d=parse_args->Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_3d(v);}
    explicit_diffusion=parse_args->Get_Option_Value("-ed");
    cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    omega=parse_args->Get_Integer_Value("-omega");
    cg_restart_iterations=parse_args->Get_Integer_Value("-cg_restart_iterations");
    cg_tolerance=(T)parse_args->Get_Double_Value("-cg_tolerance");
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Write_Output_Files(const int frame) const
{
    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/frame_title",frame_title);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_density",density_channel);
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_u",face_velocity_channels(0));
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_v",face_velocity_channels(1));
    if(d==3) hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_w",face_velocity_channels(2));
    Write_To_File_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,output_directory+"/density_data/"+std::to_string(frame)+".txt");
}
//######################################################################
// Read_Output_Files
//######################################################################
template<class T,int d> void PC_Example<T,d>::
Read_Output_Files(const int frame)
{
}
//######################################################################
template class Nova::PC_Example<float,2>;
template class Nova::PC_Example<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::PC_Example<double,2>;
template class Nova::PC_Example<double,3>;
#endif

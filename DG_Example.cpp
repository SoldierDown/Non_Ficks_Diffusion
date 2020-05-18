//!#####################################################################
//! \file DG_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Advection/Grid_Hierarchy_Advection.h>
#include <nova/Tools/Grids/Grid_Iterator_Face.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include "Ficks_Solver/Ficks_CG_System.h"
#include "Ficks_RHS_Helper.h"
#include "DG_Ficks_T_RHS_Helper.h"
#include "Non_Ficks_Solver/Non_Ficks_CG_System.h"
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
#include "Initialize_Neumann_Cells.h"
#include "Poisson_Solver/Multigrid_Data.h"
#include "Range_Checker.h"
#include "DG_Example.h"
#include "SF_M_Calculator.h"
#include "Clamp_Helper.h"
#include "Lap_Calculator.h"
#include "Density_Backup_Helper.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Advection_Helper.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Averaging_Helper.h"
#include "Density_Clamp_Helper.h"
#include "Flip_Helper.h"
#include "Gradient_Calculator.h"
#include "Sigma_Calculator.h"
#include "Write_To_File_Helper.h"
#include "DG_M_Calculator.h"
#include "DG_Epsilon_Squared_Grad_S_Helper.h"
#include "Face_Vector_Copy.h"
#include "Density_Plus_Poly_Term_Helper.h"
#include "Axis_Finite_Differential_Helper.h"
#include "Psi_Evaluation_Helper.h"
#include "Add_Constant.h"
#include "Traverse_Helper.h"
#include "K_Gradient_T_Helper.h"
#include "Add_Random_Term.h"
#include "Boundary_Check.h"
#include "Masked_Multiply.h"
#include <omp.h>
using namespace Nova;
namespace Nova{
extern int number_of_threads;
}
//######################################################################
// Constructor
//######################################################################
template<class T,int d> DG_Example<T,d>::
DG_Example()
    :Base(),hierarchy(nullptr),rasterizer(nullptr)
{
    random.Set_Seed(0);
    // // for updating density
    tau_s=(T)3.e-4;
    // for updating face_qsc_channels
    SR=(T)0.;
    // for updating face_qtc_channels
    // for updating m_channel
    // for updating epsilon_channel

    density_channel                     = &Struct_type::ch0;            // intermedia 
    density_backup_channel              = &Struct_type::ch1;            // S^n
    T_channel                           = &Struct_type::ch2;
    T_backup_channel                    = &Struct_type::ch3;

    face_qsc_channels(0)                 = &Struct_type::ch4;
    face_qsc_channels(1)                 = &Struct_type::ch5;
    if(d==3) face_qsc_channels(2)        = &Struct_type::ch6;
    face_qsc_backup_channels(0)          = &Struct_type::ch7;
    face_qsc_backup_channels(1)          = &Struct_type::ch8;
    if(d==3) face_qsc_backup_channels(2) = &Struct_type::ch9;
    
    face_qtc_channels(0)                 = &Struct_type::ch10;
    face_qtc_channels(1)                 = &Struct_type::ch11;
    if(d==3) face_qtc_channels(2)        = &Struct_type::ch12;
    face_qtc_backup_channels(0)          = &Struct_type::ch13;
    face_qtc_backup_channels(1)          = &Struct_type::ch14;
    if(d==3) face_qtc_backup_channels(2) = &Struct_type::ch15;
    
    face_velocity_channels(0)           = &Struct_type::ch16;
    face_velocity_channels(1)           = &Struct_type::ch17;
    if(d==3) face_velocity_channels(2)  = &Struct_type::ch18;
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Initialize()
{
    Initialize_SPGrid();
    Initialize_State();
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T,int d> void DG_Example<T,d>::
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
    // Initialize_Dirichlet_Cells<Struct_type,T,d>(*hierarchy,domain_walls);
    Initialize_Neumann_Cells<Struct_type,T,d>(*hierarchy,domain_walls);
    //Set_Neumann_Faces_Inside_Sources();
    hierarchy->Update_Block_Offsets();
    hierarchy->Initialize_Red_Black_Partition(2*number_of_threads);
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Limit_Dt(T& dt,const T time)
{
    T dt_convection=(T)0.;

    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);

    Compute_Time_Step<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(0),face_velocity_channels,
                                           other_face_offsets,0,dt_convection);
    Log::cout<<"dt_convection: "<<dt_convection<<std::endl;
    if(dt_convection>(T)1e-5) dt=cfl/dt_convection;
    Log::cout<<"Time Step: "<<dt<<std::endl;
}
//######################################################################
// Advect_Scalar_Field
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Advect_Scalar_Field(const T dt)
{
    Advect_Density(dt);
    Advect_Temperature(dt);
}
//######################################################################
// Advect_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Advect_Density(const T dt)
{
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
template    <class T,int d> void DG_Example<T,d>::
Advect_Temperature(const T dt)
{
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
template    <class T,int d> void DG_Example<T,d>::
Advect_Face_Vector_Field(const T dt)
{
    Advect_Face_Qsc(dt);
    Advect_Face_Qtc(dt);
}
//######################################################################
// Advect_Face_Qsc
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Advect_Face_Qsc(const T dt)
{
    // Advect face vector by axis
    Channel_Vector interpolated_face_velocity_channels;
    interpolated_face_velocity_channels(0)              = &Struct_type::ch20;
    interpolated_face_velocity_channels(1)              = &Struct_type::ch21;
    if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch22; 
    T Struct_type::* temp_channel                       = &Struct_type::ch23;  
    // Clear
    for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),interpolated_face_velocity_channels(v));
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),temp_channel); 
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Vector(*hierarchy,face_qsc_channels,face_velocity_channels,interpolated_face_velocity_channels,temp_channel,dt);
}
//######################################################################
// Advect_Face_Qtc
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Advect_Face_Qtc(const T dt)
{
    // Advect face vector by axis
    Channel_Vector interpolated_face_velocity_channels;
    interpolated_face_velocity_channels(0)              = &Struct_type::ch20;
    interpolated_face_velocity_channels(1)              = &Struct_type::ch21;
    if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch22; 
    T Struct_type::* temp_channel                       = &Struct_type::ch23;  
    // Clear
    for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),interpolated_face_velocity_channels(v));
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),temp_channel);
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Vector(*hierarchy,face_qtc_channels,face_velocity_channels,interpolated_face_velocity_channels,temp_channel,dt);
}
//######################################################################
// Update_Density
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Update_Density(const T dt)
{
    if(explicit_diffusion) Explicitly_Update_Density(dt);
    else Implicitly_Update_Density(dt);    
}
//######################################################################
// Explicitly_Update_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Explicitly_Update_Density(const T dt)
{
    Add_Poly_Term_To_Density(dt);
    Add_Random_Term_To_Density(dt);
    Add_Laplacian_Term_To_Density(dt);
    if(!FICKS) Add_Divergence_Term_To_Density(dt);
}
//######################################################################
// Add_Laplacian_Term_To_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Laplacian_Term_To_Density(const T dt)
{
    using Hierarchy_Projection                          = Grid_Hierarchy_Projection<Struct_type,T,d>;
    T Struct_type::* laplacian_channel                  = &Struct_type::ch19;
    if(FICKS) fc_1=(T)1.;
    const T one_over_dx2=Nova_Utilities::Sqr(hierarchy->Lattice(0).one_over_dX(0)); const T dt_over_tau_s=dt/tau_s;
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),laplacian_channel);
    Lap_Calculator<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_backup_channel,laplacian_channel,one_over_dx2);
    SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),fc_1*Nova_Utilities::Sqr(epsilon)*dt_over_tau_s,laplacian_channel,
                                            density_channel,density_channel,Cell_Type_Interior);
}
//######################################################################
// Add_Divergence_Term_To_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Divergence_Term_To_Density(const T dt)
{
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;
    T Struct_type::* div_qs_channel         = &Struct_type::ch19;
    const T dt_over_tau_s=(T)dt/tau_s;
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),div_qs_channel);
    // compute -div(Qc) at time n
    Hierarchy_Projection::Compute_Divergence(*hierarchy,face_qsc_backup_channels,div_qs_channel); 
    SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),dt_over_tau_s,div_qs_channel,
                                            density_channel,density_channel,Cell_Type_Interior);
}
//######################################################################
// Add_Poly_Term_To_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Poly_Term_To_Density(const T dt)
{
    // calculate v
    Channel_Vector dSdX_channels;
    dSdX_channels(0)                                = &Struct_type::ch19;
    dSdX_channels(1)                                = &Struct_type::ch20;
    if(d==3) dSdX_channels(2)                       = &Struct_type::ch21;
    T Struct_type::* sigma_channel                  = &Struct_type::ch22;
    // clear
    for(int axis=0;axis<d;++axis) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),dSdX_channels(axis));
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),sigma_channel);
    
    const TV one_over_2dx=(T).5*hierarchy->Lattice(0).one_over_dX;
    // compute grad(S) at cell center
    for(int axis=0;axis<d;++axis) Axis_Finite_Differential_Helper<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),
                                    density_backup_channel,dSdX_channels(axis),one_over_2dx(axis),axis);
    // compute sigma
    Sigma_Calculator<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),dSdX_channels,sigma_channel,delta);
    
    T Struct_type::* m_channel                      = &Struct_type::ch23;
    const T dt_over_tau_s=dt/tau_s;
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),m_channel);
    DG_M_Calculator<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),m_channel,sigma_channel,T_backup_channel,m_alpha,K*gamma);
    Range_Checker<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),m_channel,(T)-.5,(T).5);
    Density_Plus_Poly_Term_Helper<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,density_backup_channel,m_channel,dt_over_tau_s);    
}
//######################################################################
// Add_Random_Term_To_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Random_Term_To_Density(const T dt)
{
    const T random_value=random.Get_Uniform_Number(0,1);
    Add_Random_Term<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,density_backup_channel,random_factor*dt,random_value);
}
//######################################################################
// Implicitly_Update_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Implicitly_Update_Density(const T dt)
{
    Log::cout<<"Updating density: Fick's"<<std::endl;
    using Multigrid_struct_type 	                = Multigrid_Data<T>;
    Add_Poly_Term_To_Density(dt);
    Add_Random_Term_To_Density(dt);
    // so far S^n->S*
	const Grid<T,d>& grid=hierarchy->Lattice(0);
    const T eps2_dt_over_taus=Nova_Utilities::Sqr(epsilon)*dt/tau_s;
	Ficks_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,eps2_dt_over_taus,3,1,200);
    T Struct_type::* q_channel              = &Struct_type::ch19;
    T Struct_type::* r_channel              = &Struct_type::ch20;
    T Struct_type::* s_channel              = &Struct_type::ch21;
    T Struct_type::* k_channel              = &Struct_type::ch22;
    T Struct_type::* z_channel              = &Struct_type::ch22;
    T Struct_type::* b_channel              = &Struct_type::ch24;
    
    // clear all channels
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),q_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),r_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),s_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),z_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),b_channel);

    const T eps2_dt_over_tau_s_dx2=Nova_Utilities::Sqr(epsilon*grid.one_over_dX(0))*dt/tau_s; 
	Ficks_RHS_Helper<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,b_channel,eps2_dt_over_tau_s_dx2);
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
    const T tolerance=std::max((T)1e-4*b_norm,(T)1e-4);
    cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);
    Clamp_Heler<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel);
}
//######################################################################
// Update_Face_Qsc
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Update_Face_Qsc(const T dt)
{
    if(explicit_diffusion) Explicitly_Update_Face_Qsc(dt);
    else Implicitly_Update_Face_Qsc(dt);    
}
//######################################################################
// Explicitly_Update_Face_Qsc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Explicitly_Update_Face_Qsc(const T dt)
{
    Add_Linear_Term_To_Face_Qsc(dt);
    Add_Gradient_Term_To_Face_Qsc(dt);
}
//######################################################################
// Add_Linear_Term_To_Face_Qsc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Linear_Term_To_Face_Qsc(const T dt)
{
    const T minus_dt_over_tau_1=-dt/tau_1;
    for(int axis=0;axis<d;++axis)  
        SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),minus_dt_over_tau_1,
                                face_qsc_backup_channels(axis),face_qsc_channels(axis),face_qsc_channels(axis),Topology_Helper::Face_Active_Mask(axis));
}
//######################################################################
// Add_Gradient_Term_To_Face_Qsc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Gradient_Term_To_Face_Qsc(const T dt)
{
    const T minus_dt_over_tau_1=-dt/tau_1; const TV one_over_dx=hierarchy->Lattice(0).one_over_dX;
    DG_Epsilon_Squared_Grad_S_Helper<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),
                                    face_qsc_channels,density_backup_channel,epsilon*epsilon,(1-fc_1)*minus_dt_over_tau_1,one_over_dx);
}
//######################################################################
// Implicitly_Update_Face_Qsc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Implicitly_Update_Face_Qsc(const T dt)
{
}
//######################################################################
// Update_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Update_Temperature(const T dt)
{
    if(explicit_diffusion) Explicitly_Update_Temperature(dt);
    else Implicitly_Update_Temperature(dt);    
}
//######################################################################
// Explicitly_Update_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Explicitly_Update_Temperature(const T dt)
{
    Add_Constant_Term_To_Temperature(dt);
    Add_dSdt_Term_To_Temperature(dt);
    Add_Laplacian_Term_To_Temperature(dt);
    if(!FICKS) Add_Divergence_Term_To_Temperature(dt);
}
//######################################################################
// Add_Constant_Term_To_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Add_Constant_Term_To_Temperature(const T dt)
{
    Add_Constant<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),T_channel,SR*dt);
}
//######################################################################
// Add_dSdt_Term_To_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Add_dSdt_Term_To_Temperature(const T dt)
{
    T Struct_type::* dS_channel             = &Struct_type::ch19;
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),dS_channel);
    // compute S^(n+1) - S^(n)
    SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,density_backup_channel,dS_channel,Cell_Type_Interior);
    // add dS to T
    SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),(T)1.,dS_channel,T_channel,T_channel,Cell_Type_Interior);
}
//######################################################################
// Add_Laplacian_Term_To_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Add_Laplacian_Term_To_Temperature(const T dt)
{
    T Struct_type::* lap_T_channel          = &Struct_type::ch19;
    if(FICKS) fc_2=(T)1.;
    const T one_over_dx2=Nova_Utilities::Sqr(hierarchy->Lattice(0).one_over_dX(0)); 
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),lap_T_channel);
    Lap_Calculator<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),T_backup_channel,lap_T_channel,one_over_dx2);
    SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),fc_2*dt,lap_T_channel,T_channel,T_channel,Cell_Type_Interior);
}
//######################################################################
// Add_Divergence_Term_To_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Add_Divergence_Term_To_Temperature(const T dt)
{
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;
    T Struct_type::* div_qt_channel         = &Struct_type::ch19;
    
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),div_qt_channel);
    // compute -div(Qc) at time n
    Hierarchy_Projection::Compute_Divergence(*hierarchy,face_qtc_backup_channels,div_qt_channel); 
    SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),dt,div_qt_channel,T_channel,T_channel,Cell_Type_Interior);
}
//######################################################################
// Implicitly_Update_Temperature
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Implicitly_Update_Temperature(const T dt)
{
    Log::cout<<"Updating temperature: Fick's"<<std::endl;
    using Multigrid_struct_type 	        = Multigrid_Data<T>;
	const Grid<T,d>& grid=hierarchy->Lattice(0);
	Ficks_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,dt,3,1,200);
    T Struct_type::* q_channel              = &Struct_type::ch19;
    T Struct_type::* r_channel              = &Struct_type::ch20;
    T Struct_type::* s_channel              = &Struct_type::ch21;
    T Struct_type::* k_channel              = &Struct_type::ch22;
    T Struct_type::* z_channel              = &Struct_type::ch22;
    T Struct_type::* b_channel              = &Struct_type::ch24;
    
    // clear all channels
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),q_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),r_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),s_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),z_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),b_channel);

    const T dt_over_dx2=dt*Nova_Utilities::Sqr(grid.one_over_dX(0)); 
	DG_Ficks_T_RHS_Helper<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),
                        density_channel,density_backup_channel,T_channel,b_channel,dt_over_dx2,SR*dt);
	CG_Vector<Struct_type,T,d> x_V(*hierarchy,T_channel),b_V(*hierarchy,b_channel),q_V(*hierarchy,q_channel),
                                s_V(*hierarchy,s_channel),r_V(*hierarchy,r_channel),k_V(*hierarchy,k_channel),z_V(*hierarchy,z_channel);   
	Conjugate_Gradient<T> cg;
    cg_system.Multiply(x_V,r_V);
    r_V-=b_V;
    const T b_norm=cg_system.Convergence_Norm(r_V);
    Log::cout<<"Norm: "<<b_norm<<std::endl;
    // cg.print_residuals=true;
    // cg.print_diagnostics=true;
    cg.restart_iterations=cg_restart_iterations;
    const T tolerance=std::max((T)1e-4*b_norm,(T)1e-4);
    cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);
    // Clamp_Heler<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),T_channel);
}
//######################################################################
// Update_Face_Qtc
//######################################################################
template    <class T,int d> void DG_Example<T,d>::
Update_Face_Qtc(const T dt)
{
    if(explicit_diffusion) Explicitly_Update_Face_Qtc(dt);
    else Implicitly_Update_Face_Qtc(dt);    
}
//######################################################################
// Explicitly_Update_Face_Qtc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Explicitly_Update_Face_Qtc(const T dt)
{
    Add_Linear_Term_To_Face_Qtc(dt);
    Add_Gradient_Term_To_Face_Qtc(dt);
}
//######################################################################
// Add_Linear_Term_To_Face_Qtc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Linear_Term_To_Face_Qtc(const T dt)
{
    const T minus_dt_over_tau_2=-dt/tau_2;
    for(int axis=0;axis<d;++axis)  
        SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),
            minus_dt_over_tau_2,face_qtc_backup_channels(axis),face_qtc_channels(axis),face_qtc_channels(axis),Topology_Helper::Face_Active_Mask(axis));
}
//######################################################################
// Add_Gradient_Term_To_Face_Qtc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Gradient_Term_To_Face_Qtc(const T dt)
{
    const T coeff=-dt*((T)1.-fc_2)/tau_2; const TV one_over_dx=hierarchy->Lattice(0).one_over_dX;
    K_Gradient_T_Helper<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),
                                        face_qtc_channels,T_backup_channel,coeff,one_over_dx);
}
//######################################################################
// Implicitly_Update_Face_Qtc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Implicitly_Update_Face_Qtc(const T dt)
{
}
//######################################################################
// Backup
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Backup()
{
    if(FICKS){Backup_Density();Backup_Temperature();}
    else{Backup_Density();Backup_Temperature();Backup_Qsc();Backup_Qtc();}
}
//######################################################################
// Backup_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Backup_Density()
{
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_backup_channel);
    SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,
                                            density_backup_channel,(unsigned)Cell_Type_Interior);
}
//######################################################################
// Backup_Temperature
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Backup_Temperature()
{
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),T_backup_channel);
    SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),T_channel,
                                        T_backup_channel,(unsigned)Cell_Type_Interior);
}
//######################################################################
// Backup_Qsc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Backup_Qsc()
{
    for(int axis=0;axis<d;++axis) { SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),face_qsc_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),face_qsc_channels(axis),
                                                face_qsc_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));}
}
//######################################################################
// Backup_Qtc
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Backup_Qtc()
{
    for(int axis=0;axis<d;++axis) {SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),face_qtc_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),face_qtc_channels(axis),
                                                face_qtc_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));}        
}
//######################################################################
// Modify_Density_With_Sources
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Modify_Density_With_Sources()
{
    Density_Modifier<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(0),density_channel,density_sources,0);
}
//######################################################################
// Add_Source
//######################################################################
// template<class T,int d> void DG_Example<T,d>::
// Add_Source(const T dt)
// {
    //     Source_Adder<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(0),density_channel,sources,source_rate,dt,0);
// }
//######################################################################
// Advect_Face_Velocities
//######################################################################
// template<class T,int d> void DG_Example<T,d>::
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
// template<class T,int d> void DG_Example<T,d>::
// Set_Neumann_Faces_Inside_Sources()
// {
//     
//         auto flags=hierarchy->Allocator(0).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
//         for(Grid_Iterator_Face<T,d> iterator(hierarchy->Lattice(0));iterator.Valid();iterator.Next()){
//             const int axis=iterator.Axis();const T_INDEX& face_index=iterator.Face_Index();
//             uint64_t face_offset=Flag_array_mask::Linear_Offset(face_index._data);
//             const uint64_t face_active_mask=Topology_Helper::Face_Active_Mask(axis);
//             const uint64_t face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
//             const TV X=hierarchy->Lattice(0).Face(axis,face_index);
//             for(size_t i=0;i<sources.size();++i) if(sources(i)->Inside(X)){
//                 if(hierarchy->template Set<unsigned>(0,&Struct_type::flags).Is_Set(face_offset,face_valid_mask))
//                     flags(face_offset)&=~face_active_mask;
//                 break;}}
// }
//######################################################################
// Initialize_Velocity_Field
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Initialize_Velocity_Field()
{
    Uniform_Velocity_Field_Initializer<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(0),face_velocity_channels,bv,0);
}
//######################################################################
// Project
//######################################################################
// template<class T,int d> void DG_Example<T,d>::
// Project(const T dt)
// {
//     using Multigrid_struct_type             = Multigrid_Data<T>;
//     using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;

//     // set up divergence channel
//     T Struct_type::* divergence_channel     = &Struct_type::ch8;

//     // clear
//     
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),pressure_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),divergence_channel);

//     // set boundary conditions
//     
//         Boundary_Condition_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(0),face_velocity_channels,pressure_channel,bv,0);

//     // compute divergence
//     Hierarchy_Projection::Compute_Divergence(*hierarchy,face_velocity_channels,divergence_channel);

//     Poisson_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,3,1,200);

//     T Struct_type::* q_channel              = &Struct_type::ch9;
//     T Struct_type::* r_channel              = &Struct_type::ch10;
//     T Struct_type::* s_channel              = &Struct_type::ch11;
//     T Struct_type::* k_channel              = &Struct_type::ch11;
//     T Struct_type::* z_channel              = &Struct_type::ch12;

//     // clear all channels
//     
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),q_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),r_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),s_channel);
//         SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),z_channel);

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

//     
//         Apply_Pressure<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(0),face_velocity_channels,pressure_channel,0);
// }
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Register_Options()
{
    Base::Register_Options();
    // for phase field
    parse_args->Add_Option_Argument("-ficks","Fick's diffusion.");
    parse_args->Add_Double_Argument("-cdv",(T)1.,"Constant density value.");
    parse_args->Add_Double_Argument("-taus",(T)3.e-4,"tau s.");
    parse_args->Add_Double_Argument("-tau1",(T).0001,"tau 1.");
    parse_args->Add_Double_Argument("-fc1",(T)0.,"Fc 1.");
    parse_args->Add_Double_Argument("-tau2",(T).0001,"tau 2.");
    parse_args->Add_Double_Argument("-fc2",(T)0.,"Fc 2.");
    parse_args->Add_Double_Argument("-gamma",(T)10.,"gamma.");
    parse_args->Add_Double_Argument("-delta",(T).2,"delta.");
    parse_args->Add_Double_Argument("-eps",(T)1.e-2,"epsilon.");
    parse_args->Add_Double_Argument("-K",(T)2.,"K.");
    parse_args->Add_Double_Argument("-bv",(T)1.,"background velocity.");

    parse_args->Add_Double_Argument("-ma",(T).9,"m_alpha.");
    parse_args->Add_Double_Argument("-SR",(T)0.,"SR.");
    parse_args->Add_Double_Argument("-rf",(T).05,"random factor.");

    parse_args->Add_Option_Argument("-ed","Explicit diffusion.");
    parse_args->Add_Option_Argument("-uvf","Uniform velocity field.");

    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL number.");
    parse_args->Add_Double_Argument("-dt",(T)1.e-6,"Constant time step.");
    parse_args->Add_Integer_Argument("-test_number",1,"Test number.");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(200.),"n","Grid resolution.");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(200.),"n","Grid resolution.");
    parse_args->Add_Double_Argument("-cw",2e-2,"cell width.");

    // for CG
    parse_args->Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args->Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args->Add_Double_Argument("-cg_tolerance",1e-4,"CG tolerance.");
}
//######################################################################
// Parse_Options
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Parse_Options()
{
    Base::Parse_Options();

    const_density_value=parse_args->Get_Double_Value("-cdv");
    tau_s=parse_args->Get_Double_Value("-taus");
    tau_1=parse_args->Get_Double_Value("-tau1");
    fc_1=parse_args->Get_Double_Value("-fc1");
    tau_2=parse_args->Get_Double_Value("-tau2");
    fc_2=parse_args->Get_Double_Value("-fc2");
    gamma=parse_args->Get_Double_Value("-gamma");
    delta=parse_args->Get_Double_Value("-delta");
    K=parse_args->Get_Double_Value("-K");
    bv=parse_args->Get_Double_Value("-bv");
    epsilon=parse_args->Get_Double_Value("-eps");
    m_alpha=parse_args->Get_Double_Value("-ma");
    random_factor=parse_args->Get_Double_Value("-rf");
    test_number=parse_args->Get_Integer_Value("-test_number");
    FICKS=parse_args->Get_Option_Value("-ficks");
    explicit_diffusion=parse_args->Get_Option_Value("-ed");
    uvf=parse_args->Get_Option_Value("-uvf");

    cfl=(T)parse_args->Get_Double_Value("-cfl");
    cell_width=(T)parse_args->Get_Double_Value("-cw");
    const_time_step=(T)parse_args->Get_Double_Value("-dt");
    levels=parse_args->Get_Integer_Value("-levels");
    mg_levels=parse_args->Get_Integer_Value("-mg_levels");
    number_of_threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(number_of_threads);
    if(d==2){auto cell_counts_2d=parse_args->Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_2d(v);}
    else{auto cell_counts_3d=parse_args->Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_3d(v);}
    cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    cg_restart_iterations=parse_args->Get_Integer_Value("-cg_restart_iterations");
    cg_tolerance=(T)parse_args->Get_Double_Value("-cg_tolerance");
    if(FICKS){fc_1=(T)1.; fc_2=(T)1.;}
}
//######################################################################
// Log_Parameters
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Log_Parameters() const
{
    Base::Log_Parameters();
    *output<<"tau_s="<<tau_s<<std::endl;
    *output<<"tau_1="<<tau_1<<std::endl;
    *output<<"fc_1="<<fc_1<<std::endl;
    *output<<"tau_2="<<tau_2<<std::endl;
    *output<<"fc_2="<<fc_2<<std::endl;
    *output<<"gamma="<<gamma<<std::endl;
    *output<<"K="<<K<<std::endl;
    *output<<"delta="<<delta<<std::endl;
    *output<<"epsilon="<<epsilon<<std::endl;
    *output<<"m_alpha="<<m_alpha<<std::endl;
    *output<<"SR="<<SR<<std::endl;
    *output<<"rf="<<random_factor<<std::endl;
    *output<<"dt="<<const_time_step<<std::endl;
    *output<<"dt="<<cell_width<<std::endl;
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Write_Output_Files(const int frame) const
{
    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/frame_title",frame_title);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_density",density_channel);
}
//######################################################################
// Read_Output_Files
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Read_Output_Files(const int frame)
{
}
//######################################################################
template class Nova::DG_Example<float,2>;
template class Nova::DG_Example<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::DG_Example<double,2>;
template class Nova::DG_Example<double,3>;
#endif
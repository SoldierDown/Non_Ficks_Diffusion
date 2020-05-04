//!#####################################################################
//! \file DG_Example.cpp
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
    gamma=(T)10.; Teq=(T)1.;
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
    // Add_Random_Term_To_Density(dt);
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
    
    // using Hierarchy_Projection                          = Grid_Hierarchy_Projection<Struct_type,T,d>;
    // using Flags_type                                    = typename Struct_type::Flags_type;
    // using T_Range_Iterator                              = Range_Iterator<d,T_INDEX>;

    // if(FICKS) fc_1=(T)1.;
    // const T one_over_dx2=Nova_Utilities::Sqr(hierarchy->Lattice(0).one_over_dX(0)); const T dt_over_tau_s=dt/tau_s;

    // auto blocks=hierarchy->Blocks(0);
    // auto block_size=hierarchy->Allocator(0).Block_Size();
    // auto density_backup_data=hierarchy->Allocator(0).template Get_Const_Array<Struct_type,T>(density_backup_channel);
    // auto density_data=hierarchy->Allocator(0).template Get_Array<Struct_type,T>(density_channel);
    // auto flags=hierarchy->Allocator(0).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
    
    // for(unsigned block=0;block<blocks.second;++block){
    //     uint64_t offset=blocks.first[block];
    //     Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
    //     T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
    //     for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
    //         if(hierarchy->template Set<unsigned>(0,&Struct_type::flags).Is_Set(offset,Cell_Type_Interior)){
    //             const T_INDEX index=base_index+range_iterator.Index(); T laplacian=(T)0.;
    //             for(T_Range_Iterator iterator(T_INDEX(-1),T_INDEX(1));iterator.Valid();iterator.Next()){
    //                 if(iterator.Index()==T_INDEX()) ; //skip itself
    //                 else{T_INDEX neighbor_index=index+iterator.Index(); uint64_t neighbor_offset(Flag_array_mask::Linear_Offset(neighbor_index._data));
    //                     if((hierarchy->template Set<unsigned>(0,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Interior))
    //                         ||(hierarchy->template Set<unsigned>(0,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Dirichlet))){
    //                         laplacian+=density_backup_data(neighbor_offset)-density_backup_data(offset);}}}
    //                         laplacian*=one_third*one_over_dx2*fc_1*Nova_Utilities::Sqr(epsilon)*dt_over_tau_s;
    //                         density_data(offset)+=laplacian;} range_iterator.Next();}}
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
    DG_M_Calculator<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),m_channel,sigma_channel,T_backup_channel,m_alpha,gamma,Teq);
    Range_Checker<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),m_channel,(T)-.5,(T).5);
    Density_Plus_Poly_Term_Helper<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,density_backup_channel,m_channel,dt_over_tau_s);    
}
//######################################################################
// Add_Random_Term_To_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Add_Random_Term_To_Density(const T dt)
{
    const T random_factor=(T).2;
    const T random_value=random.Get_Uniform_Number(0,1);
    Add_Random_Term<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel,density_backup_channel,random_factor,random_value);
}
//######################################################################
// Implicitly_Update_Density
//######################################################################
template<class T,int d> void DG_Example<T,d>::
Implicitly_Update_Density(const T dt)
{
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
    // add K*dS to T
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

    // using Hierarchy_Projection                          = Grid_Hierarchy_Projection<Struct_type,T,d>;
    // using Flags_type                                    = typename Struct_type::Flags_type;
    // using T_Range_Iterator                              = Range_Iterator<d,T_INDEX>;

    // if(FICKS) fc_2=(T)1.;
    // const T one_over_dx2=Nova_Utilities::Sqr(hierarchy->Lattice(0).one_over_dX(0)); const T dt_over_tau_s=dt/tau_s;

    // auto blocks=hierarchy->Blocks(0);
    // auto block_size=hierarchy->Allocator(0).Block_Size();
    // auto T_backup_data=hierarchy->Allocator(0).template Get_Const_Array<Struct_type,T>(T_backup_channel);
    // auto T_data=hierarchy->Allocator(0).template Get_Array<Struct_type,T>(T_channel);
    // auto flags=hierarchy->Allocator(0).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
    
    // for(unsigned block=0;block<blocks.second;++block){
    //     uint64_t offset=blocks.first[block];
    //     Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
    //     T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
    //     for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
    //         if(hierarchy->template Set<unsigned>(0,&Struct_type::flags).Is_Set(offset,Cell_Type_Interior)){
    //             const T_INDEX index=base_index+range_iterator.Index(); T laplacian=(T)0.;
    //             for(T_Range_Iterator iterator(T_INDEX(-1),T_INDEX(1));iterator.Valid();iterator.Next()){
    //                 if(iterator.Index()==T_INDEX()) ; //skip itself
    //                 else{T_INDEX neighbor_index=index+iterator.Index(); uint64_t neighbor_offset(Flag_array_mask::Linear_Offset(neighbor_index._data));
    //                     if((hierarchy->template Set<unsigned>(0,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Interior))
    //                         ||(hierarchy->template Set<unsigned>(0,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Dirichlet))){
    //                         laplacian+=T_backup_data(neighbor_offset)-T_backup_data(offset);}}}
    //                         laplacian*=one_third*one_over_dx2*fc_2*k*dt;
    //                         T_data(offset)+=laplacian;} range_iterator.Next();}}
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
    for(int level=0;level<levels;++level)
        Density_Modifier<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,density_sources,level);
}
//######################################################################
// Add_Source
//######################################################################
// template<class T,int d> void DG_Example<T,d>::
// Add_Source(const T dt)
// {
    // for(int level=0;level<levels;++level)
    //     Source_Adder<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,sources,source_rate,dt,level);
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
//     for(int level=0;level<levels;++level){
//         auto flags=hierarchy->Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
//         for(Grid_Iterator_Face<T,d> iterator(hierarchy->Lattice(level));iterator.Valid();iterator.Next()){
//             const int axis=iterator.Axis();const T_INDEX& face_index=iterator.Face_Index();
//             uint64_t face_offset=Flag_array_mask::Linear_Offset(face_index._data);
//             const uint64_t face_active_mask=Topology_Helper::Face_Active_Mask(axis);
//             const uint64_t face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
//             const TV X=hierarchy->Lattice(level).Face(axis,face_index);
//             for(size_t i=0;i<sources.size();++i) if(sources(i)->Inside(X)){
//                 if(hierarchy->template Set<unsigned>(level,&Struct_type::flags).Is_Set(face_offset,face_valid_mask))
//                     flags(face_offset)&=~face_active_mask;
//                 break;}}}
// }
//######################################################################
// Initialize_Velocity_Field
//######################################################################
// template<class T,int d> void DG_Example<T,d>::
// Initialize_Velocity_Field()
// {
//     for(int level=0;level<levels;++level)
//         Uniform_Velocity_Field_Initializer<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,bv,level);
// }
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
template<class T,int d> void DG_Example<T,d>::
Register_Options()
{
    Base::Register_Options();
    // for phase field
    parse_args->Add_Option_Argument("-ficks","Fick's diffusion");
    parse_args->Add_Double_Argument("-cdv",(T)1.,"Constant density value.");
    parse_args->Add_Double_Argument("-taus",(T)3.e-4,"tau s.");
    parse_args->Add_Double_Argument("-tau1",(T).0001,"tau 1.");
    parse_args->Add_Double_Argument("-fc1",(T)0.,"Fc 1.");
    parse_args->Add_Double_Argument("-tau2",(T).0001,"tau 2.");
    parse_args->Add_Double_Argument("-fc2",(T)0.,"Fc 2.");
    parse_args->Add_Double_Argument("-Teq",(T)1.,"Teq.");
    parse_args->Add_Double_Argument("-gamma",(T)40.,"gamma.");
    parse_args->Add_Double_Argument("-delta",(T).2,"delta.");
    parse_args->Add_Double_Argument("-eps",(T)1.e-2,"epsilon.");

    parse_args->Add_Double_Argument("-K",(T)4.,"K.");
    parse_args->Add_Double_Argument("-ma",(T).9,"m_alpha.");
    parse_args->Add_Double_Argument("-SR",(T)0.,"SR.");

    parse_args->Add_Option_Argument("-ed","Explicit diffusion");
    parse_args->Add_Option_Argument("-uvf","Uniform velocity field.");

    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL number.");
    parse_args->Add_Double_Argument("-dt",(T)1.e-6,"Constant time step.");
    parse_args->Add_Integer_Argument("-test_number",1,"Test number.");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(200.),"n","Grid resolution");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(200.),"n","Grid resolution");

    
    // for CG
    parse_args->Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args->Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args->Add_Double_Argument("-cg_tolerance",1e-4,"CG tolerance");
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
    epsilon=parse_args->Get_Double_Value("-eps");
    K=parse_args->Get_Double_Value("-K");
    m_alpha=parse_args->Get_Double_Value("-ma");

    FICKS=parse_args->Get_Option_Value("-ficks");
    explicit_diffusion=parse_args->Get_Option_Value("-ed");
    uvf=parse_args->Get_Option_Value("-uvf");


    cfl=(T)parse_args->Get_Double_Value("-cfl");
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
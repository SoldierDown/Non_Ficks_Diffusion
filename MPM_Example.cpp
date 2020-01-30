//!#####################################################################
//! \file MPM_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Grids/Grid_Iterator_Cell.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Utilities.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <omp.h>
#include "MPM_Example.h"
#include "Traverse_Helper.h"
#include "Velocity_Normalization_Helper.h"
#include "Explicit_Force_Helper.h"
#include "Grid_Based_Collision_Helper.h"
#include "Flag_Helper.h"
#include "MPM_RHS_Helper.h"
#include "Channel_Vector_Norm_Helper.h"
#include "Compare_Helper.h"
#include "MPM_Flags.h"

#include <nova/SPGrid/Tools/SPGrid_Arithmetic.h>

#include "./Implicit_Force_Helper/MPM_CG_Vector.h"
#include "./Implicit_Force_Helper/MPM_CG_System.h"


using namespace Nova;
using namespace SPGrid;
//######################################################################
// Constructor
//######################################################################
template<class T,int d> MPM_Example<T,d>::
MPM_Example()
    :Base(),hierarchy(nullptr)
{
    solver_tolerance=(T)1e-7;
    solver_iterations=10000;

    gravity=TV::Axis_Vector(1)*(T)-2.;
    flip=(T).9;
     
    flags_channel                           = &Struct_type::flags;
    mass_channel                            = &Struct_type::ch0;
    velocity_channels(0)                    = &Struct_type::ch1;
    velocity_channels(1)                    = &Struct_type::ch2;
    if(d==3) velocity_channels(2)           = &Struct_type::ch3;
    velocity_star_channels(0)               = &Struct_type::ch4;
    velocity_star_channels(1)               = &Struct_type::ch5;
    if(d==3) velocity_star_channels(2)      = &Struct_type::ch6;
    f_channels(0)                           = &Struct_type::ch7;
    f_channels(1)                           = &Struct_type::ch8;
    if(d==3) f_channels(2)                  = &Struct_type::ch9;

    collide_nodes_channel                   = &Struct_type::ch10;

    rhs_channels(0)                         = &Struct_type::ch11;
    rhs_channels(1)                         = &Struct_type::ch12;
    if(d==3) rhs_channels(2)                = &Struct_type::ch13;

    q_channels(0)                           = &Struct_type::ch14;
    q_channels(1)                           = &Struct_type::ch15;
    if(d==3) q_channels(2)                  = &Struct_type::ch16;

    s_channels(0)                           = &Struct_type::ch17;
    s_channels(1)                           = &Struct_type::ch18;
    if(d==3) s_channels(2)                  = &Struct_type::ch19;

    r_channels(0)                           = &Struct_type::ch20;
    r_channels(1)                           = &Struct_type::ch21;
    if(d==3) r_channels(2)                  = &Struct_type::ch22;

    z_channels(0)                           = &Struct_type::ch23;
    z_channels(1)                           = &Struct_type::ch24;
    if(d==3) z_channels(2)                  = &Struct_type::ch25;


}
//######################################################################
// Destructor
//######################################################################
template<class T,int d> MPM_Example<T,d>::
~MPM_Example()
{
    if(hierarchy!=nullptr) delete hierarchy;
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Initialize()
{
    Initialize_Particles(Base::test_number);
    Populate_Simulated_Particles();
    Initialize_SPGrid();
    particle_bins.Resize(threads,threads);
    Log::cout<<"barrier size: "<<barriers.size()<<std::endl;
}
//######################################################################
// Reset_Grid_Based_Variables
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Reset_Grid_Based_Variables()
{
}
//######################################################################
// Reset_Solver_Channels
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Reset_Solver_Channels()
{
    for(int level=0;level<levels;++level) for(int v=0;v<d;++v){
        Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channels(v));
        Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channels(v));
        Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channels(v));
        Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channels(v));}

}
//######################################################################
// Compute_Bounding_Box
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Compute_Bounding_Box(Range<T,d>& bbox)
{
    Array<TV> min_corner_per_thread(threads,domain.max_corner);
    Array<TV> max_corner_per_thread(threads,domain.min_corner);

#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int tid=omp_get_thread_num();
        const int id=simulated_particles(i);
        T_Particle p=particles(id);
        TV& current_min_corner=min_corner_per_thread(tid);
        TV& current_max_corner=max_corner_per_thread(tid);
        for(int v=0;v<d;++v){
            T dd=(T)3./counts(v);
            current_min_corner(v)=std::min(current_min_corner(v),p.X(v)-dd);
            current_max_corner(v)=std::max(current_max_corner(v),p.X(v)+dd);}}

    for(int v=0;v<d;++v){
        bbox.min_corner(v)=min_corner_per_thread(0)(v);
        bbox.max_corner(v)=max_corner_per_thread(0)(v);
    }

    for(int tid=1;tid<threads;++tid) for(int v=0;v<d;++v){
        bbox.min_corner(v)=std::min(bbox.min_corner(v),min_corner_per_thread(tid)(v));
        bbox.max_corner(v)=std::max(bbox.max_corner(v),max_corner_per_thread(tid)(v));}
    
    for(int v=0;v<d;++v){ 
        bbox.min_corner(v)=std::max((T)0.,bbox.min_corner(v));
        bbox.max_corner(v)=std::min((T)1.,bbox.max_corner(v));}
}
//######################################################################
// Rasterize_Voxels
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Rasterize_Voxels(const Range<T,d>& bbox)
{
    using Cell_Iterator             = Grid_Iterator_Cell<T,d>;

    const Grid<T,d>& grid=hierarchy->Lattice(0);
    Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));

    for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next())
        hierarchy->Activate_Cell(0,iterator.Cell_Index(),Cell_Type_Interior);
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Initialize_SPGrid()
{
    using Hierarchy_Initializer     = Grid_Hierarchy_Initializer<Struct_type,T,d>;

    Compute_Bounding_Box(bbox);

    if(hierarchy!=nullptr) delete hierarchy;
    hierarchy=new Hierarchy(counts,domain,levels);

    Rasterize_Voxels(bbox);

    Hierarchy_Initializer::Flag_Ghost_Cells(*hierarchy);
    Hierarchy_Initializer::Flag_Valid_Faces(*hierarchy);
    Hierarchy_Initializer::Flag_Active_Faces(*hierarchy);
    Hierarchy_Initializer::Flag_Active_Nodes(*hierarchy);
    Hierarchy_Initializer::Flag_Shared_Nodes(*hierarchy);
    Hierarchy_Initializer::Flag_Ghost_Nodes(*hierarchy);
    Hierarchy_Initializer::Flag_T_Junction_Nodes(*hierarchy);
    hierarchy->Update_Block_Offsets();
    hierarchy->Initialize_Red_Black_Partition(2*threads);
}
//######################################################################
// Populated_Simulated_Particles
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Populate_Simulated_Particles()
{
    simulated_particles.Clear();
    for(int i=0;i<particles.size();++i)
        if(particles(i).valid)
            simulated_particles.Append(i);
    Log::cout<<"\nSimulated particles: "<<simulated_particles.size()<<std::endl;
}
//######################################################################
// Max_Particle_Velocity
//######################################################################
template<class T,int d> T MPM_Example<T,d>::
Max_Particle_Velocity() const
{
    Array<T> result_per_thread(threads);
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int tid=omp_get_thread_num();
        T& r=result_per_thread(tid); const int id=simulated_particles(i);
        r=std::max(r,particles(id).V.Norm_Squared());}
    T result=(T)0.;
    for(int tid=0;tid<threads;++tid) result=std::max(result,result_per_thread(tid));
    std::cout<<"max v: "<<std::sqrt(result)<<std::endl;
    return std::sqrt(result);
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Limit_Dt(T& dt,const T time)
{
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Update_Particle_Weights()
{
    const Grid<T,d>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
        for(int i=0;i<simulated_particles.size();i++){
            const int id=simulated_particles(i);
            T_Particle& p=particles(id);
            p.Update_Weights(grid);}
}
//######################################################################
// Group_Particles
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Group_Particles()
{
    const Grid<T,d>& grid=hierarchy->Lattice(0);
    Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));
    T_INDEX min_corner=grid.Clamp_To_Cell(bbox.min_corner),max_corner=grid.Clamp_To_Cell(bbox.max_corner)+1;
    x_intervals.resize(threads); x_intervals(0).min_corner=min_corner(0); x_intervals(threads-1).max_corner=max_corner(0);
    for(int i=1;i<threads;++i){
        int n=min_corner(0)+(max_corner(0)-min_corner(0)+1)*i/threads;
        x_intervals(i-1).max_corner=n-1;
        x_intervals(i).min_corner=n;}
    for(int i=0;i<threads;++i) for(int j=0;j<threads;++j) particle_bins(i,j).Clear();
#pragma omp parallel for
    for(int i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i);
        const T_INDEX& closest_node=particles(id).closest_node;
        const int tid_collect=omp_get_thread_num();
        const Interval<int> particle_x_interval=Interval<int>(closest_node(0)-1,closest_node(0)+1);
        for(int tid_process=0;tid_process<threads;++tid_process){
            const Interval<int> thread_x_interval=x_intervals(tid_process);
            if(particle_x_interval.Intersection(thread_x_interval)) particle_bins(tid_process,tid_collect).Append(id);}
    }
    
    // for(int tid_process=0;tid_process<threads;++tid_process){
    //     for(int tid_collect=0;tid_collect<threads;++tid_collect){
    //         const Array<int> index=particle_bins(tid_process,tid_collect);
    //         for(int i=0;i<index.size();++i) Log::cout<<i<<": "<<index(i)<<std::endl;
    //     }
    // }
}
//######################################################################
// Rasterize
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Rasterize()
{
    Log::cout.precision(10);
    auto mass=hierarchy->Channel(0,mass_channel); auto flags=hierarchy->Channel(0,flags_channel);
    auto collide_nodes=hierarchy->Channel(0,collide_nodes_channel);
    const Grid<T,d>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(int tid_process=0;tid_process<threads;++tid_process){
        const Interval<int> thread_x_interval=x_intervals(tid_process);
        for(int tid_collect=0;tid_collect<threads;++tid_collect){
            const Array<int> index=particle_bins(tid_process,tid_collect);
            for(int i=0;i<index.size();++i){
                T_Particle& p=particles(index(i));T_INDEX& closest_node=p.closest_node;
                const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_node(0),thread_x_interval.max_corner-closest_node(0));
        for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){T_INDEX current_node=iterator.Current_Node();
            T weight=iterator.Weight(); TV current_node_location=grid.Node(current_node);
            if(weight>(T)0.){                                    
                mass(current_node._data)+=weight*p.mass; flags(current_node._data)|=Node_Saturated;
                for(int id=0;id<barriers.size();++id) if((current_node_location-barriers(id).surface).Dot_Product(barriers(id).normal)<(T)0.) collide_nodes(current_node._data)=(T)1.;                 
                for(int v=0;v<d;++v) hierarchy->Channel(0,velocity_channels(v))(current_node._data)+=weight*p.mass*p.V(v);}}}}}

    // normalize weights for velocity (to conserve momentum)
    for(int level=0;level<levels;++level) Velocity_Normalization_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_channels);     
}
//######################################################################
// Update_Constitutive_Model_State
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Update_Constitutive_Model_State()
{
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i); 
        T_Particle &particle=particles(id);    
        particle.constitutive_model.Precompute();}
}
//######################################################################
// Update_Particle_Velocities_And_Positions
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Update_Particle_Velocities_And_Positions(const T dt)
{
    Apply_Force(dt);
    const Grid<T,d>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i); 
        T_Particle &p=particles(id);
        TV V_pic=TV(),V_flip=p.V; 
        Matrix<T,d> grad_Vp=Matrix<T,d>();
        for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){
            T_INDEX current_node=iterator.Current_Node(); T weight=iterator.Weight();
            if(weight>(T)0.){
                TV V_grid,delta_V_grid;
                for(int v=0;v<d;++v){
                    V_grid(v)=hierarchy->Channel(0,velocity_star_channels(v))(current_node._data); 
                    delta_V_grid(v)=hierarchy->Channel(0,velocity_star_channels(v))(current_node._data)-hierarchy->Channel(0,velocity_channels(v))(current_node._data);}
                V_pic+=weight*V_grid; 
                V_flip+=weight*delta_V_grid;
                grad_Vp+=Matrix<T,d>::Outer_Product(V_grid,iterator.Weight_Gradient());}}
            p.constitutive_model.Fe+=dt*grad_Vp*p.constitutive_model.Fe;
            p.V=V_flip*flip+V_pic*((T)1.-flip);
            p.X+=V_pic*dt;
        if(!grid.domain.Inside(p.X)) p.valid=false;
    }        
}
//######################################################################
// Apply_Force
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Apply_Force(const T dt)
{
    Apply_Explicit_Force(dt);
    Grid_Based_Collison(true);
    if(true){
    Conjugate_Gradient<T> cg;
    Krylov_Solver<T>* solver=(Krylov_Solver<T>*)&cg;
    MPM_CG_System<Struct_type,T,d> mpm_system(*hierarchy,simulated_particles,particles,particle_bins,x_intervals,barriers,(T)0.,dt,threads);
    mpm_system.use_preconditioner=false;
    // set rhs here
    for(int level=0;level<levels;++level) MPM_RHS_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_star_channels,rhs_channels);     
    // clear channels
    Reset_Solver_Channels();
    MPM_CG_Vector<Struct_type,T,d> solver_vp(*hierarchy,velocity_star_channels),solver_rhs(*hierarchy,rhs_channels),solver_q(*hierarchy,q_channels),solver_s(*hierarchy,s_channels),solver_r(*hierarchy,r_channels),solver_k(*hierarchy,z_channels),solver_z(*hierarchy,z_channels);
    solver->Solve(mpm_system,solver_vp,solver_rhs,solver_q,solver_s,solver_r,solver_k,solver_z,solver_tolerance,0,solver_iterations);}

}
//######################################################################
// Apply_Explicit_Force
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Apply_Explicit_Force(const T dt)
{
    const Grid<T,d>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(int tid_process=0;tid_process<threads;++tid_process){
        const Interval<int>& thread_x_interval=x_intervals(tid_process);
        for(int tid_collect=0;tid_collect<threads;++tid_collect){
            const Array<int> index=particle_bins(tid_process,tid_collect);
            for(int i=0;i<index.size();++i){
                T_Particle& p=particles(index(i));T_INDEX& closest_node=p.closest_node;
                Matrix<T,d> P=p.constitutive_model.P(),F=p.constitutive_model.Fe; T V0=p.volume;
                Matrix<T,d> V0_P_FT=P.Times_Transpose(F)*V0;                    
                V0_P_FT=P.Times_Transpose(F)*V0;
                const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_node(0),thread_x_interval.max_corner-closest_node(0));
                for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){auto data=iterator.Current_Node()._data; T weight=iterator.Weight();                            
                // if(weight>(T)0.) 
                {TV weight_grad=iterator.Weight_Gradient();
                for(int v=0;v<d;++v)  hierarchy->Channel(0,f_channels(v))(data)+=gravity(v)*p.mass*weight-(V0_P_FT*weight_grad)(v);}}}}}

    for(int level=0;level<levels;++level) Explicit_Force_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),f_channels,velocity_channels,velocity_star_channels,dt);
}
//######################################################################
// Grid_Based_Collision
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Grid_Based_Collison(const bool detect_collision)
{   
    if(detect_collision) for(int level=0;level<levels;++level) Grid_Based_Collision_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_star_channels,collide_nodes_channel,barriers);
}
//######################################################################
// Estimate_Particle_Volumes
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Estimate_Particle_Volumes()
{   
    auto mass=hierarchy->Channel(0,mass_channel);
    const Grid<T,d>& grid=hierarchy->Lattice(0); const T one_over_volume_per_cell=(T)1./grid.dX.Product();
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int id=simulated_particles(i); T_Particle& p=particles(id);     
        T particle_density=(T)0.;
        for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){
            particle_density+=iterator.Weight()*mass(iterator.Current_Node()._data);}
        particle_density*=one_over_volume_per_cell;
        p.volume=p.mass/particle_density;} 
}
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Register_Options()
{
    Base::Register_Options();
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Double_Argument("-cfl",(T)0.1,"CFL number.");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(32.),"n","Grid resolution");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(32.),"n","Grid resolution");
}
//######################################################################
// Test
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Test()
{

}
//######################################################################
// Parse_Options
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Parse_Options()
{
    Base::Parse_Options();
    threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(threads);
    Base::test_number=parse_args->Get_Integer_Value("-test_number");
    cfl=parse_args->Get_Double_Value("-cfl");
    levels=parse_args->Get_Integer_Value("-levels");
    if(d==2){auto cell_counts_2d=parse_args->Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_2d(v);}
    else{auto cell_counts_3d=parse_args->Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_3d(v);}
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Write_Output_Files(const int frame) const
{
    if(frame==first_frame){std::string deformables_filename=(d==2)?output_directory+"/common/metadata.mpm2d":output_directory+"/common/metadata.mpm3d";
        File_Utilities::Write_To_Text_File(deformables_filename,std::to_string(frame));}

    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/frame_title",frame_title);

    File_Utilities::Write_To_File(output_directory+"/"+std::to_string(frame)+"/particles",particles);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);
}
//######################################################################
// Read_Output_Files
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Read_Output_Files(const int frame)
{
    File_Utilities::Read_From_File(output_directory+"/"+std::to_string(frame)+"/particles",particles);
}
//######################################################################
template class Nova::MPM_Example<float,2>;
template class Nova::MPM_Example<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::MPM_Example<double,2>;
template class Nova::MPM_Example<double,3>;
#endif

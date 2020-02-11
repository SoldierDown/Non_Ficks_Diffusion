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
#include "Flag_Setup_Helper.h"


#include <nova/SPGrid/Tools/SPGrid_Arithmetic.h>

#include "./Implicit_Force_Helper/MPM_CG_Vector.h"
#include "./Implicit_Force_Helper/MPM_CG_System.h"

using namespace Nova;
using namespace SPGrid;
//######################################################################
// Constructor
//######################################################################
template<class T> MPM_Example<T,2>::
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
    velocity_star_channels(0)               = &Struct_type::ch4;
    velocity_star_channels(1)               = &Struct_type::ch5;
    f_channels(0)                           = &Struct_type::ch7;
    f_channels(1)                           = &Struct_type::ch8;

    rhs_channels(0)                         = &Struct_type::ch10;
    rhs_channels(1)                         = &Struct_type::ch11;

    q_channels(0)                           = &Struct_type::ch13;
    q_channels(1)                           = &Struct_type::ch14;

    s_channels(0)                           = &Struct_type::ch16;
    s_channels(1)                           = &Struct_type::ch17;

    r_channels(0)                           = &Struct_type::ch19;
    r_channels(1)                           = &Struct_type::ch20;

    z_channels(0)                           = &Struct_type::ch22;
    z_channels(1)                           = &Struct_type::ch23;
}
//######################################################################
// Constructor
//######################################################################
template<class T> MPM_Example<T,3>::
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
    velocity_channels(2)                    = &Struct_type::ch3;
    
    velocity_star_channels(0)               = &Struct_type::ch4;
    velocity_star_channels(1)               = &Struct_type::ch5;
    velocity_star_channels(2)               = &Struct_type::ch6;

    f_channels(0)                           = &Struct_type::ch7;
    f_channels(1)                           = &Struct_type::ch8;
    f_channels(2)                           = &Struct_type::ch9;

    rhs_channels(0)                         = &Struct_type::ch10;
    rhs_channels(1)                         = &Struct_type::ch11;
    rhs_channels(2)                         = &Struct_type::ch12;

    q_channels(0)                           = &Struct_type::ch13;
    q_channels(1)                           = &Struct_type::ch14;
    q_channels(2)                           = &Struct_type::ch15;

    s_channels(0)                           = &Struct_type::ch16;
    s_channels(1)                           = &Struct_type::ch17;
    s_channels(2)                           = &Struct_type::ch18;

    r_channels(0)                           = &Struct_type::ch19;
    r_channels(1)                           = &Struct_type::ch20;
    r_channels(2)                           = &Struct_type::ch21;

    z_channels(0)                           = &Struct_type::ch22;
    z_channels(1)                           = &Struct_type::ch23;
    z_channels(2)                           = &Struct_type::ch24;
}
//######################################################################
// Destructor
//######################################################################
template<class T> MPM_Example<T,2>::
~MPM_Example()
{
    if(hierarchy!=nullptr) delete hierarchy;
}
//######################################################################
// Destructor
//######################################################################
template<class T> MPM_Example<T,3>::
~MPM_Example()
{
    if(hierarchy!=nullptr) delete hierarchy;
}
//######################################################################
// Initialize
//######################################################################
template<class T> void MPM_Example<T,2>::
Initialize()
{
    Initialize_Particles(Base::test_number);
    Populate_Simulated_Particles();
    Initialize_SPGrid();
    particle_bins.Resize(threads,threads);
    Log::cout<<"barrier size: "<<barriers.size()<<std::endl;
}
//######################################################################
// Initialize
//######################################################################
template<class T> void MPM_Example<T,3>::
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
template<class T> void MPM_Example<T,2>::
Reset_Grid_Based_Variables()
{
}
//######################################################################
// Reset_Grid_Based_Variables
//######################################################################
template<class T> void MPM_Example<T,3>::
Reset_Grid_Based_Variables()
{
}
//######################################################################
// Reset_Solver_Channels
//######################################################################
template<class T> void MPM_Example<T,2>::
Reset_Solver_Channels()
{
    for(int level=0;level<levels;++level) for(int v=0;v<2;++v){
        Clear<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channels(v));
        Clear<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channels(v));
        Clear<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channels(v));
        Clear<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channels(v));}
}
//######################################################################
// Reset_Solver_Channels
//######################################################################
template<class T> void MPM_Example<T,3>::
Reset_Solver_Channels()
{
    for(int level=0;level<levels;++level) for(int v=0;v<3;++v){
        Clear<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channels(v));
        Clear<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channels(v));
        Clear<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channels(v));
        Clear<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channels(v));}
}
//######################################################################
// Compute_Bounding_Box
//######################################################################
template<class T> void MPM_Example<T,2>::
Compute_Bounding_Box(Range<T,2>& bbox)
{
    Array<TV> min_corner_per_thread(threads,domain.max_corner);
    Array<TV> max_corner_per_thread(threads,domain.min_corner);

#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int tid=omp_get_thread_num();
        const int id=simulated_particles(i);
        T_Particle p=particles(id);
        TV& current_min_corner=min_corner_per_thread(tid);
        TV& current_max_corner=max_corner_per_thread(tid);
        for(int v=0;v<2;++v){
            const T dd=(T)2./counts(v);
            current_min_corner(v)=std::min(current_min_corner(v),p.X(v)-dd);
            current_max_corner(v)=std::max(current_max_corner(v),p.X(v)+dd);}}

    for(int v=0;v<2;++v){
        bbox.min_corner(v)=min_corner_per_thread(0)(v);
        bbox.max_corner(v)=max_corner_per_thread(0)(v);
    }

    for(int tid=1;tid<threads;++tid) for(int v=0;v<2;++v){
        bbox.min_corner(v)=std::min(bbox.min_corner(v),min_corner_per_thread(tid)(v));
        bbox.max_corner(v)=std::max(bbox.max_corner(v),max_corner_per_thread(tid)(v));}
    
    for(int v=0;v<2;++v){ 
        bbox.min_corner(v)=std::max((T)0.,bbox.min_corner(v));
        bbox.max_corner(v)=std::min((T)1.,bbox.max_corner(v));}
}
//######################################################################
// Compute_Bounding_Box
//######################################################################
template<class T> void MPM_Example<T,3>::
Compute_Bounding_Box(Range<T,3>& bbox)
{
    Array<TV> min_corner_per_thread(threads,domain.max_corner);
    Array<TV> max_corner_per_thread(threads,domain.min_corner);

#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int tid=omp_get_thread_num();
        const int id=simulated_particles(i);
        T_Particle p=particles(id);
        TV& current_min_corner=min_corner_per_thread(tid);
        TV& current_max_corner=max_corner_per_thread(tid);
        for(int v=0;v<3;++v){
            const T dd=(T)2./counts(v);
            current_min_corner(v)=std::min(current_min_corner(v),p.X(v)-dd);
            current_max_corner(v)=std::max(current_max_corner(v),p.X(v)+dd);}}

    for(int v=0;v<3;++v){
        bbox.min_corner(v)=min_corner_per_thread(0)(v);
        bbox.max_corner(v)=max_corner_per_thread(0)(v);
    }

    for(int tid=1;tid<threads;++tid) for(int v=0;v<3;++v){
        bbox.min_corner(v)=std::min(bbox.min_corner(v),min_corner_per_thread(tid)(v));
        bbox.max_corner(v)=std::max(bbox.max_corner(v),max_corner_per_thread(tid)(v));}
    
    for(int v=0;v<3;++v){ 
        bbox.min_corner(v)=std::max((T)0.,bbox.min_corner(v));
        bbox.max_corner(v)=std::min((T)1.,bbox.max_corner(v));}
}
//######################################################################
// Rasterize_Voxels: set inside cell
//######################################################################
template<class T> void MPM_Example<T,2>::
Rasterize_Voxels()
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    using Cell_Iterator             = Grid_Iterator_Cell<T,2>;
    using Hierarchy_Initializer     = Grid_Hierarchy_Initializer<Struct_type,T,2>;
    const Grid<T,2>& grid=hierarchy->Lattice(0);
    Range<int,2> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));
    for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next())
        hierarchy->Activate_Cell(0,iterator.Cell_Index(),Cell_Type_Interior);
    
    hierarchy->Update_Block_Offsets();
    hierarchy->Initialize_Red_Black_Partition(2*threads);

    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
	ras_vx_cnt++;
    ras_vx_rt+=dur.count();
}
//######################################################################
// Rasterize_Voxels: set inside cell
//######################################################################
template<class T> void MPM_Example<T,3>::
Rasterize_Voxels()
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    using Cell_Iterator             = Grid_Iterator_Cell<T,3>;
    using Hierarchy_Initializer     = Grid_Hierarchy_Initializer<Struct_type,T,3>;
    const Grid<T,3>& grid=hierarchy->Lattice(0);
    Range<int,3> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));
    for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next())
        hierarchy->Activate_Cell(0,iterator.Cell_Index(),Cell_Type_Interior);
    
    hierarchy->Update_Block_Offsets();
    hierarchy->Initialize_Red_Black_Partition(2*threads);

    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
	ras_vx_cnt++;
    ras_vx_rt+=dur.count();
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T> void MPM_Example<T,2>::
Initialize_SPGrid()
{
    Compute_Bounding_Box(bbox);
    if(hierarchy!=nullptr) delete hierarchy;
    hierarchy=new Hierarchy(counts,domain,levels);
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T> void MPM_Example<T,3>::
Initialize_SPGrid()
{
    Compute_Bounding_Box(bbox);
    if(hierarchy!=nullptr) delete hierarchy;
    hierarchy=new Hierarchy(counts,domain,levels);
}
//######################################################################
// Populated_Simulated_Particles
//######################################################################
template<class T> void MPM_Example<T,2>::
Populate_Simulated_Particles()
{
    simulated_particles.Clear();
    for(int i=0;i<particles.size();++i)
        if(particles(i).valid)
            simulated_particles.Append(i);
    // Log::cout<<"\nSimulated particles: "<<simulated_particles.size()<<std::endl;
}
//######################################################################
// Populated_Simulated_Particles
//######################################################################
template<class T> void MPM_Example<T,3>::
Populate_Simulated_Particles()
{
    simulated_particles.Clear();
    for(int i=0;i<particles.size();++i)
        if(particles(i).valid)
            simulated_particles.Append(i);
    // Log::cout<<"\nSimulated particles: "<<simulated_particles.size()<<std::endl;
}
//######################################################################
// Max_Particle_Velocity
//######################################################################
template<class T> T MPM_Example<T,2>::
Max_Particle_Velocity() const
{
    Array<T> result_per_thread(threads);
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int tid=omp_get_thread_num();
        T& r=result_per_thread(tid); const int id=simulated_particles(i);
        r=std::max(r,particles(id).V.Norm_Squared());}
    T result=(T)0.;
    for(int tid=0;tid<threads;++tid) result=std::max(result,result_per_thread(tid));
    Log::cout<<"max v: "<<std::sqrt(result)<<std::endl;
    return std::sqrt(result);
}
//######################################################################
// Max_Particle_Velocity
//######################################################################
template<class T> T MPM_Example<T,3>::
Max_Particle_Velocity() const
{
    Array<T> result_per_thread(threads);
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int tid=omp_get_thread_num();
        T& r=result_per_thread(tid); const int id=simulated_particles(i);
        r=std::max(r,particles(id).V.Norm_Squared());}
    T result=(T)0.;
    for(int tid=0;tid<threads;++tid) result=std::max(result,result_per_thread(tid));
    Log::cout<<"max v: "<<std::sqrt(result)<<std::endl;
    return std::sqrt(result);
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T> void MPM_Example<T,2>::
Limit_Dt(T& dt,const T time)
{
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T> void MPM_Example<T,3>::
Limit_Dt(T& dt,const T time)
{
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T> void MPM_Example<T,2>::
Update_Particle_Weights()
{
    const Grid<T,2>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
        for(int i=0;i<simulated_particles.size();i++){
            const int id=simulated_particles(i);
            T_Particle& p=particles(id);
            p.Update_Weights(grid);}
}
//######################################################################
// Update_Particle_Weights
//######################################################################
template<class T> void MPM_Example<T,3>::
Update_Particle_Weights()
{
    const Grid<T,3>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
        for(int i=0;i<simulated_particles.size();i++){
            const int id=simulated_particles(i);
            T_Particle& p=particles(id);
            p.Update_Weights(grid);}
}
//######################################################################
// Group_Particles
//######################################################################
template<class T> void MPM_Example<T,2>::
Group_Particles()
{
    const Grid<T,2>& grid=hierarchy->Lattice(0);
    T_INDEX min_corner=grid.Cell_Indices().min_corner, max_corner=grid.Cell_Indices().max_corner;
    x_intervals.resize(threads); x_intervals(0).min_corner=min_corner(0); x_intervals(threads-1).max_corner=max_corner(0);
    const T ratio=(max_corner(0)-min_corner(0)+1)/threads;
    for(int i=0;i<threads-1;++i){
        int n=min_corner(0)+(max_corner(0)-min_corner(0)+1)*(i+1)/threads;
        x_intervals(i).max_corner=n-1;
        x_intervals(i+1).min_corner=n;}
    for(int i=0;i<threads;++i) for(int j=0;j<threads;++j) particle_bins(i,j).Clear();
#pragma omp parallel for
    for(int i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i);
        const T_INDEX& closest_cell=particles(id).closest_cell;
        const int tid_collect=omp_get_thread_num();
        const Interval<int> particle_x_interval=Interval<int>(closest_cell(0)-1,closest_cell(0)+1);
        for(int tid_process=0;tid_process<threads;++tid_process){
            const Interval<int> thread_x_interval=x_intervals(tid_process);
            if(particle_x_interval.Intersection(thread_x_interval)) particle_bins(tid_process,tid_collect).Append(id);}}    

    // for(int tid_process=0;tid_process<threads;++tid_process){
    //     Log::cout<<"interval: "<<x_intervals(tid_process).min_corner<<","<<x_intervals(tid_process).max_corner<<std::endl;
    //     for(int tid_collect=0;tid_collect<threads;++tid_collect){
    //         const Array<int> index=particle_bins(tid_process,tid_collect);
    //         Log::cout<<"size: "<<index.size()<<std::endl;
    //         // for(int i=0;i<index.size();++i) Log::cout<<i<<": "<<index(i)<<std::endl;
    //     }}
}
//######################################################################
// Group_Particles
//######################################################################
template<class T> void MPM_Example<T,3>::
Group_Particles()
{
    const Grid<T,3>& grid=hierarchy->Lattice(0);
    T_INDEX min_corner=grid.Cell_Indices().min_corner, max_corner=grid.Cell_Indices().max_corner;
    x_intervals.resize(threads); x_intervals(0).min_corner=min_corner(0); x_intervals(threads-1).max_corner=max_corner(0);
    const T ratio=(max_corner(0)-min_corner(0)+1)/threads;
    for(int i=0;i<threads-1;++i){
        int n=min_corner(0)+(max_corner(0)-min_corner(0)+1)*(i+1)/threads;
        x_intervals(i).max_corner=n-1;
        x_intervals(i+1).min_corner=n;}
    for(int i=0;i<threads;++i) for(int j=0;j<threads;++j) particle_bins(i,j).Clear();
#pragma omp parallel for
    for(int i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i);
        const T_INDEX& closest_cell=particles(id).closest_cell;
        const int tid_collect=omp_get_thread_num();
        const Interval<int> particle_x_interval=Interval<int>(closest_cell(0)-1,closest_cell(0)+1);
        for(int tid_process=0;tid_process<threads;++tid_process){
            const Interval<int> thread_x_interval=x_intervals(tid_process);
            if(particle_x_interval.Intersection(thread_x_interval)) particle_bins(tid_process,tid_collect).Append(id);}}    
}
//######################################################################
// Rasterize
//######################################################################
template<class T> void MPM_Example<T,2>::
Rasterize()
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    auto flags=hierarchy->Channel(0,flags_channel); auto mass=hierarchy->Channel(0,mass_channel);       
    auto v0=hierarchy->Channel(0,velocity_channels(0)); auto v1=hierarchy->Channel(0,velocity_channels(1));
    const Grid<T,2>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(int tid_process=0;tid_process<threads;++tid_process){
        const Interval<int> thread_x_interval=x_intervals(tid_process);
        for(int tid_collect=0;tid_collect<threads;++tid_collect){
            Array<int>& index=particle_bins(tid_process,tid_collect);
            for(int i=0;i<index.size();++i){
                T_Particle& p=particles(index(i));T_INDEX& closest_cell=p.closest_cell;
                const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_cell(0),thread_x_interval.max_corner-closest_cell(0));
                for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){
                    TV V=p.V; T p_mass=p.mass; T weight=iterator.Weight(); auto data=iterator.Current_Cell()._data;
                    mass(data)+=weight*p_mass; v0(data)+=weight*(p_mass*V(0)); v1(data)+=weight*(p_mass*V(1));}}}}

    // set flags
    for(int level=0;level<levels;++level) Flag_Setup_Helper<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level));     
    // normalize weights for velocity (to conserve momentum)
    for(int level=0;level<levels;++level) Velocity_Normalization_Helper<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_channels);     
    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    ras_cnt++;
    ras_rt+=dur.count();
}
//######################################################################
// Rasterize
//######################################################################
template<class T> void MPM_Example<T,3>::
Rasterize()
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    auto flags=hierarchy->Channel(0,flags_channel); auto mass=hierarchy->Channel(0,mass_channel);       
    auto v0=hierarchy->Channel(0,velocity_channels(0)); auto v1=hierarchy->Channel(0,velocity_channels(1)); auto v2=hierarchy->Channel(0,velocity_channels(2)); 
    const Grid<T,3>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(int tid_process=0;tid_process<threads;++tid_process){
        const Interval<int> thread_x_interval=x_intervals(tid_process);
        for(int tid_collect=0;tid_collect<threads;++tid_collect){
            Array<int>& index=particle_bins(tid_process,tid_collect);
            for(int i=0;i<index.size();++i){
                T_Particle& p=particles(index(i));T_INDEX& closest_cell=p.closest_cell;
                const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_cell(0),thread_x_interval.max_corner-closest_cell(0));
                for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){
                    TV V=p.V; T p_mass=p.mass; T weight=iterator.Weight(); auto data=iterator.Current_Cell()._data;
                    mass(data)+=weight*p_mass; v0(data)+=weight*(p_mass*V(0)); v1(data)+=weight*(p_mass*V(1)); v2(data)+=weight*(p_mass*V(2));}}}}

    // set flags
    for(int level=0;level<levels;++level) Flag_Setup_Helper<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level));     
    // normalize weights for velocity (to conserve momentum)
    for(int level=0;level<levels;++level) Velocity_Normalization_Helper<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_channels);     
    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    ras_cnt++;
    ras_rt+=dur.count();
}
//######################################################################
// Update_Constitutive_Model_State
//######################################################################
template<class T> void MPM_Example<T,2>::
Update_Constitutive_Model_State()
{
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i); 
        T_Particle &particle=particles(id);    
        particle.constitutive_model.Precompute();}
}
//######################################################################
// Update_Constitutive_Model_State
//######################################################################
template<class T> void MPM_Example<T,3>::
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
template<class T> void MPM_Example<T,2>::
Update_Particle_Velocities_And_Positions(const T dt)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    Array<Array<int> > remove_indices(threads);
    auto vs0=hierarchy->Channel(0,velocity_star_channels(0));   auto vs1=hierarchy->Channel(0,velocity_star_channels(1));
    auto v0=hierarchy->Channel(0,velocity_channels(0));         auto v1=hierarchy->Channel(0,velocity_channels(1));
    Apply_Force(dt);
    const Grid<T,2>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i); 
        T_Particle &p=particles(id);
        TV V_pic=TV(),V_flip=p.V; 
        T_Mat grad_Vp=T_Mat();
        for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){
            auto data=iterator.Current_Cell()._data; T weight=iterator.Weight();
            if(weight>(T)0.){
                TV V_grid({vs0(data),vs1(data)}),delta_V_grid({vs0(data)-v0(data),vs1(data)-v1(data)});
                V_pic+=weight*V_grid; 
                V_flip+=weight*delta_V_grid;
                grad_Vp+=T_Mat::Outer_Product(V_grid,iterator.Weight_Gradient());}}
            p.constitutive_model.Fe+=dt*grad_Vp*p.constitutive_model.Fe;
            p.V=V_flip*flip+V_pic*((T)1.-flip);
            p.X+=V_pic*dt;
        if(!grid.domain.Inside(p.X)){
            remove_indices(omp_get_thread_num()).Append(i);
            p.valid=false;}}

    for(int i=1;i<remove_indices.size();++i)
            remove_indices(0).Append_Elements(remove_indices(i));
        Array<int>::Sort(remove_indices(0));
        for(int i=remove_indices(0).size()-1;i>=0;i--){
            int k=remove_indices(0)(i);
            invalid_particles.Append(simulated_particles(k));
            simulated_particles.Remove_Index(k);}

    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    update_x_v_cnt++;
    update_x_v_rt+=dur.count();
}
//######################################################################
// Update_Particle_Velocities_And_Positions
//######################################################################
template<class T> void MPM_Example<T,3>::
Update_Particle_Velocities_And_Positions(const T dt)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    Array<Array<int> > remove_indices(threads);
    auto vs0=hierarchy->Channel(0,velocity_star_channels(0));   auto vs1=hierarchy->Channel(0,velocity_star_channels(1)); auto vs2=hierarchy->Channel(0,velocity_star_channels(2)); 
    auto v0=hierarchy->Channel(0,velocity_channels(0));         auto v1=hierarchy->Channel(0,velocity_channels(1)); auto v2=hierarchy->Channel(0,velocity_channels(2));         
    Apply_Force(dt);
    const Grid<T,3>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i); 
        T_Particle &p=particles(id);
        TV V_pic=TV(),V_flip=p.V; 
        T_Mat grad_Vp=T_Mat();
        for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){
            auto data=iterator.Current_Cell()._data; T weight=iterator.Weight();
            if(weight>(T)0.){
                TV V_grid({vs0(data),vs1(data),vs2(data)}),delta_V_grid({vs0(data)-v0(data),vs1(data)-v1(data),vs2(data)-v2(data)});
                V_pic+=weight*V_grid; 
                V_flip+=weight*delta_V_grid;
                grad_Vp+=T_Mat::Outer_Product(V_grid,iterator.Weight_Gradient());}}
            p.constitutive_model.Fe+=dt*grad_Vp*p.constitutive_model.Fe;
            p.V=V_flip*flip+V_pic*((T)1.-flip);
            p.X+=V_pic*dt;
        if(!grid.domain.Inside(p.X)){
            remove_indices(omp_get_thread_num()).Append(i);
            p.valid=false;}}

    for(int i=1;i<remove_indices.size();++i)
            remove_indices(0).Append_Elements(remove_indices(i));
        Array<int>::Sort(remove_indices(0));
        for(int i=remove_indices(0).size()-1;i>=0;i--){
            int k=remove_indices(0)(i);
            invalid_particles.Append(simulated_particles(k));
            simulated_particles.Remove_Index(k);}

    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    update_x_v_cnt++;
    update_x_v_rt+=dur.count();
}
//######################################################################
// Apply_Force
//######################################################################
template<class T> void MPM_Example<T,2>::
Apply_Force(const T dt)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    Apply_Explicit_Force(dt);
    Grid_Based_Collision(true);
    if(true){
    Conjugate_Gradient<T> cg;
    Krylov_Solver<T>* solver=(Krylov_Solver<T>*)&cg;
    MPM_CG_System<Struct_type,T,2> mpm_system(*hierarchy,simulated_particles,particles,particle_bins,x_intervals,barriers,(T)0.,dt,threads);
    mpm_system.use_preconditioner=false;
    // set rhs here
    for(int level=0;level<levels;++level) MPM_RHS_Helper<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_star_channels,rhs_channels);     
    // clear channels
    Reset_Solver_Channels();
    MPM_CG_Vector<Struct_type,T,2> solver_vp(*hierarchy,velocity_star_channels),solver_rhs(*hierarchy,rhs_channels),solver_q(*hierarchy,q_channels),solver_s(*hierarchy,s_channels),solver_r(*hierarchy,r_channels),solver_k(*hierarchy,z_channels),solver_z(*hierarchy,z_channels);
    solver->Solve(mpm_system,solver_vp,solver_rhs,solver_q,solver_s,solver_r,solver_k,solver_z,solver_tolerance,0,solver_iterations);}
    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    apply_force_cnt++;
    apply_force_rt+=dur.count();
}
//######################################################################
// Apply_Force
//######################################################################
template<class T> void MPM_Example<T,3>::
Apply_Force(const T dt)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    Apply_Explicit_Force(dt);
    Grid_Based_Collision(true);
    if(true){
    Conjugate_Gradient<T> cg;
    Krylov_Solver<T>* solver=(Krylov_Solver<T>*)&cg;
    MPM_CG_System<Struct_type,T,3> mpm_system(*hierarchy,simulated_particles,particles,particle_bins,x_intervals,barriers,(T)0.,dt,threads);
    mpm_system.use_preconditioner=false;
    // set rhs here
    for(int level=0;level<levels;++level) MPM_RHS_Helper<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_star_channels,rhs_channels);     
    // clear channels
    Reset_Solver_Channels();
    MPM_CG_Vector<Struct_type,T,3> solver_vp(*hierarchy,velocity_star_channels),solver_rhs(*hierarchy,rhs_channels),solver_q(*hierarchy,q_channels),solver_s(*hierarchy,s_channels),solver_r(*hierarchy,r_channels),solver_k(*hierarchy,z_channels),solver_z(*hierarchy,z_channels);
    solver->Solve(mpm_system,solver_vp,solver_rhs,solver_q,solver_s,solver_r,solver_k,solver_z,solver_tolerance,0,solver_iterations);}
    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    apply_force_cnt++;
    apply_force_rt+=dur.count();
}
//######################################################################
// Apply_Explicit_Force
//######################################################################
template<class T> void MPM_Example<T,2>::
Apply_Explicit_Force(const T dt)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    const Grid<T,2>& grid=hierarchy->Lattice(0);
    auto f0=hierarchy->Channel(0,f_channels(0)); auto f1=hierarchy->Channel(0,f_channels(1));
    Array<T> thread_rt(threads);
    Array<int> thread_cnt(threads);
#pragma omp parallel for
    for(int tid_process=0;tid_process<threads;++tid_process){
        const Interval<int>& thread_x_interval=x_intervals(tid_process);
        T& rt=thread_rt(tid_process);
        int& cnt=thread_cnt(tid_process);
        for(int tid_collect=0;tid_collect<threads;++tid_collect){
            const Array<int>& index=particle_bins(tid_process,tid_collect);
            for(int i=0;i<index.size();++i){
                T_Particle& p=particles(index(i));T_INDEX& closest_cell=p.closest_cell;
                T_Mat P=p.constitutive_model.P(),F=p.constitutive_model.Fe;T V0=p.volume;
                T_Mat V0_P_FT=P.Times_Transpose(F)*V0;                    
                const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_cell(0),thread_x_interval.max_corner-closest_cell(0));  
                for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){
                    // high_resolution_clock::time_point tbt=high_resolution_clock::now();
                    auto data=iterator.Current_Cell()._data;         
                    TV body_force=gravity*p.mass*iterator.Weight(); TV inner_force=V0_P_FT*iterator.Weight_Gradient();
                    f0(data)-=inner_force(0); f1(data)-=inner_force(1);
                    f0(data)+=body_force(0); f1(data)+=body_force(1);
                }}}}
    
    // T total_rt=0;
    // int total_cnt=0;
    // for(int i=0;i<threads;++i){ total_cnt+=thread_cnt(i); total_rt+=thread_rt(i);}
    // Log::cout<<"per particle: "<<total_rt/total_cnt<<std::endl;

    for(int level=0;level<levels;++level) Explicit_Force_Helper<Struct_type,T,2>(hierarchy->Allocator(level),hierarchy->Blocks(level),f_channels,velocity_channels,velocity_star_channels,dt);
    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    explicit_force_cnt++;
    explicit_force_rt+=dur.count();
}
//######################################################################
// Apply_Explicit_Force
//######################################################################
template<class T> void MPM_Example<T,3>::
Apply_Explicit_Force(const T dt)
{
    high_resolution_clock::time_point tb=high_resolution_clock::now();
    const Grid<T,3>& grid=hierarchy->Lattice(0);
    auto f0=hierarchy->Channel(0,f_channels(0)); auto f1=hierarchy->Channel(0,f_channels(1)); auto f2=hierarchy->Channel(0,f_channels(2)); 
    Array<T> thread_rt(threads);
    Array<int> thread_cnt(threads);
#pragma omp parallel for
    for(int tid_process=0;tid_process<threads;++tid_process){
        const Interval<int>& thread_x_interval=x_intervals(tid_process);
        T& rt=thread_rt(tid_process);
        int& cnt=thread_cnt(tid_process);
        for(int tid_collect=0;tid_collect<threads;++tid_collect){
            const Array<int>& index=particle_bins(tid_process,tid_collect);
            for(int i=0;i<index.size();++i){
                T_Particle& p=particles(index(i));T_INDEX& closest_cell=p.closest_cell;
                T_Mat P=p.constitutive_model.P(),F=p.constitutive_model.Fe;T V0=p.volume;
                T_Mat V0_P_FT=P.Times_Transpose(F)*V0;                    
                const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_cell(0),thread_x_interval.max_corner-closest_cell(0));  
                for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){
                    // high_resolution_clock::time_point tbt=high_resolution_clock::now();
                    auto data=iterator.Current_Cell()._data;         
                    TV body_force=gravity*p.mass*iterator.Weight(); TV inner_force=V0_P_FT*iterator.Weight_Gradient();
                    f0(data)-=inner_force(0); f1(data)-=inner_force(1); f2(data)-=inner_force(2); 
                    f0(data)+=body_force(0); f1(data)+=body_force(1); f2(data)+=body_force(2);}}}}
    
    // T total_rt=0;
    // int total_cnt=0;
    // for(int i=0;i<threads;++i){ total_cnt+=thread_cnt(i); total_rt+=thread_rt(i);}
    // Log::cout<<"per particle: "<<total_rt/total_cnt<<std::endl;

    for(int level=0;level<levels;++level) Explicit_Force_Helper<Struct_type,T,3>(hierarchy->Allocator(level),hierarchy->Blocks(level),f_channels,velocity_channels,velocity_star_channels,dt);
    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
    explicit_force_cnt++;
    explicit_force_rt+=dur.count();
}
//######################################################################
// Grid_Based_Collision
//######################################################################
template<class T> void MPM_Example<T,2>::
Grid_Based_Collision(const bool detect_collision)
{   
    if(detect_collision) for(int level=0;level<levels;++level) Grid_Based_Collision_Helper<Struct_type,T,2>(hierarchy->Lattice(level),hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_star_channels,barriers);
}
//######################################################################
// Grid_Based_Collision
//######################################################################
template<class T> void MPM_Example<T,3>::
Grid_Based_Collision(const bool detect_collision)
{   
    if(detect_collision) for(int level=0;level<levels;++level) Grid_Based_Collision_Helper<Struct_type,T,3>(hierarchy->Lattice(level),hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_star_channels,barriers);
}
//######################################################################
// Estimate_Particle_Volumes
//######################################################################
template<class T> void MPM_Example<T,2>::
Estimate_Particle_Volumes()
{   
    auto mass=hierarchy->Channel(0,mass_channel);
    const Grid<T,2>& grid=hierarchy->Lattice(0); const T one_over_volume_per_cell=(T)1./grid.dX.Product();
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int id=simulated_particles(i); T_Particle& p=particles(id);     
        T particle_density=(T)0.;
        for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){
            particle_density+=iterator.Weight()*mass(iterator.Current_Cell()._data);}
        particle_density*=one_over_volume_per_cell;
        p.volume=p.mass/particle_density;} 
}
//######################################################################
// Estimate_Particle_Volumes
//######################################################################
template<class T> void MPM_Example<T,3>::
Estimate_Particle_Volumes()
{   
    auto mass=hierarchy->Channel(0,mass_channel);
    const Grid<T,3>& grid=hierarchy->Lattice(0); const T one_over_volume_per_cell=(T)1./grid.dX.Product();
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int id=simulated_particles(i); T_Particle& p=particles(id);     
        T particle_density=(T)0.;
        for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){
            particle_density+=iterator.Weight()*mass(iterator.Current_Cell()._data);}
        particle_density*=one_over_volume_per_cell;
        p.volume=p.mass/particle_density;} 
}
//######################################################################
// Register_Options
//######################################################################
template<class T> void MPM_Example<T,2>::
Register_Options()
{
    Base::Register_Options();
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Double_Argument("-cfl",(T)0.1,"CFL number.");
    parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(32.),"n","Grid resolution");
}
//######################################################################
// Register_Options
//######################################################################
template<class T> void MPM_Example<T,3>::
Register_Options()
{
    Base::Register_Options();
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Double_Argument("-cfl",(T)0.1,"CFL number.");
    parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(32.),"n","Grid resolution");
}
//######################################################################
// Test
//######################################################################
template<class T> void MPM_Example<T,2>::
Test()
{
    // const Grid<T,2>& grid=hierarchy->Lattice(0); 
    // Log::cout<<"node: "<<grid.Node(T_INDEX(33))<<std::endl;
    // Log::cout<<"node: "<<grid.Node(T_INDEX(34))<<std::endl;
    // Log::cout<<"min corner: "<<grid.domain.min_corner<<std::endl;
    // Log::cout<<"max corner: "<<grid.domain.max_corner<<std::endl;
    // const int iterations=1e8;
    // auto mass=hierarchy->Channel(0,mass_channel);
    // T_INDEX index({0,0}); auto data=index._data;
    // // read
    // high_resolution_clock::time_point tb1=high_resolution_clock::now();
    // for(int i=0;i<iterations;i++) mass(data);
    // high_resolution_clock::time_point te1=high_resolution_clock::now();
	// duration<double> dur1=duration_cast<duration<double>>(te1-tb1);
    // Log::cout<<"read from channel: "<<dur1.count()/iterations<<std::endl;

    // // write
    // high_resolution_clock::time_point tb2=high_resolution_clock::now();
    // for(int i=0;i<iterations;i++) mass(data)=(T)0.;
    // high_resolution_clock::time_point te2=high_resolution_clock::now();
	// duration<double> dur2=duration_cast<duration<double>>(te2-tb2);
    // Log::cout<<"write to channel: "<<dur2.count()/iterations<<std::endl;

    
    // // influence iterator
    // high_resolution_clock::time_point tb3=high_resolution_clock::now();
    // for(int it=0;it<iterations;++it){
    //     for(unsigned i=0;i<particles.size();++i){T_Particle& p=particles(i);     
    //         for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){}}}
    // high_resolution_clock::time_point te3=high_resolution_clock::now();
	// duration<double> dur3=duration_cast<duration<double>>(te3-tb3);
    // Log::cout<<"iterator: "<<dur3.count()/iterations<<std::endl;




}
//######################################################################
// Test
//######################################################################
template<class T> void MPM_Example<T,3>::
Test()
{
    const int iterations=1e8;
    auto mass=hierarchy->Channel(0,mass_channel);
    T_INDEX index({0,0,0}); auto data=index._data;
    // read
    high_resolution_clock::time_point tb1=high_resolution_clock::now();
    for(int i=0;i<iterations;i++) mass(data);
    high_resolution_clock::time_point te1=high_resolution_clock::now();
	duration<double> dur1=duration_cast<duration<double>>(te1-tb1);
    Log::cout<<"read from channel: "<<dur1.count()/iterations<<std::endl;

    // write
    high_resolution_clock::time_point tb2=high_resolution_clock::now();
    for(int i=0;i<iterations;i++) mass(data)=(T)0.;
    high_resolution_clock::time_point te2=high_resolution_clock::now();
	duration<double> dur2=duration_cast<duration<double>>(te2-tb2);
    Log::cout<<"write to channel: "<<dur2.count()/iterations<<std::endl;

    
    // influence iterator
    high_resolution_clock::time_point tb3=high_resolution_clock::now();
    for(int it=0;it<iterations;++it){
        for(unsigned i=0;i<particles.size();++i){T_Particle& p=particles(i);     
            for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){}}}
    high_resolution_clock::time_point te3=high_resolution_clock::now();
	duration<double> dur3=duration_cast<duration<double>>(te3-tb3);
    Log::cout<<"iterator: "<<dur3.count()/iterations<<std::endl;


}
//######################################################################
// Parse_Options
//######################################################################
template<class T> void MPM_Example<T,2>::
Parse_Options()
{
    Base::Parse_Options();
    threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(threads);
    Base::test_number=parse_args->Get_Integer_Value("-test_number");
    cfl=parse_args->Get_Double_Value("-cfl");
    levels=parse_args->Get_Integer_Value("-levels");
    auto cell_counts_2d=parse_args->Get_Vector_2D_Value("-size");
    for(int v=0;v<2;++v) counts(v)=cell_counts_2d(v);
}
//######################################################################
// Parse_Options
//######################################################################
template<class T> void MPM_Example<T,3>::
Parse_Options()
{
    Base::Parse_Options();
    threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(threads);
    Base::test_number=parse_args->Get_Integer_Value("-test_number");
    cfl=parse_args->Get_Double_Value("-cfl");
    levels=parse_args->Get_Integer_Value("-levels");
    auto cell_counts_3d=parse_args->Get_Vector_3D_Value("-size");for(int v=0;v<3;++v) counts(v)=cell_counts_3d(v);
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T> void MPM_Example<T,2>::
Write_Output_Files(const int frame) const
{
    if(frame==first_frame){std::string deformables_filename=output_directory+"/common/metadata.mpm2d";
        File_Utilities::Write_To_Text_File(deformables_filename,std::to_string(frame));}

    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/frame_title",frame_title);

    File_Utilities::Write_To_File(output_directory+"/"+std::to_string(frame)+"/particles",particles);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T> void MPM_Example<T,3>::
Write_Output_Files(const int frame) const
{
    if(frame==first_frame){std::string deformables_filename=output_directory+"/common/metadata.mpm3d";
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
template<class T> void MPM_Example<T,2>::
Read_Output_Files(const int frame)
{
    File_Utilities::Read_From_File(output_directory+"/"+std::to_string(frame)+"/particles",particles);
}
//######################################################################
// Read_Output_Files
//######################################################################
template<class T> void MPM_Example<T,3>::
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

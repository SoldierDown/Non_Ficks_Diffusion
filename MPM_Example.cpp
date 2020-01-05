//!#####################################################################
//! \file MPM_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Grids/Grid_Iterator_Cell.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <omp.h>
#include "MPM_Example.h"
#include "Velocity_Normalization_Helper.h"
using namespace Nova;
using namespace SPGrid;
//######################################################################
// Constructor
//######################################################################
template<class T,int d> MPM_Example<T,d>::
MPM_Example()
    :Base(),hierarchy(nullptr)
{
    gravity=-TV::Axis_Vector(1)*(T)2.;
    flip=(T)0.9;
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
    valid_nodes                             = &Struct_type::ch10;
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
    valid_grid_indices_thread.resize(threads);
    Initialize_Particles();
    Initialize_SPGrid();
}
//######################################################################
// Reset_Grid_Based_Variables
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Reset_Grid_Based_Variables()
{
    // clear mass, velocity and force channels
    for(int level=0;level<levels;++level){
        Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),mass_channel);
        Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),valid_nodes);
        for(int v=0;v<d;++v) {
            Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_channels(v));
            Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_star_channels(v));
            Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),f_channels(v));}}
    valid_grid_indices.Clear();
    for(int i=0;i<valid_grid_indices_thread.size();++i) valid_grid_indices_thread(i).Clear();
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
            T dd=(T)2./counts(v);
            current_min_corner(v)=std::min(current_min_corner(v),p.X(v)-dd);
            current_max_corner(v)=std::max(current_max_corner(v),p.X(v)+dd);}}

    for(int v=0;v<d;++v){
        bbox.min_corner(v)=min_corner_per_thread(0)(v);
        bbox.max_corner(v)=max_corner_per_thread(0)(v);
    }

    for(int tid=1;tid<threads;++tid) for(int v=0;v<d;++v){
        bbox.min_corner(v)=std::min(bbox.min_corner(v),min_corner_per_thread(tid)(v));
        bbox.max_corner(v)=std::max(bbox.max_corner(v),max_corner_per_thread(tid)(v));}
    

    printf("\nafter: min corner: %f, %f, max corner: %f, %f\n",bbox.min_corner(0),bbox.min_corner(1),bbox.max_corner(0),bbox.max_corner(1));
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

    printf("\n before min corner: %f, %f, max corner: %f, %f\n",bbox.min_corner(0),bbox.min_corner(1),bbox.max_corner(0),bbox.max_corner(1));
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
    std::cout<<"simulated particles: "<<simulated_particles.size()<<std::endl;
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
// Rasterize
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Rasterize()
{
    gravity=Vector<T,d>::Axis_Vector(1)*(T)-2.;
    const Grid<T,d>& grid=hierarchy->Lattice(0);
    std::cout<<"cell width:"<<grid.dX(0)<<std::endl;
    for(unsigned i=0;i<simulated_particles.size();++i){ T weight_sum=(T)0.; const int id=simulated_particles(i); T_Particle& p=particles(id); T_INDEX closest_node=grid.Closest_Node(p.X);
        for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){T_INDEX current_node=closest_node+iterator.Index();
            if(grid.Node_Indices().Inside(current_node)){const TV current_node_location=grid.Node(current_node);T weight=N(p.X-current_node_location);
                weight_sum+=weight;
                hierarchy->Channel(0,mass_channel)(current_node._data)+=weight*p.mass;
                if(weight>(T)0.) hierarchy->Channel(0,valid_nodes)(current_node._data)=(T)1.;
                for(int v=0;v<d;++v) {hierarchy->Channel(0,velocity_channels(v))(current_node._data)+=weight*p.mass*p.V(v);
                    // hierarchy->Channel(0,f_channels(v))(current_node._data)+=weight*p.mass*gravity(v);
                }}}

            //std::cout<<"weight sum: "<<weight_sum<<std::endl;    
        }
    
    // normalize weights for velocity (to conserve momentum)
    for(int level=0;level<levels;++level)
        Velocity_Normalization_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),velocity_channels,mass_channel);

    T max_mass=(T)-100.;
    T min_mass=FLT_MAX;
    T max_v=(T)-100.;
    T min_v=FLT_MAX;
    using Cell_Iterator     = Grid_Iterator_Cell<T,d>;
    Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));

    for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next()){
        T_INDEX current_node=iterator.Cell_Index(); T mass=hierarchy->Channel(0,mass_channel)(current_node._data);
        if(mass>max_mass) max_mass=mass;
        if(mass<min_mass) min_mass=mass;
        for(int v=0;v<d;++v){
            T vnorm=abs(hierarchy->Channel(0,velocity_channels(v))(current_node._data));
            if(vnorm>max_v) max_v=vnorm;
            if(vnorm<min_v) min_v=vnorm;
        }
    }
    printf("mass range: %f, %f\n", min_mass, max_mass);
    printf("v range: %f, %f\n", min_v, max_v);

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

    Array<Array<int> > remove_indices(threads);
    const Grid<T,d>& grid=hierarchy->Lattice(0);
//#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){
        const int id=simulated_particles(i); 
        T_Particle &particle=particles(id);
        particle.V=TV();
        T_INDEX closest_node=grid.Closest_Node(particle.X); 
        TV V_pic=TV(),V_flip=particle.V; 
        Matrix<T,d> grad_Vp=Matrix<T,d>();
        for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){
            T_INDEX current_node=closest_node+iterator.Index();
            
            if(grid.Node_Indices().Inside(current_node)){
                const TV current_node_location=grid.Node(current_node);
                T weight=N(particle.X-current_node_location);
                for(int v=0;v<d;v++) particle.V(v)+=weight*hierarchy->Channel(0,velocity_star_channels(v))(current_node._data);
                if(weight>(T)0.){ 
                    TV weight_grad=dN(particle.X-current_node_location),V_grid, delta_V_grid;
                    for(int v=0;v<d;++v) { V_grid(v)=hierarchy->Channel(0,velocity_star_channels(v))(current_node._data);
                        delta_V_grid(v)=hierarchy->Channel(0,velocity_star_channels(v))(current_node._data)-hierarchy->Channel(0,velocity_channels(v))(current_node._data);}
                    //printf("weight: %f, V_grid: %f\n", weight, V_grid.Norm())      ;
                    V_pic+=weight*V_grid; 
                    V_flip+=weight*delta_V_grid;
                    grad_Vp+=Matrix<T,d>::Outer_Product(V_grid,weight_grad);}}}

            particle.constitutive_model.Fe+=dt*grad_Vp*particle.constitutive_model.Fe;
            //printf("V_flip: %f, V_pic: %f\n", V_flip.Norm(), V_pic.Norm());
            particle.V=V_flip*flip+V_pic*(1-flip);
            
            particle.X+=V_pic*dt;

        if(!grid.domain.Inside(particle.X)) particle.valid=false;

    }
    // for(int i=1;i<remove_indices.size();++i)
    //     remove_indices(0).Append_Elements(remove_indices(i));
    // Array<int>::Sort(remove_indices(1));
    // for(int i=remove_indices(0).size()-1;i>=0;--i){
    //     int k=remove_indices(0)(i);
    //     invalid_particles.Append(simulated_particles(k));
    //     simulated_particles.Remove_Index(k);}

}
//######################################################################
// Apply_Force
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Apply_Force(const T dt)
{
    Apply_Explicit_Force(dt);
    using Cell_Iterator             = Grid_Iterator_Cell<T,d>;

    const Grid<T,d>& grid=hierarchy->Lattice(0);
    Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));

    for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next()){
        T_INDEX current_node=iterator.Cell_Index(); T mass=hierarchy->Channel(0,mass_channel)(current_node._data);
        if(hierarchy->Channel(0,valid_nodes)(current_node._data)>(T)0.){
            for(int v=0;v<d;++v) { hierarchy->Channel(0,velocity_star_channels(v))(current_node._data)=hierarchy->Channel(0,velocity_channels(v))(current_node._data)
                                                                                            +dt/mass*hierarchy->Channel(0,f_channels(v))(current_node._data);
            //printf("v*: %f, v: %f, dt/mass: %f, f: %f\n",hierarchy->Channel(0,velocity_star_channels(v))(current_node._data), hierarchy->Channel(0,velocity_channels(v))(current_node._data), dt/mass, hierarchy->Channel(0,f_channels(v))(current_node._data));
    }}
    }
    Grid_Based_Collison();
}
//######################################################################
// Apply_Explicit_Force
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Apply_Explicit_Force(const T dt)
{
    const Grid<T,d>& grid=hierarchy->Lattice(0);
#pragma omp parallel for
    for(unsigned i=0;i<simulated_particles.size();++i){const int id=simulated_particles(i); T_Particle &particle=particles(id); T V0=particle.volume;
        Matrix<T,d> P=particle.constitutive_model.P(),F=particle.constitutive_model.Fe, V0_P_FT=P.Times_Transpose(F)*V0;
        T_INDEX closest_node=grid.Closest_Node(particle.X);
        for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){T_INDEX current_node=closest_node+iterator.Index();
            if(grid.Node_Indices().Inside(current_node)){const TV current_node_location=grid.Node(current_node);T weight=N(particle.X-current_node_location);
                if(weight>(T)0.){ TV weight_grad=dN(particle.X-current_node_location);
                    for(int v=0;v<d;++v) { hierarchy->Channel(0,f_channels(v))(current_node._data)-=(V0_P_FT*weight_grad)(v);
                        hierarchy->Channel(0,f_channels(v))(current_node._data)+=gravity(v)*particle.mass*weight;
    }}}}}
}
//######################################################################
// Grid_Based_Collision
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Grid_Based_Collison()
{   
    T mu=(T)0.;
    using Cell_Iterator     = Grid_Iterator_Cell<T,d>;
    const Grid<T,d>& grid=hierarchy->Lattice(0);
    Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));

    for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next()){
        T_INDEX current_node=iterator.Cell_Index(); 
        const TV current_node_location=grid.Node(current_node);
        if(current_node_location(1)<(T).1){
            TV normal=TV::Axis_Vector(1)*(T)1.;
            TV vel=TV();
            for(int v=0;v<d;++v) vel(v)=hierarchy->Channel(0,velocity_star_channels(v))(current_node._data);
            T projection=vel.Dot_Product(normal);
            if(projection<(T)0.){
                vel-=normal*projection;
                if(-projection*mu<vel.Norm()) vel+=vel.Normalized()*projection*mu;
                else vel=TV();}
            for(int v=0;v<d;++v) hierarchy->Channel(0,velocity_star_channels(v))(current_node._data);
        }
    } 
}
//######################################################################
// Estimate_Particle_Volumes
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Estimate_Particle_Volumes()
{   
    const Grid<T,d>& grid=hierarchy->Lattice(0);
    T one_over_volume_per_cell=(T)1.;
    for(int v=0;v<d;++v) one_over_volume_per_cell*=grid.dX(v);
    for(unsigned i=0;i<simulated_particles.size();++i){const int id=simulated_particles(i); T_Particle& particle=particles(id); 
        T_INDEX closest_node=grid.Closest_Node(particle.X); T particle_density=(T)0.;
        for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){T_INDEX current_node=closest_node+iterator.Index();
            if(grid.Node_Indices().Inside(current_node)){const TV current_node_location=grid.Node(current_node);T weight=N(particles(id).X-current_node_location);
                particle_density+=weight*hierarchy->Channel(0,mass_channel)(current_node._data);}}
                particle.volume=particle.mass/particle_density;}   
}
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Register_Options()
{
    Base::Register_Options();
    parse_args->Add_Integer_Argument("-threads",8,"Number of threads for OpenMP to use");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL number.");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(32.),"n","Grid resolution");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(32.),"n","Grid resolution");
}
//######################################################################
// Test
//######################################################################
template<class T,int d> void MPM_Example<T,d>::
Test()
{
    using Cell_Iterator             = Grid_Iterator_Cell<T,d>;

    const Grid<T,d>& grid=hierarchy->Lattice(0);
    Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(bbox.min_corner),grid.Clamp_To_Cell(bbox.max_corner));

    T min_mass=(T)100.;
    T max_mass=(T)-100.;
    for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next()){
        T_INDEX current_node=iterator.Cell_Index(); T mass=hierarchy->Channel(0,mass_channel)(current_node._data);
        if(mass>max_mass) max_mass=mass;
        if(mass<min_mass) min_mass=mass;   
    }
    std::cout<<"min mass: "<<min_mass<<", max mass: "<<max_mass<<std::endl;
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

    cfl=(T)parse_args->Get_Double_Value("-cfl");
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

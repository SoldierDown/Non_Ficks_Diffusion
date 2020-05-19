//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Dynamics/Hierarchy/Visualization/Grid_Hierarchy_Visualization.h>
#include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Log/Log.h>
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/Tools/Utilities/Constants.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Pthread_Queue.h>
#include "Poisson_CG_System.h"
#include "Copy_Interior_Flags.h"
#include "Copy_Unsubdivided_Cells.h"
#include "Initialize_Dirichlet_Cells.h"
#include "Initialize_Levelset.h"
#include "Multigrid_Data.h"
#include "../Poisson_Data.h"
#include "../Rasterizers/Adaptive_Sphere_Rasterizer.h"
#include "../Rasterizers/Central_Rasterizer.h"
#include "../Rasterizers/Randomized_Rasterizer.h"
#include "Sphere_Levelset.h"
#include <omp.h>

using namespace Nova;
using namespace SPGrid;

extern Pthread_Queue* pthread_queue;

namespace Nova{
int number_of_threads=0;
}

//template<class Struct_type,class T,int d>
//void Initialize_Dirichlet_Cells(Grid_Hierarchy<Struct_type,T,d>& hierarchy,T Struct_type::* levelset_channel)
//{
//    using Flags_type                            = typename Struct_type::Flags_type;
//    using Allocator_type                        = SPGrid_Allocator<Struct_type,d>;
//    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
//
//    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
//        auto levelset=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(levelset_channel);
//        auto flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
//
//        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
//            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
//                if(flags(offset)&Cell_Type_Interior && levelset(offset)>(T)0.){
//                    flags(offset)|=Cell_Type_Dirichlet;
//                    flags(offset)&=~Cell_Type_Interior;}}}
//}

template<class Struct_type,class T,int d>
void Initialize_Guess(Grid_Hierarchy<Struct_type,T,d>& hierarchy,T Struct_type::* u_channel,const bool random_guess)
{
    using TV                                    = Vector<T,d>;
    using T_INDEX                               = Vector<int,d>;
    using Flags_type                            = typename Struct_type::Flags_type;
    using Allocator_type                        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;

    for(int level=0;level<hierarchy.Levels();++level)
        SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),u_channel);

    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(u_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            std::array<int,d> base_index_s=Flag_array_mask::LinearToCoord(offset);
            T_INDEX base_index=*reinterpret_cast<T_INDEX*>(&base_index_s);

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior){const TV X=hierarchy.Lattice(level).Center(index);
                    data(offset)=(T)sin(2.*pi*X(0))*sin(2.*pi*X(1));}
                range_iterator.Next();}}}
}

template<class Struct_type,class T,int d>
void Compute_Error(Grid_Hierarchy<Struct_type,T,d>& hierarchy,T Struct_type::* channel)
{
    using TV                                    = Vector<T,d>;
    using T_INDEX                               = Vector<int,d>;
    using Flags_type                            = typename Struct_type::Flags_type;
    using Allocator_type                        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;

    T max_error=0.;
    T_INDEX max_index;

    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            std::array<int,d> base_index_s=Flag_array_mask::LinearToCoord(offset);
            T_INDEX base_index=*reinterpret_cast<T_INDEX*>(&base_index_s);

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior){const TV X=hierarchy.Lattice(level).Center(index);T sum=0;
                    for(int axis=0;axis<d;++axis) sum+=Nova_Utilities::Sqr(X(axis));
                    max_error=std::max(max_error,std::fabs(data(offset)-sum+(T).0625));}
                range_iterator.Next();}}}

    Log::cout<<"Max Error: "<<max_error<<std::endl;
}

template<class Struct_type,class T,int d>
void Compute_Right_Hand_Side(Grid_Hierarchy<Struct_type,T,d>& hierarchy,T Struct_type::* channel,const T& value)
{
    using Flags_type                            = typename Struct_type::Flags_type;
    using Allocator_type                        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;

    double scaling_factor=hierarchy.Lattice(0).one_over_dX.Product();

    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        const T cell_volume=hierarchy.Lattice(level).dX.Product();

        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) data(offset)=scaling_factor*value*cell_volume;}}
}

template<class Struct_type,class T,int d>
void Compute_Graded_Octree(Grid_Hierarchy<Struct_type,T,d>*& current_hierarchy)
{

    using TV                                    = Vector<T,d>;
    using T_INDEX                               = Vector<int,d>;
    using Flags_type                            = typename Struct_type::Flags_type;
    using Hierarchy                             = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Lookup                      = Grid_Hierarchy_Lookup<Struct_type,T,d>;
    using Allocator_type                        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper                       = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator                        = SPGrid_Block_Iterator<Flag_array_mask>;

    enum {number_of_nodes_per_cell              = Topology_Helper::number_of_nodes_per_cell};

    uint64_t nodes_of_cell_offsets[number_of_nodes_per_cell];
    Topology_Helper::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
    const int levels=current_hierarchy->Levels();

    for(int level=0;level<levels-2;++level){
        Hierarchy *graded_hierarchy=new Hierarchy(current_hierarchy->Lattice(0),levels);

        // copy cells from all levels <= level+1
        for(int current_level=0;current_level<=level+1;++current_level){
            for(Block_Iterator iterator(current_hierarchy->Blocks(current_level));iterator.Valid();iterator.Next_Block())
                graded_hierarchy->Page_Map(current_level).Set_Page(iterator.Offset());
            Copy_Interior_Flags<Struct_type,T,d>(current_hierarchy->Allocator(current_level),current_hierarchy->Blocks(current_level),*graded_hierarchy,current_level);}

        // balance level and level+1
        {
            auto blocks=current_hierarchy->Blocks(level);
            auto block_size=current_hierarchy->Allocator(level).Block_Size();
            auto flags=current_hierarchy->Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

            for(unsigned block=0;block<blocks.second;++block){uint64_t offset=blocks.first[block];
                Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
                std::array<int,d> base_index_s=Flag_array_mask::LinearToCoord(offset);
                T_INDEX base_index=*reinterpret_cast<T_INDEX*>(&base_index_s);

                for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){const T_INDEX index=base_index+range_iterator.Index();
                    if(flags(offset)&Cell_Type_Ghost){uint64_t new_offset;int new_level;const TV X=current_hierarchy->Lattice(level).Center(index);
                        bool lookup=Hierarchy_Lookup::Cell_Lookup(*current_hierarchy,X,new_offset,new_level);
                        if(lookup && new_level>level+1){T_INDEX new_index(Flag_array_mask::LinearToCoord(new_offset));
                            while(new_level>level+1){uint64_t base_offset=Flag_array_mask::UpsampleOffset(new_offset);new_index=2*new_index-1;int node=0;
                                for(Range_Iterator<d> node_iterator(T_INDEX(),T_INDEX(1));node_iterator.Valid();node_iterator.Next(),++node){
                                    uint64_t node_offset=Flag_array_mask::Packed_Add(base_offset,nodes_of_cell_offsets[node]);T_INDEX node_index=new_index+node_iterator.Index();
                                    if(current_hierarchy->Lattice(new_level-1).Cell_Domain(node_index).Inside(X)){new_index=node_index;new_offset=node_offset;}
                                    else if(!(graded_hierarchy->template Set<unsigned>(new_level-1,&Struct_type::flags).Is_Set(node_offset,Cell_Type_Interior)))
                                        graded_hierarchy->Activate_Cell(new_level-1,node_offset,Cell_Type_Interior);}
                                --new_level;}
                            if(!(graded_hierarchy->template Set<unsigned>(new_level,&Struct_type::flags).Is_Set(new_offset,Cell_Type_Interior)))
                                graded_hierarchy->Activate_Cell(new_level,new_offset,Cell_Type_Interior);}}
                    range_iterator.Next();}}
        }

        // copy all cells from levels >= level+2 that were not subdivided
        for(int current_level=level+2;current_level<levels;++current_level)
            Copy_Unsubdivided_Cells<Struct_type,T,d>(current_hierarchy->Allocator(current_level),current_hierarchy->Blocks(current_level),
                                                     *graded_hierarchy,nodes_of_cell_offsets,current_level);
        delete current_hierarchy;

        // flag ghost cells and update block offsets
        Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Cells(*graded_hierarchy);
        graded_hierarchy->Update_Block_Offsets();

        // create minimal hierarchy that does not contain a cell and all its children
        Hierarchy *minimal_hierarchy=new Hierarchy(graded_hierarchy->Lattice(0),levels);
        for(Block_Iterator iterator(graded_hierarchy->Blocks(0));iterator.Valid();iterator.Next_Block())
            minimal_hierarchy->Page_Map(0).Set_Page(iterator.Offset());
        Copy_Interior_Flags<Struct_type,T,d>(graded_hierarchy->Allocator(0),graded_hierarchy->Blocks(0),*minimal_hierarchy,0);

        for(int current_level=1;current_level<levels;++current_level)
            Copy_Unsubdivided_Cells<Struct_type,T,d>(graded_hierarchy->Allocator(current_level),graded_hierarchy->Blocks(current_level),
                                                     *minimal_hierarchy,nodes_of_cell_offsets,current_level);
        // flag ghost cells and update block offsets
        Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Cells(*minimal_hierarchy);
        minimal_hierarchy->Update_Block_Offsets();

        delete graded_hierarchy;
        current_hierarchy=minimal_hierarchy;}
}

template<class Struct_type,class Multigrid_struct_type,class T,int d>
void Check_Symmetry(Poisson_CG_System<Struct_type,Multigrid_struct_type,T,d>& cg_system)
{
    using Flags_type                            = typename Struct_type::Flags_type;
    using Allocator_type                        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;

    T Struct_type::* r1_channel                 = &Struct_type::ch6;
    T Struct_type::* r2_channel                 = &Struct_type::ch7;
    T Struct_type::* q1_channel                 = &Struct_type::ch8;
    T Struct_type::* q2_channel                 = &Struct_type::ch9;
    Grid_Hierarchy<Struct_type,T,d>& hierarchy  = cg_system.hierarchy;

    CG_Vector<Struct_type,T,d> r1(hierarchy,r1_channel);
    CG_Vector<Struct_type,T,d> r2(hierarchy,r2_channel);
    CG_Vector<Struct_type,T,d> q1(hierarchy,q1_channel);
    CG_Vector<Struct_type,T,d> q2(hierarchy,q2_channel);

    Random_Numbers<T> random;
    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto data1=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(r1_channel);
        auto data2=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(r2_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){
                    data1(offset)=random.Get_Uniform_Number((T)0.,(T)1.);
                    data2(offset)=random.Get_Uniform_Number((T)0.,(T)1.);}}}

    cg_system.Multiply(r1,q1);
    cg_system.Multiply(r2,q2);

    Log::cout<<"r1^T*A*r2: "<<cg_system.Inner_Product(r1,q2)<<std::endl;
    Log::cout<<"r2^T*A*r1: "<<cg_system.Inner_Product(r2,q1)<<std::endl;
}

template<class Base_struct_type,class Multigrid_struct_type,class T,int d>
void Check_Multigrid_Symmetry(Poisson_CG_System<Base_struct_type,Multigrid_struct_type,T,d>& cg_system,const int mg_levels)
{
    using Base_hierarchy                        = Grid_Hierarchy<Base_struct_type,T,d>;
    using Base_flags_type                       = typename Base_struct_type::Flags_type;
    using Base_allocator_type                   = SPGrid::SPGrid_Allocator<Base_struct_type,d>;
    using Base_flag_array_mask                  = typename Base_allocator_type::template Array_mask<unsigned>;

    T Base_struct_type::* r1_channel            = &Base_struct_type::ch6;
    T Base_struct_type::* r2_channel            = &Base_struct_type::ch7;
    T Base_struct_type::* q1_channel            = &Base_struct_type::ch8;
    T Base_struct_type::* q2_channel            = &Base_struct_type::ch9;
    Base_hierarchy& hierarchy                   = cg_system.hierarchy;

    CG_Vector<Base_struct_type,T,d> r1(hierarchy,r1_channel);
    CG_Vector<Base_struct_type,T,d> r2(hierarchy,r2_channel);
    CG_Vector<Base_struct_type,T,d> q1(hierarchy,q1_channel);
    CG_Vector<Base_struct_type,T,d> q2(hierarchy,q2_channel);

    Random_Numbers<T> random;
    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto data1=hierarchy.Allocator(level).template Get_Array<Base_struct_type,T>(r1_channel);
        auto data2=hierarchy.Allocator(level).template Get_Array<Base_struct_type,T>(r2_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Base_struct_type,unsigned>(&Base_struct_type::flags);

        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(unsigned e=0;e<Base_flag_array_mask::elements_per_block;++e,offset+=sizeof(Base_flags_type))
                if(flags(offset)&Cell_Type_Interior){
                    data1(offset)=random.Get_Uniform_Number((T)0.,(T)1.);
                    data2(offset)=random.Get_Uniform_Number((T)0.,(T)1.);}}}

    {
        Poisson_Multigrid_Solver<Base_struct_type,Multigrid_struct_type,T,d> multigrid_solver(hierarchy,mg_levels);
        multigrid_solver.Initialize();
        multigrid_solver.Initialize_Right_Hand_Side(r1_channel);
        multigrid_solver.Initialize_Guess();

        multigrid_solver.V_Cycle(5,1,200);
        multigrid_solver.Copy_Channel_Values(q1_channel,multigrid_solver.u_channel,false);
    }

    {
        Poisson_Multigrid_Solver<Base_struct_type,Multigrid_struct_type,T,d> multigrid_solver(hierarchy,mg_levels);
        multigrid_solver.Initialize();
        multigrid_solver.Initialize_Right_Hand_Side(r2_channel);
        multigrid_solver.Initialize_Guess();

        multigrid_solver.V_Cycle(5,1,200);
        multigrid_solver.Copy_Channel_Values(q2_channel,multigrid_solver.u_channel,false);
    }

    Log::cout<<"r1^T*A*r2: "<<cg_system.Inner_Product(r1,q2)<<std::endl;
    Log::cout<<"r2^T*A*r1: "<<cg_system.Inner_Product(r2,q1)<<std::endl;
}

int main(int argc,char** argv)
{
    enum {d=2};
    typedef float T;typedef Vector<T,d> TV;
    typedef Vector<int,d> T_INDEX;

    using Struct_type                           = Poisson_Data<T>;
    using Multigrid_struct_type                 = Multigrid_Data<T>;
    using Allocator_type                        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper                       = Grid_Topology_Helper<Flag_array_mask>;
    using Channel_Vector                        = Vector<T Struct_type::*,d>;
    using Hierarchy                             = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Lookup                      = Grid_Hierarchy_Lookup<Struct_type,T,d>;
    using Hierarchy_Rasterizer                  = Hierarchical_Rasterizer<Struct_type,T,d>;
    using Hierarchy_Visualization               = Grid_Hierarchy_Visualization<Struct_type,T>;
    enum {number_of_faces_per_cell              = Topology_Helper::number_of_faces_per_cell};

    Log::Initialize_Logging();

    Parse_Args parse_args;
    parse_args.Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args.Add_Integer_Argument("-test_number",1,"Test number.");
    parse_args.Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
    parse_args.Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args.Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args.Add_Option_Argument("-random_guess","Use random initial guess.");
    parse_args.Add_String_Argument("-solver","cg","Choice of solver.");
    parse_args.Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    if(d==2) parse_args.Add_Vector_2D_Argument("-size",Vector<double,2>(64.),"n","Grid resolution");
    else if(d==3) parse_args.Add_Vector_3D_Argument("-size",Vector<double,3>(64.),"n","Grid resolution");
    parse_args.Parse(argc,argv);

    int levels=parse_args.Get_Integer_Value("-levels");
    int test_number=parse_args.Get_Integer_Value("-test_number"),frame=0;
    number_of_threads=parse_args.Get_Integer_Value("-threads");
    if(number_of_threads) pthread_queue=new Pthread_Queue(number_of_threads);

    int mg_levels=parse_args.Get_Integer_Value("-mg_levels");
    int cg_iterations=parse_args.Get_Integer_Value("-cg_iterations");
    int cg_restart_iterations=parse_args.Get_Integer_Value("-cg_restart_iterations");
    bool random_guess=parse_args.Is_Value_Set("-random_guess");
    std::string solver=parse_args.Get_String_Value("-solver");
    omp_set_num_threads(number_of_threads);

    T_INDEX cell_counts;
    if(d==2){auto cell_counts_2d=parse_args.Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) cell_counts(v)=cell_counts_2d(v);}
    else{auto cell_counts_3d=parse_args.Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) cell_counts(v)=cell_counts_3d(v);}

    std::string output_directory="Test_"+std::to_string(test_number)+"_Resolution_"+std::to_string(cell_counts(0));
    File_Utilities::Create_Directory(output_directory);
    File_Utilities::Create_Directory(output_directory+"/common");
    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/info.nova-animation",std::to_string(frame));
    Log::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);

    std::string surface_directory="Surface_"+std::to_string(test_number)+"_Resolution_"+std::to_string(cell_counts(0));
    File_Utilities::Create_Directory(surface_directory);
    File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));

    Hierarchy *hierarchy = new Hierarchy(cell_counts,Range<T,d>(TV(-.5),TV(.5)),levels);
    Adaptive_Sphere_Rasterizer<Struct_type,T,d> rasterizer(*hierarchy,TV(),(T).1);
    //Central_Rasterizer<Struct_type,T,d> rasterizer(*hierarchy);
    //Randomized_Rasterizer<Struct_type,T,d> rasterizer(*hierarchy);
    for(Grid_Hierarchy_Iterator<d,Hierarchy_Rasterizer> iterator(hierarchy->Lattice(levels-1).Cell_Indices(),levels-1,rasterizer);iterator.Valid();iterator.Next());
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Cells(*hierarchy);
    hierarchy->Update_Block_Offsets();

    // grade the octree
    //Compute_Graded_Octree(hierarchy);

    Vector<Vector<bool,2>,d> domain_walls;
    if(test_number==1) for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=false;
    else if(test_number==2){for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=true;
        domain_walls(1)(1)=false;}
    Initialize_Dirichlet_Cells<Struct_type,T,d>(*hierarchy,domain_walls);

    T Struct_type::* x_channel          = &Struct_type::ch0;

    // compute volume-weighted divergence
    T Struct_type::* b_channel          = &Struct_type::ch1;

#if 1
    TV X=hierarchy->Lattice(0).domain.Center()*(T).25;
    X(1)=(T)-.4;
    const T value=(T)hierarchy->Lattice(0).one_over_dX.Product();
    uint64_t offset;int level;
    Hierarchy_Lookup::Cell_Lookup(*hierarchy,X,offset,level);
    hierarchy->Allocator(level).template Get_Array<Struct_type,T>(b_channel)(offset)=value;
#else
    Compute_Right_Hand_Side<Struct_type,T,d>(*hierarchy,b_channel,(T)-2*d);

    // initialize levelset channel
    T Struct_type::* levelset_channel   = &Struct_type::ch9;
    Analytic_Levelset<T,d> *levelset = new Sphere_Levelset<T,d>(TV(),(T).25);
    for(unsigned level=0;level<levels;++level)
        Initialize_Levelset<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),levelset_channel,levelset,level);
    delete levelset;

    Initialize_Dirichlet_Cells(*hierarchy,levelset_channel);
#endif

    Channel_Vector gradient_channels;
    gradient_channels(0)                = &Struct_type::ch2;
    gradient_channels(1)                = &Struct_type::ch3;
    if(d==3) gradient_channels(2)       = &Struct_type::ch4;

    Poisson_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,3,1,200);

    Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);

#if 0
    Check_Multigrid_Symmetry<Struct_type,Multigrid_struct_type,T,d>(cg_system,mg_levels);
#endif

    if(solver=="mgpcg"){
        T Struct_type::* q_channel      = &Struct_type::ch5;
        T Struct_type::* r_channel      = &Struct_type::ch6;
        T Struct_type::* s_channel      = &Struct_type::ch7;
        T Struct_type::* k_channel      = &Struct_type::ch7;
        T Struct_type::* z_channel      = &Struct_type::ch8;

#if 0
        Check_Symmetry(cg_system);
#endif
        // clear all channels
        for(int level=0;level<levels;++level){
            Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
            Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
            Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
            Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);}

        CG_Vector<Struct_type,T,d> x_V(*hierarchy,x_channel);
        CG_Vector<Struct_type,T,d> b_V(*hierarchy,b_channel);
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
        cg.print_residuals=true;
        cg.print_diagnostics=true;
        cg.restart_iterations=cg_restart_iterations;
        const T tolerance=std::max((T)1e-6*b_norm,(T)1e-6);
        cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);

        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(++frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);}
    else if(solver=="smoother") for(int i=1;i<=cg_iterations;++i){++frame;
        T Struct_type::* temp_channel   = &Struct_type::ch5;
        Poisson_Multigrid_Smoother<Struct_type,T,d>::Exact_Solve(*hierarchy,gradient_channels,x_channel,b_channel,
                                                         temp_channel,1,(unsigned)Cell_Type_Interior);
        Poisson_Multigrid_Smoother<Struct_type,T,d>::Compute_Residual(*hierarchy,gradient_channels,x_channel,b_channel,
                                                              temp_channel,(unsigned)Cell_Type_Interior);
        CG_Vector<Struct_type,T,d> r_V(*hierarchy,temp_channel);
        Log::cout<<cg_system.Convergence_Norm(r_V)<<std::endl;
        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);}
    else if(solver=="mg"){Initialize_Guess(*hierarchy,x_channel,random_guess);
        Poisson_Multigrid_Solver<Struct_type,Multigrid_struct_type,T,d> multigrid_solver(*hierarchy,mg_levels);
        multigrid_solver.Initialize();
        multigrid_solver.Initialize_Right_Hand_Side(b_channel);
        multigrid_solver.Initialize_Guess();
        if(random_guess) multigrid_solver.Copy_Channel_Values(x_channel,multigrid_solver.u_channel);

        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);
        multigrid_solver.Compute_Residual(0);
        Log::cout<<multigrid_solver.Convergence_Norm(0,multigrid_solver.temp_channel)<<std::endl;

        for(int i=1;i<=cg_iterations;++i){multigrid_solver.V_Cycle(5,1,200);
            multigrid_solver.Copy_Channel_Values(x_channel,multigrid_solver.u_channel,false);
            File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(++frame));
            File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
            Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);

            multigrid_solver.Compute_Residual(0);
            Log::cout<<multigrid_solver.Convergence_Norm(0,multigrid_solver.temp_channel)<<std::endl;}}

    //Compute_Error<Struct_type,T,d>(*hierarchy,x_channel);

    delete hierarchy;

    Log::Finish_Logging();
    return 0;
}

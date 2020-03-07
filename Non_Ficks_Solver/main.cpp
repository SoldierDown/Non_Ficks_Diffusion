//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Hierarchy/Visualization/Grid_Hierarchy_Visualization.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Log/Log.h>
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/Tools/Utilities/Constants.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Pthread_Queue.h>
#include "Non_Ficks_CG_System.h"
#include "../Initialize_Dirichlet_Cells.h"
#include "../MPM_Data.h"
#include "../Multigrid_Solver/Multigrid_Data.h"
#include "../Rasterizers/Adaptive_Sphere_Rasterizer.h"
#include <omp.h>

using namespace Nova;
using namespace SPGrid;

extern Pthread_Queue* pthread_queue;

namespace Nova{
int number_of_threads=0;
}

int main(int argc,char** argv)
{
    enum {d=2};
    typedef float T;typedef Vector<T,d> TV;
    typedef Vector<int,d> T_INDEX;

    using Struct_type                           = MPM_Data<T>;
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

    std::string surface_directory="Surface_"+std::to_string(test_number)+"_Resolution_"+std::to_string(cell_counts(0));
    File_Utilities::Create_Directory(surface_directory);
    File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));

    Hierarchy *hierarchy = new Hierarchy(cell_counts,Range<T,d>(TV(-.5),TV(.5)),levels);
    Adaptive_Sphere_Rasterizer<Struct_type,T,d> rasterizer(*hierarchy,TV(),(T).1);
    for(Grid_Hierarchy_Iterator<d,Hierarchy_Rasterizer> iterator(hierarchy->Lattice(levels-1).Cell_Indices(),levels-1,rasterizer);iterator.Valid();iterator.Next());
    hierarchy->Update_Block_Offsets();

    Vector<Vector<bool,2>,d> domain_walls;
    if(test_number==1) for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=false;
    else if(test_number==2){for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=true;
        domain_walls(1)(1)=false;}
    Initialize_Dirichlet_Cells<Struct_type,T,d>(*hierarchy,domain_walls);

    T Struct_type::* x_channel          = &Struct_type::ch0;
    T Struct_type::* b_channel          = &Struct_type::ch1;
    
    const T diff_coeff=(T)1.; const T dt=(T)1; const T Fc=(T)0.; const T tau=(T)1.;
    const T coeff1=dt*diff_coeff*(Fc*tau+dt)/(dt+tau);
    // initialize right hand side
    TV X=hierarchy->Lattice(0).domain.Center()*(T).25;
    X(1)=(T)-.4;
    const T value=hierarchy->Lattice(0).one_over_dX.Product();
    uint64_t offset;int level;
    Hierarchy_Lookup::Cell_Lookup(*hierarchy,X,offset,level);
    hierarchy->Allocator(level).template Get_Array<Struct_type,T>(b_channel)(offset)=value;

    Non_Ficks_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,coeff1,3,1,200);

    Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);

    if(solver=="mgpcg"){
        T Struct_type::* q_channel      = &Struct_type::ch5;
        T Struct_type::* r_channel      = &Struct_type::ch6;
        T Struct_type::* s_channel      = &Struct_type::ch7;
        T Struct_type::* k_channel      = &Struct_type::ch7;
        T Struct_type::* z_channel      = &Struct_type::ch8;

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
        Non_Ficks_Smoother<Struct_type,T,d>::Exact_Solve(*hierarchy,x_channel,b_channel,temp_channel,coeff1,1,(unsigned)Cell_Type_Interior);
        Non_Ficks_Smoother<Struct_type,T,d>::Compute_Residual(*hierarchy,x_channel,b_channel,temp_channel,coeff1,(unsigned)Cell_Type_Interior);
        CG_Vector<Struct_type,T,d> r_V(*hierarchy,temp_channel);
        Log::cout<<cg_system.Convergence_Norm(r_V)<<std::endl;
        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);}

    delete hierarchy;

    Log::Finish_Logging();
    return 0;
}

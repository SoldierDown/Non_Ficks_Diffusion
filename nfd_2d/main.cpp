//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <chrono>
#include "../MPM_Driver.h"
#include "Standard_Tests/Standard_Tests.h"

// for test
#include <nova/Tools/Utilities/Pthread_Queue.h>
#include <nova/Dynamics/Hierarchy/Visualization/Grid_Hierarchy_Visualization.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Tools/Grids/Grid_Iterator_Cell.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <omp.h>
#include "Helper.h"

extern Pthread_Queue* pthread_queue;
using namespace Nova;
using namespace std::chrono;

namespace Nova{
    int number_of_threads=0;
}


int main(int argc,char** argv)
{
    Log::cout.precision(13);
    bool run_test=true;
    if(run_test){
        typedef float T;
        enum {d=2};
        typedef Vector<T,d> TV;
        typedef Vector<int,d> T_INDEX;

        using Cell_Iterator                         = Grid_Iterator_Cell<T,d>;
        using Struct_type                           = MPM_Data<T>;
        using Multigrid_struct_type                 = Multigrid_Data<T>;
        using Allocator_type                        = SPGrid::SPGrid_Allocator<Struct_type,d>;
        using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
        using Topology_Helper                       = Grid_Topology_Helper<Flag_array_mask>;
        using Channel_Vector                        = Vector<T Struct_type::*,d>;
        using Hierarchy                             = Grid_Hierarchy<Struct_type,T,d>;
        using Hierarchy_Visualization               = Grid_Hierarchy_Visualization<Struct_type,T>;
        
        Log::Initialize_Logging();

        Parse_Args parse_args;
        parse_args.Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
        parse_args.Add_Integer_Argument("-sm_iterations",1,"Number of smoother iterations.");
        parse_args.Add_Integer_Argument("-boundary_iterations",5,"Number of boundary iterations.");
        parse_args.Add_Integer_Argument("-interior_iterations",1,"Number of interior iterations.");
        parse_args.Add_Integer_Argument("-bottom_iterations",200,"Number of bottom iterations.");
        parse_args.Add_Integer_Argument("-test_number",1,"Test number.");
        parse_args.Add_Integer_Argument("-mg_levels",2,"Number of levels in the Multigrid hierarchy.");
        parse_args.Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
        parse_args.Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
        parse_args.Add_Option_Argument("-random_guess","Use random initial guess.");
        parse_args.Add_Option_Argument("-bs","Run for boundary only");
        parse_args.Add_Option_Argument("-is","Run for interior only");
        parse_args.Add_Option_Argument("-simple_case","Ran simple case.");
        parse_args.Add_Option_Argument("-ficks","Fick's diffusion.");
        parse_args.Add_String_Argument("-solver","cg","Choice of solver.");
        parse_args.Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
        if(d==2) parse_args.Add_Vector_2D_Argument("-size",Vector<double,2>(64.),"n","Grid resolution");
        else parse_args.Add_Vector_3D_Argument("-size",Vector<double,3>(64.),"n","Grid resolution");
        parse_args.Parse(argc,argv);


        int sm_iterations=parse_args.Get_Integer_Value("-sm_iterations");
        int boundary_iterations=parse_args.Get_Integer_Value("-boundary_iterations");
        int interior_iterations=parse_args.Get_Integer_Value("-interior_iterations");
        int bottom_iterations=parse_args.Get_Integer_Value("-bottom_iterations");
        bool simple_case=parse_args.Get_Option_Value("-simple_case");
        bool random_guess=parse_args.Get_Option_Value("-random_guess");
        bool bs=parse_args.Get_Option_Value("-bs");
        bool is=parse_args.Get_Option_Value("-is");
        bool FICKS=parse_args.Get_Option_Value("-ficks");
        int levels=parse_args.Get_Integer_Value("-levels");
        int test_number=parse_args.Get_Integer_Value("-test_number"),frame=0;
        number_of_threads=parse_args.Get_Integer_Value("-threads");
        if(number_of_threads) pthread_queue=new Pthread_Queue(number_of_threads);

        int mg_levels=parse_args.Get_Integer_Value("-mg_levels");
        int cg_iterations=parse_args.Get_Integer_Value("-cg_iterations");
        int cg_restart_iterations=parse_args.Get_Integer_Value("-cg_restart_iterations");
        omp_set_num_threads(number_of_threads);

        T_INDEX cell_counts;
        if(d==2) {auto cell_counts_2d=parse_args.Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) cell_counts(v)=cell_counts_2d(v);}
        else if(d==3) {auto cell_counts_3d=parse_args.Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) cell_counts(v)=cell_counts_3d(v);}


        std::string output_directory="Test_"+std::to_string(test_number)+"_Resolution_"+std::to_string(cell_counts(0));
        File_Utilities::Create_Directory(output_directory);
        File_Utilities::Create_Directory(output_directory+"/common");
        File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(output_directory+"/info.nova-animation",std::to_string(frame));
        Log::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);

        std::string surface_directory=std::to_string(d)+"d_"+std::to_string(mg_levels)+(FICKS?"levels_F_":"levels_NF_")+(simple_case?"simple_case_":"complex_case_")+(random_guess?"random_init_":"0_init_")+"Resolution_"+std::to_string(cell_counts(0));
        surface_directory="V_test";
        File_Utilities::Create_Directory(surface_directory);
        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));

        T Struct_type::* x_channel                = &Struct_type::ch0;
        T Struct_type::* b_channel                = &Struct_type::ch1;
        T Struct_type::* r_channel                = &Struct_type::ch2;

        const T diff_coeff=(T)1e-3; 
        const T dt=(T)1e-3; const T Fc=(T)0.; const T tau=(T)1.; 
        
        Hierarchy *hierarchy=new Hierarchy(cell_counts,Range<T,d>(TV(-1),TV(1)),levels);            

        
        const Grid<T,d>& grid=hierarchy->Lattice(0);
        Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(TV(-1)),grid.Clamp_To_Cell(TV(1)));
        for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next()){
            T_INDEX cell_index=iterator.Cell_Index();
            if(simple_case){
                if(cell_index(0)==1||cell_index(1)==1||cell_index(0)==cell_counts(0)||cell_index(1)==cell_counts(1)) hierarchy->Activate_Cell(0,cell_index,Cell_Type_Dirichlet);
                else hierarchy->Activate_Cell(0,cell_index,Cell_Type_Interior);}
            else{
                if(cell_index(0)!=1&&cell_index(0)!=cell_counts(0)&&cell_index(1)!=cell_counts(1)) hierarchy->Activate_Cell(0,cell_index,Cell_Type_Interior);// set as exteriror do nothing
                else if(cell_index(1)==1) hierarchy->Activate_Cell(0,cell_index,Cell_Type_Dirichlet);}}
        

        hierarchy->Update_Block_Offsets();
        hierarchy->Initialize_Red_Black_Partition(2*number_of_threads);
        Sphere_Levelset<T,d> *levelset=new Sphere_Levelset<T,d>(TV(),(T).25);
        Levelset_Initializer<Struct_type,T,d>(*hierarchy,hierarchy->Allocator(0),hierarchy->Blocks(0),x_channel,*levelset);
        delete levelset;
        if(!simple_case)
            for(int level=0;level<levels;++level) Neumann_BC_Initializer<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel);
        for(int level=0;level<levels;++level){
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel);
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),b_channel);
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);}

        Initialize_Guess(*hierarchy,x_channel,false);
        const T_INDEX pin_cell=T_INDEX(int(cell_counts(0)/6));
        for(int level=0;level<levels;++level) hierarchy->Channel(level,x_channel)(pin_cell._data)=(T)10.;
        Compute_Right_Hand_Side(*hierarchy,x_channel,b_channel,FICKS,dt,diff_coeff,Fc,tau);

        Initialize_Guess(*hierarchy,x_channel,random_guess);
        Multigrid_Solver<Struct_type,Multigrid_struct_type,T,d> multigrid_solver(*hierarchy,mg_levels,FICKS,dt,diff_coeff,Fc,tau);
        multigrid_solver.Initialize();
        multigrid_solver.Initialize_Right_Hand_Side(b_channel);
        multigrid_solver.Initialize_Guess();
        multigrid_solver.Copy_Channel_Values(x_channel,multigrid_solver.x_channel);
        
        Log::cout<<"rhs norm: "<<multigrid_solver.Convergence_Norm(0,multigrid_solver.b_channel)<<std::endl;
        Log::cout<<"init guess norm: "<<multigrid_solver.Convergence_Norm(0,multigrid_solver.x_channel)<<std::endl;

        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);
        multigrid_solver.Compute_Residual(0);
        Log::cout<<multigrid_solver.Convergence_Norm(0,multigrid_solver.temp_channel)<<std::endl;

        frame++;
  
        // for(int level=0;level<levels;++level) Multigrid_Smoother<Multigrid_struct_type,T,d>::Jacobi_Iteration(*(multigrid_solver.multigrid_hierarchy(level)),multigrid_solver.multigrid_hierarchy(level)->Blocks(level),level,multigrid_solver.x_channel,multigrid_solver.b_channel,multigrid_solver.temp_channel,bottom_iterations,Cell_Type_Interior,FICKS,dt,diff_coeff,Fc,tau);  
        multigrid_solver.V_Cycle(boundary_iterations,interior_iterations,bottom_iterations);
        multigrid_solver.Copy_Channel_Values(x_channel,multigrid_solver.x_channel,false);
        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);
        multigrid_solver.Compute_Residual(0);
        Log::cout<<multigrid_solver.Convergence_Norm(0,multigrid_solver.temp_channel)<<std::endl;

        // Log::cout<<"sm iterations: "<<sm_iterations<<std::endl;
        // for(int i=1;i<=cg_iterations;++i){ ++frame;
            // if(bs) for(int level=0;level<levels;++level) Multigrid_Smoother<Multigrid_Struct_Type,T,d>::Jacobi_Iteration(*hierarchy,hierarchy->Boundary_Blocks(level),level,x_channel,b_channel,r_channel,sm_iterations,(unsigned)MG_Boundary,FICKS,a,twod_a_plus_one,coeff1);
            // if(is) for(int level=0;level<levels;++level) Multigrid_Smoother<Multigrid_Struct_Type,T,d>::Jacobi_Iteration(*hierarchy,hierarchy->Blocks(level),level,x_channel,b_channel,r_channel,sm_iterations,(unsigned)Cell_Type_Interior,FICKS,a,twod_a_plus_one,coeff1);
            // Multigrid_Smoother<Multigrid_Struct_Type,T,d>::Compute_Residual(*hierarchy,x_channel,b_channel,r_channel,(unsigned)Cell_Type_Interior,FICKS,a,twod_a_plus_one,coeff1);
            // for(int level=0;level<levels;++level) Saturation_Clamp_Heler<Multigrid_Struct_Type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel);     
            
            // Diffusion_CG_Vector<Multigrid_Struct_Type,T,d> r_V(*hierarchy,r_channel);
            // Diffusion_CG_Vector<Multigrid_Struct_Type,T,d> x_V(*hierarchy,x_channel);
            // Log::cout<<"residual norm: "<<cg_system.Convergence_Norm(r_V)<<std::endl;
            // Log::cout<<"x norm: "<<cg_system.Convergence_Norm(x_V)<<std::endl;
            // File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
            // File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
            // Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);
        // }
        delete hierarchy;
        Log::Finish_Logging();
        return 0;
    }
    high_resolution_clock::time_point tb = high_resolution_clock::now();
    enum {d=2};
    typedef float T;
    MPM_Example<T,d> *example=new Standard_Tests<T,d>();
    example->Parse(argc,argv);
    example->bbox=Range<T,d>(example->domain.max_corner,example->domain.min_corner);
    File_Utilities::Create_Directory(example->output_directory);
    File_Utilities::Create_Directory(example->output_directory+"/common");
    Log::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
        
    MPM_Driver<T,d> driver(*example);
    driver.Execute_Main_Program();
    high_resolution_clock::time_point te=high_resolution_clock::now();
	duration<double> dur=duration_cast<duration<double>>(te-tb);
	Log::cout<<"Total duration: "<<dur.count()<<std::endl;        
    delete example;

    return 0;
}

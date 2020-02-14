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
#include "../Smoother/Multigrid_Data.h"
#include "../Smoother/Multigrid_Smoother.h"
#include "../Diffusion_Helper/Diffusion_CG_System.h"
#include "../Diffusion_Helper/Diffusion_CG_Vector.h"
#include "../Ficks_RHS_Helper.h"
#include "../Non_Ficks_RHS_Helper.h"
#include "../Saturation_Clamp_Helper.h"
#include "../Smoother/Div_Qc_Updater.h"
#include "../Smoother/Initial_Guess_Helper.h"
#include "../Smoother/Sphere_Levelset.h"
#include "../Smoother/Levelset_Initializer.h"
#include "../Smoother/Neumann_BC_Initializer.h"
#include "../Smoother/Mark_Boundary.h"
#include <omp.h>

extern Pthread_Queue* pthread_queue;
using namespace Nova;
using namespace std::chrono;

int main(int argc,char** argv)
{
    bool run_test=true;
    if(run_test){
        typedef float T;
        enum {d=3};
        typedef Vector<T,d> TV;
        typedef Vector<int,d> T_INDEX;

        using Cell_Iterator                         = Grid_Iterator_Cell<T,d>;
        using Multigrid_struct_type                 = Multigrid_Data<T>;
        using Multigrid_allocator_type              = SPGrid::SPGrid_Allocator<Multigrid_struct_type,d>;
        using Flag_array_mask                       = typename Multigrid_allocator_type::template Array_mask<unsigned>;
        using Topology_Helper                       = Grid_Topology_Helper<Flag_array_mask>;
        using Channel_Vector                        = Vector<T Multigrid_struct_type::*,d>;
        using Hierarchy                             = Grid_Hierarchy<Multigrid_struct_type,T,d>;
        using Hierarchy_Initializer                 = Grid_Hierarchy_Initializer<Multigrid_struct_type,T,d>;
        using Hierarchy_Visualization               = Grid_Hierarchy_Visualization<Multigrid_struct_type,T>;
        using Multigrid_flag_array_mask             = typename Multigrid_allocator_type::template Array_mask<unsigned>;
        
        enum {number_of_faces_per_cell              = Topology_Helper::number_of_faces_per_cell};

        Log::Initialize_Logging();

        Parse_Args parse_args;
        parse_args.Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
        parse_args.Add_Integer_Argument("-sm_iterations",1,"Number of smoother iterations.");
        parse_args.Add_Integer_Argument("-test_number",1,"Test number.");
        parse_args.Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
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
        bool simple_case=parse_args.Get_Option_Value("-simple_case");
        bool random_guess=parse_args.Get_Option_Value("-random_guess");
        bool bs=parse_args.Get_Option_Value("-bs");
        bool is=parse_args.Get_Option_Value("-is");
        bool FICKS=parse_args.Get_Option_Value("-ficks");
        int levels=parse_args.Get_Integer_Value("-levels");
        int test_number=parse_args.Get_Integer_Value("-test_number"),frame=0;
        int number_of_threads=parse_args.Get_Integer_Value("-threads");
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

        std::string surface_directory=std::to_string(d)+(FICKS?"d_F_":"d_NF_")+(simple_case?"simple_case_":"complex_case_")+(random_guess?"random_init_":"0_init_")+"Resolution_"+std::to_string(cell_counts(0));
        surface_directory="Slice";
        File_Utilities::Create_Directory(surface_directory);
        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));

        T Multigrid_struct_type::* x_channel                = &Multigrid_struct_type::ch0;
        T Multigrid_struct_type::* b_channel                = &Multigrid_struct_type::ch1;
        T Multigrid_struct_type::* r_channel                = &Multigrid_struct_type::ch2;

        Hierarchy *hierarchy=new Hierarchy(cell_counts,Range<T,d>(TV(-1),TV(1)),levels);
        Sphere_Levelset<T,d>* levelset=new Sphere_Levelset<T,d>(TV(),(T).25);        
        // reuse x_channel to set up level set
        for(int level=0;level<levels;++level) Levelset_Initializer<Multigrid_struct_type,T,d>(hierarchy->Lattice(level),hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel,levelset);
        delete levelset;
        const Grid<T,d>& grid=hierarchy->Lattice(0);
        Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(TV(-1)),grid.Clamp_To_Cell(TV(1)));
        for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next()){
            T_INDEX cell_index=iterator.Cell_Index();
            if(cell_index(0)==1||cell_index(1)==1||cell_index(0)==cell_counts(0)||cell_index(1)==cell_counts(1)) hierarchy->Activate_Cell(0,cell_index,Cell_Type_Dirichlet);
            else hierarchy->Activate_Cell(0,cell_index,Cell_Type_Interior);}
        if(!simple_case){for(int level=0;level<levels;++level) Neumann_BC_Initializer<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel);}
        hierarchy->Update_Block_Offsets();
        hierarchy->Initialize_Red_Black_Partition(2*number_of_threads);

        for(int level=0;level<levels;++level){
            SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel);
            SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),b_channel);
            SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
        }

        int boundary_radius=3;
        Array<uint64_t> neighbor_offsets;
        const T_INDEX boundary_radius_vector(boundary_radius),zero_vector=T_INDEX();
        for(Range_Iterator<d> iterator(-boundary_radius_vector,boundary_radius_vector);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Index();
            if(index!=zero_vector) neighbor_offsets.Append(Multigrid_flag_array_mask::Linear_Offset(index._data));}
        for(int level=0;level<levels;++level)
                Mark_Boundary<Multigrid_struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),neighbor_offsets,level,(unsigned)MG_Boundary);
        hierarchy->Initialize_Boundary_Blocks((unsigned)MG_Boundary);

        // saturation_channel
        // velocity_channels as gradient_channels
        // mass_channel as result_channel
        const T diff_coeff=(T)1e-3; const T one_over_dx2=grid.one_over_dX(0)*grid.one_over_dX(1); 
        const T dt=(T)1e-3; const T Fc=(T)0.; const T tau=(T)0.; const T a=diff_coeff*dt*one_over_dx2; 
        const T coeff1=dt*diff_coeff*(Fc*tau+dt)*one_over_dx2/(dt+tau);
        const T coeff2=dt*tau/(dt+tau);
        const T coeff3=dt*diff_coeff*(1-Fc)/(dt+tau);
        const T coeff4=tau/(dt+tau);
        const T twod_a_plus_one=(T)2.*d*a+(T)1.;
        const T_INDEX pin_cell=T_INDEX(10);
        Log::cout<<"dt: "<<dt<<", diff_coeff: "<<diff_coeff<<", one_over_dx2: "<<one_over_dx2<<", tau: "<<tau<<", a: "<<a<<", coeff1: "<<coeff1<<", twod_a_plus_one: "<<twod_a_plus_one<<std::endl;
        for(int level=0;level<levels;++level) hierarchy->Channel(level,b_channel)(pin_cell._data)=(T)1.;
        for(int level=0;level<levels;++level) Initial_Guess_Helper<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel,random_guess);     
        
        Diffusion_CG_System<Multigrid_struct_type,T,d> cg_system(*hierarchy,FICKS);
        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);

        // write hierarchy
        File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
        hierarchy->Write_Hierarchy(output_directory,frame);

        Diffusion_CG_Vector<Multigrid_struct_type,T,d> r_V_before(*hierarchy,b_channel);
        Log::cout<<"rhs norm: "<<cg_system.Convergence_Norm(r_V_before)<<std::endl;
        
        Log::cout<<"sm iterations: "<<sm_iterations<<std::endl;
        for(int i=1;i<=cg_iterations;++i){ ++frame;
            if(bs) for(int level=0;level<levels;++level) Multigrid_Smoother<Multigrid_struct_type,T,d>::Jacobi_Iteration(*hierarchy,hierarchy->Boundary_Blocks(level),level,x_channel,b_channel,r_channel,sm_iterations,(unsigned)MG_Boundary,FICKS,a,twod_a_plus_one,coeff1);
            if(is) for(int level=0;level<levels;++level) Multigrid_Smoother<Multigrid_struct_type,T,d>::Jacobi_Iteration(*hierarchy,hierarchy->Blocks(level),level,x_channel,b_channel,r_channel,sm_iterations,(unsigned)Cell_Type_Interior,FICKS,a,twod_a_plus_one,coeff1);
            Multigrid_Smoother<Multigrid_struct_type,T,d>::Compute_Residual(*hierarchy,x_channel,b_channel,r_channel,(unsigned)Cell_Type_Interior,FICKS,a,twod_a_plus_one,coeff1);
            for(int level=0;level<levels;++level) Saturation_Clamp_Heler<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),x_channel);     
            
            Diffusion_CG_Vector<Multigrid_struct_type,T,d> r_V(*hierarchy,r_channel);
            Diffusion_CG_Vector<Multigrid_struct_type,T,d> x_V(*hierarchy,x_channel);
            Log::cout<<"residual norm: "<<cg_system.Convergence_Norm(r_V)<<std::endl;
            Log::cout<<"x norm: "<<cg_system.Convergence_Norm(x_V)<<std::endl;
            File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
            File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
            Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,x_channel,surface_directory,frame);
        }
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

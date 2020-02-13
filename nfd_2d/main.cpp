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
#include <omp.h>

extern Pthread_Queue* pthread_queue;
using namespace Nova;
using namespace std::chrono;

int main(int argc,char** argv)
{
    bool run_test=false;
    if(run_test){
        typedef float T;
        enum {d=3};
        typedef Vector<T,d> TV;
        typedef Vector<int,d> T_INDEX;

        using Cell_Iterator                         = Grid_Iterator_Cell<T,d>;
        using Multigrid_struct_type                 = Multigrid_Data<T>;
        using Allocator_type                        = SPGrid::SPGrid_Allocator<Multigrid_struct_type,d>;
        using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
        using Topology_Helper                       = Grid_Topology_Helper<Flag_array_mask>;
        using Channel_Vector                        = Vector<T Multigrid_struct_type::*,d>;
        using Hierarchy                             = Grid_Hierarchy<Multigrid_struct_type,T,d>;
        using Hierarchy_Initializer                 = Grid_Hierarchy_Initializer<Multigrid_struct_type,T,d>;
        using Hierarchy_Visualization               = Grid_Hierarchy_Visualization<Multigrid_struct_type,T>;
        enum {number_of_faces_per_cell              = Topology_Helper::number_of_faces_per_cell};

        Log::Initialize_Logging();

        Parse_Args parse_args;
        parse_args.Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
        parse_args.Add_Integer_Argument("-test_number",1,"Test number.");
        parse_args.Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
        parse_args.Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
        parse_args.Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
        parse_args.Add_Option_Argument("-random_guess","Use random initial guess.");
        parse_args.Add_Option_Argument("-simple_case","Ran simple case.");
        parse_args.Add_Option_Argument("-ficks","Fick's diffusion.");
        parse_args.Add_String_Argument("-solver","cg","Choice of solver.");
        parse_args.Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
        if(d==2) parse_args.Add_Vector_2D_Argument("-size",Vector<double,2>(64.),"n","Grid resolution");
        else parse_args.Add_Vector_3D_Argument("-size",Vector<double,3>(64.),"n","Grid resolution");
        parse_args.Parse(argc,argv);

        bool simple_case=parse_args.Get_Option_Value("-simple_case");
        bool random_guess=parse_args.Get_Option_Value("-random_guess");
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
        File_Utilities::Create_Directory(surface_directory);
        File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
        File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));

        T Multigrid_struct_type::* saturation_channel       = &Multigrid_struct_type::ch0;
        T Multigrid_struct_type::* div_Qc_channel           = &Multigrid_struct_type::ch1;
        T Multigrid_struct_type::* rhs_channel              = &Multigrid_struct_type::ch2;
        T Multigrid_struct_type::* result_channel           = &Multigrid_struct_type::ch3;
        Channel_Vector gradient_channels; 
        gradient_channels(0)                                = &Multigrid_struct_type::ch4;
        gradient_channels(1)                                = &Multigrid_struct_type::ch5;
        if(d==3) gradient_channels(2)                       = &Multigrid_struct_type::ch6;

        Hierarchy *hierarchy=new Hierarchy(cell_counts,Range<T,d>(TV(-1),TV(1)),levels);
        Sphere_Levelset<T,d>* levelset=new Sphere_Levelset<T,d>(TV(),(T).25);        
        for(int level=0;level<levels;++level) Levelset_Initializer<Multigrid_struct_type,T,d>(hierarchy->Lattice(level),hierarchy->Allocator(level),hierarchy->Blocks(level),div_Qc_channel,levelset);
        delete levelset;
        const Grid<T,d>& grid=hierarchy->Lattice(0);
        Range<int,d> bounding_grid_cells(grid.Clamp_To_Cell(TV(-1)),grid.Clamp_To_Cell(TV(1)));
        for(Cell_Iterator iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next()){
            T_INDEX cell_index=iterator.Cell_Index();
            if(cell_index(0)==1||cell_index(1)==1||cell_index(0)==cell_counts(0)||cell_index(1)==cell_counts(1)) hierarchy->Activate_Cell(0,cell_index,Cell_Type_Dirichlet);
            else hierarchy->Activate_Cell(0,cell_index,Cell_Type_Interior);}
        if(!simple_case){for(int level=0;level<levels;++level) Neumann_BC_Initializer<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_Qc_channel);}
        hierarchy->Update_Block_Offsets();
        hierarchy->Initialize_Red_Black_Partition(2*number_of_threads);




        for(int level=0;level<levels;++level){
            SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel);
            SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),div_Qc_channel);
            SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),rhs_channel);
            SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),result_channel);
            for(int v=0;v<d;++v) SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),gradient_channels(v));
        }


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
        for(int level=0;level<levels;++level) {hierarchy->Channel(level,saturation_channel)(pin_cell._data)=(T)1.;hierarchy->Channel(level,rhs_channel)(pin_cell._data)=(T)1.;}
        if(FICKS) for(int level=0;level<levels;++level) Ficks_RHS_Helper<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel,rhs_channel,a);     
        else for(int level=0;level<levels;++level) Non_Ficks_RHS_Helper<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel,div_Qc_channel,rhs_channel,coeff1,coeff2);     
        


        Diffusion_CG_System<Multigrid_struct_type,T,d> cg_system(*hierarchy,FICKS);
        Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,saturation_channel,surface_directory,frame);

        // write hierarchy
        File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
        hierarchy->Write_Hierarchy(output_directory,frame);



        Diffusion_CG_Vector<Multigrid_struct_type,T,d> r_V_before(*hierarchy,rhs_channel);
        Log::cout<<"rhs norm: "<<cg_system.Convergence_Norm(r_V_before)<<std::endl;
        
        // initial guess
        for(int level=0;level<levels;++level) Initial_Guess_Helper<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel,random_guess);     
            

        for(int i=1;i<=cg_iterations;++i){ ++frame;
            for(int level=0;level<levels;++level){
                SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),result_channel);
                SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),gradient_channels(0));
                SPGrid::Clear<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),gradient_channels(1));   
            }

            Multigrid_Smoother<Multigrid_struct_type,T,d>::Exact_Solve(*hierarchy,gradient_channels,saturation_channel,rhs_channel,
                                                             result_channel,1,(unsigned)Cell_Type_Interior,FICKS,a,twod_a_plus_one,coeff1);
            Multigrid_Smoother<Multigrid_struct_type,T,d>::Compute_Residual(*hierarchy,gradient_channels,saturation_channel,rhs_channel,
                                                                  result_channel,(unsigned)Cell_Type_Interior,FICKS,a,twod_a_plus_one,coeff1);
            for(int level=0;level<levels;++level) Saturation_Clamp_Heler<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel);     
            // update Qc
            if(!FICKS) for(int level=0;level<levels;++level) Div_Qc_Updater<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel,div_Qc_channel,coeff3,coeff4,one_over_dx2);     
            
            Diffusion_CG_Vector<Multigrid_struct_type,T,d> r_V(*hierarchy,result_channel);
            Log::cout<<cg_system.Convergence_Norm(r_V)<<std::endl;
            File_Utilities::Create_Directory(surface_directory+"/"+std::to_string(frame));
            File_Utilities::Write_To_Text_File(surface_directory+"/info.nova-animation",std::to_string(frame));
            Hierarchy_Visualization::Visualize_Heightfield(*hierarchy,saturation_channel,surface_directory,frame);
            if(FICKS) for(int level=0;level<levels;++level) Ficks_RHS_Helper<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel,rhs_channel,a);     
            else for(int level=0;level<levels;++level) Non_Ficks_RHS_Helper<Multigrid_struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),saturation_channel,div_Qc_channel,rhs_channel,coeff1,coeff2);     

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

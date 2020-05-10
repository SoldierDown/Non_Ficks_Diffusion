//!#####################################################################
//! \file OpenVDB_Converter.h
//!#####################################################################
// Class OpenVDB_Converter
//######################################################################
#ifndef __OpenVDB_Converter__
#define __OpenVDB_Converter__
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Dynamics/Rigid_Bodies/Rigid_Body.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include <assert.h>
#include <omp.h>
#include <openvdb/openvdb.h>
#include "Viewer_Data.h"
#include "../Poisson_Solver/Grid_Hierarchy_Projection.h"
#include <iostream>
using namespace std;

namespace Nova{
template<class T>
inline T clamp(T a, T lower, T upper)
{   
    if(a<lower) return lower;
    else if(a>upper) return upper;
    else return a;
}
template<class T,int d>
class OpenVDB_Converter
{
    using TV                                = Vector<T,d>;
    using T_INDEX                           = Vector<int,d>;
    using Struct_type                       = Viewer_Data<T>;
    using Index_type                        = std::array<SPGrid::ucoord_t,d>;
    using Channel_Vector                    = Vector<T Struct_type::*,d>;
    using Hierarchy                         = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Lookup                  = Grid_Hierarchy_Lookup<Struct_type,T,d>;
    using Allocator_type                    = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                   = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper                   = Grid_Topology_Helper<Flag_array_mask>;
    using T_Range_Iterator                  = Range_Iterator<d,T_INDEX>;
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;

    struct Voxel
    {
        TV location;
        T density;
        Voxel() {}
        
        Voxel(const TV& location_input,const T density_input)
        {
            for(int v=0;v<d;++v) location[v]=location_input[v];
            if(d==2) location[2]=1;
            density=density_input;
        }
    };
    public: 
    bool interpolated;
    int threads;
    Parse_Args *parse_args;
    std::string output_directory;
    Hierarchy *hierarchy;
    std::vector<Voxel> voxels;
    std::vector<T> node_density;
    unsigned elements_per_block;
    int selected_voxel_level,levels;
    int first_frame, last_frame;
    int step;
    int fix_frame;
    Index_type block_size;
    T threshold;
    T density_factor;
    T Struct_type::* density_channel;
    Grid<T,d> grid;
    int xm,ym,zm;
    T dx,quarter_dx;
  public:
    OpenVDB_Converter()
        :parse_args(nullptr),hierarchy(nullptr),elements_per_block(0),levels(0)
    {}

    ~OpenVDB_Converter()
    {if(hierarchy!=nullptr) delete hierarchy;}

    void Initialize()
    {
        std::istream *input=(d==3)?File_Utilities::Safe_Open_Input(output_directory+"/common/hierarchy.struct3d"):File_Utilities::Safe_Open_Input(output_directory+"/common/hierarchy.struct2d");
        Read_Write<int>::Read(*input,levels);
        Log::cout<<"Levels: "<<levels<<std::endl;
        Read_Write<unsigned>::Read(*input,elements_per_block);
        Log::cout<<"# of elements/block: "<<elements_per_block<<std::endl;
        for(int v=0;v<d;++v) Read_Write<unsigned>::Read(*input,block_size[v]);
        delete input;
        omp_set_num_threads(threads);
        File_Utilities::Read_From_File(output_directory+"/common/fine_grid",grid);
        xm=grid.counts(0); ym=grid.counts(1); zm=grid.counts(2);
        density_channel                         = &Struct_type::ch0;
        dx=grid.dX(0);  quarter_dx=(T).25*dx;
        Log::cout<<"# cells: "<<grid.counts(0)<<"x"<<grid.counts(1)<<"x"<<grid.counts(2)<<std::endl;
        Log::cout<<"dx: "<<dx<<", quarter dx: "<<quarter_dx<<std::endl;
        Log::cout<<"density factor: "<<density_factor<<std::endl;
        Log::cout<<"density threshold: "<<threshold<<std::endl;
        if(interpolated) Log::cout<<"Interpolated"<<std::endl; else Log::cout<<"Raw cell data"<<std::endl;
    }

    int Cell_ID(int i,int j,int k){
        if(i<0) return -1; if(i>=xm) return -1;
        if(j<0) return -1; if(j>=ym) return -1;
        if(k<0) return -1; if(k>=zm) return -1;
        return i*ym*zm+j*zm+k;
    }

    int Node_ID(int i,int j,int k){
        if(i<0) return -1; if(i>=2*xm) return -1;
        if(j<0) return -1; if(j>=2*ym) return -1;
        if(k<0) return -1; if(k>=2*zm) return -1;
        return i*4*ym*zm+j*2*zm+k;
    }

    void Convert_Frame(const int current_frame)
    {
        voxels.clear();
        if(hierarchy!=nullptr) delete hierarchy;
        hierarchy=new Hierarchy(grid,levels);
        
        voxels.resize(grid.Number_Of_Cells().Product());
        
        std::stringstream ss;ss<<output_directory<<"/"<<current_frame;
        std::istream* input1=File_Utilities::Safe_Open_Input(ss.str()+"/flags");
        std::istream* input2=File_Utilities::Safe_Open_Input(ss.str()+"/block_offsets");

        // check for density
        bool draw_density=File_Utilities::File_Exists(ss.str()+"/spgrid_density",false);
        std::istream* input3=draw_density?File_Utilities::Safe_Open_Input(ss.str()+"/spgrid_density"):nullptr;

        for(int level=0;level<levels;++level){unsigned number_of_blocks=0;
            Read_Write<unsigned>::Read(*input2,number_of_blocks);

            for(unsigned block=0;block<number_of_blocks;++block){T_INDEX base_index;
                Read_Write<T_INDEX>::Read(*input2,base_index);
                T_Range_Iterator range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);

                for(unsigned e=0;e<elements_per_block;++e){unsigned flag;T density;
                    const T_INDEX cell_ijk=base_index+range_iterator.Index();
                    Read_Write<unsigned>::Read(*input1,flag); Read_Write<T>::Read(*input3,density);
                    if(flag&(Cell_Type_Interior|Cell_Type_Dirichlet)){const int cell_id=Cell_ID(cell_ijk(0)-1,cell_ijk(1)-1,cell_ijk(2)-1);
                        TV cell_location=grid.Center(cell_ijk); voxels[cell_id]=Voxel(cell_location,density);}
                    range_iterator.Next();}}}
        const T l=(T).04,u=(T)1.; const T gamma_l=(T)1.5,gamma_u=(T)3.;
        if(interpolated){
            Log::cout<<"using interpolated data"<<std::endl;
		    openvdb::initialize();
            openvdb::FloatGrid::Ptr mygrid = openvdb::FloatGrid::create(); openvdb::FloatGrid::Accessor accessor = mygrid->getAccessor();
            mygrid->setTransform(openvdb::math::Transform::createLinearTransform(256./xm*0.1f));
            Log::cout<<"cell width: "<<mygrid->voxelSize()[0]<<std::endl;
#pragma omp parallel for
            for(int cell_id=0;cell_id<voxels.size();++cell_id){
                const int i=cell_id/(ym*zm); const int j=(cell_id-i*ym*zm)/zm; const int k=cell_id-i*ym*zm-j*zm;
                TV cell_location=voxels[cell_id].location; T cell_density=voxels[cell_id].density;
                for(int ii=-1;ii<=1;ii+=2) for(int jj=-1;jj<=1;jj+=2) for(int kk=-1;kk<=1;kk+=2){
                    T interpolated_density=0.; T_INDEX node_ijk({2*i+(ii+1)/2,2*j+(jj+1)/2,2*k+(kk+1)/2}); 
                    int node_index=Node_ID(node_ijk(0),node_ijk(1),node_ijk(2)); assert(node_index!=-1);
                    TV node_location({cell_location(0)+ii*quarter_dx,cell_location(1)+jj*quarter_dx,cell_location(2)+kk*quarter_dx});
                    for(T_Range_Iterator iterator(T_INDEX(),T_INDEX({ii,jj,kk}));iterator.Valid();iterator.Next()){
                        T_INDEX iter_cell_index=T_INDEX({i+1,j+1,k+1})+iterator.Index();   
                        int iter_id=Cell_ID(iter_cell_index(0)-1,iter_cell_index(1)-1,iter_cell_index(2)-1);
                        if(iter_id==-1) ;
                        else{TV iter_cell_location=grid.Center(iter_cell_index);
                        T iter_density=voxels[iter_id].density; T factor=1.;
                        for(int axis=0;axis<3;++axis){ T delta_dis=1.-abs(iter_cell_location(axis)-node_location(axis))/dx; factor*=delta_dis;}
                        interpolated_density+=factor*iter_density;}}
                        interpolated_density*=density_factor;
                        // transform: 
                        interpolated_density=(clamp(interpolated_density,l,u)-l)/(u-l);
                        interpolated_density=pow(interpolated_density,gamma_l)*((T)1.-interpolated_density)
                                                +pow(interpolated_density,gamma_u)*interpolated_density;
                        if(interpolated_density>threshold){ openvdb::Coord xyz(node_ijk(0),node_ijk(1),node_ijk(2));
                            accessor.setValue(xyz,interpolated_density);}}
            }

            string output_filename=output_directory+"/converted_data/"+std::to_string(current_frame/step)+".vdb";
            mygrid->setName("density");
            openvdb::io::File(output_filename).write({mygrid});}
        else{
            Log::cout<<"using cell data"<<std::endl;
            openvdb::initialize();
            openvdb::FloatGrid::Ptr mygrid = openvdb::FloatGrid::create(); openvdb::FloatGrid::Accessor accessor = mygrid->getAccessor();
            mygrid->setTransform(openvdb::math::Transform::createLinearTransform(512./xm*0.1f));
            Log::cout<<"cell width: "<<mygrid->voxelSize()[0]<<std::endl;
            for(int i=0;i<xm;++i) for(int j=0;j<ym;++j) for(int k=0;k<zm;++k){
                int cell_id=Cell_ID(i,j,k); T cell_density=voxels[cell_id].density;
                // transform: 
                cell_density*=density_factor;
                cell_density=(clamp(cell_density,l,u)-l)/(u-l);
                cell_density=pow(cell_density,gamma_l)*((T)1.-cell_density)+pow(cell_density,gamma_u)*cell_density;
                if(cell_density>threshold){ openvdb::Coord xyz(i,j,k);
                    accessor.setValue(xyz,cell_density);}}
            string output_filename=output_directory+"/converted_data/"+std::to_string(current_frame/step)+".vdb";
            mygrid->setName("density");
            openvdb::io::File(output_filename).write({mygrid});}

        delete input2;delete input1;
        if(input3!=nullptr) delete input3;
    }

    void Register_Options()
    {
        if(!parse_args) return;
        parse_args->Add_Option_Argument("-inter","Use interpolation");
        parse_args->Add_String_Argument("-o","","output directory");
        parse_args->Add_Integer_Argument("-first_frame",0,"first frame");
        parse_args->Add_Integer_Argument("-last_frame",0,"last frame");
        parse_args->Add_Integer_Argument("-fix",0,"fix frame");
        parse_args->Add_Integer_Argument("-step",1,"step");
        parse_args->Add_Integer_Argument("-threads",1,"threads");
        parse_args->Add_Double_Argument("-df",4.,"density factor");
        parse_args->Add_Double_Argument("-th",0.,"density threshold");
    }

    void Parse_Options()
    {
        if(!parse_args) return;
        fix_frame=parse_args->Get_Integer_Value("-fix");
        interpolated=parse_args->Get_Option_Value("-inter");
        output_directory=parse_args->Get_String_Value("-o");
        if(fix_frame>0){last_frame=fix_frame; first_frame=fix_frame; step=1;}
        else{ first_frame=parse_args->Get_Integer_Value("-first_frame");
        last_frame=parse_args->Get_Integer_Value("-last_frame");
        step=parse_args->Get_Integer_Value("-step");}
        
        density_factor=parse_args->Get_Double_Value("-df");
        threshold=parse_args->Get_Double_Value("-th");
        threads=parse_args->Get_Integer_Value("-threads");
        if(threads>8) threads=8;
    }

    void Parse(int argc,char* argv[])
    {
        parse_args=new Parse_Args;
        Register_Options();    
        parse_args->Parse(argc,argv);
        Parse_Options();
    }
    


        




};
}
#endif

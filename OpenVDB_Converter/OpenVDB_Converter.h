//!#####################################################################
//! \file OpenVDB_Converter.h
//!#####################################################################
// Class OpenVDB_Converter
//######################################################################
#ifndef __OpenVDB_Converter__
#define __OpenVDB_Converter__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Dynamics/Rigid_Bodies/Rigid_Body.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include <assert.h>
// #include <openvdb/openvdb.h>
#include "Viewer_Data.h"
#include "../Poisson_Solver/Grid_Hierarchy_Projection.h"
#include <iostream>
using namespace std;

namespace Nova{
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
    const std::string directory_name;
    Hierarchy *hierarchy;
    std::vector<Voxel> voxels;
    std::vector<T> node_density;
    unsigned elements_per_block;
    int selected_voxel_level,levels;
    Index_type block_size;

    Grid<T,d> grid;
    int xm,ym,zm;
    T Struct_type::* density_channel;
    
  public:
    OpenVDB_Converter(const std::string& directory_name_input)
        :hierarchy(nullptr),directory_name(directory_name_input),elements_per_block(0),levels(0)
    {
        std::istream *input=(d==3)?File_Utilities::Safe_Open_Input(directory_name+"/common/hierarchy.struct3d"):File_Utilities::Safe_Open_Input(directory_name+"/common/hierarchy.struct2d");
        Read_Write<int>::Read(*input,levels);
        std::cout<<"Levels: "<<levels<<std::endl;
        Read_Write<unsigned>::Read(*input,elements_per_block);
        std::cout<<"# of elements/block: "<<elements_per_block<<std::endl;
        for(int v=0;v<d;++v) Read_Write<unsigned>::Read(*input,block_size[v]);
        delete input;

        File_Utilities::Read_From_File(directory_name+"/common/fine_grid",grid);
        xm=grid.counts(0); ym=grid.counts(1); zm=grid.counts(2);
        std::cout<<"size: "<<xm<<","<<ym<<","<<zm<<std::endl;
        density_channel                         = &Struct_type::ch0;
    }

    ~OpenVDB_Converter()
    {if(hierarchy!=nullptr) delete hierarchy;}

    int id(int i,int j,int k){
        if(i<0) return -1; if(i>=xm) return -1;
        if(j<0) return -1; if(j>=ym) return -1;
        if(k<0) return -1; if(k>=zm) return -1;
        return i*ym*zm+j*zm+k;
    }

    int node_id(int i,int j,int k){
        if(i<0) return -1; if(i>=2*xm) return -1;
        if(j<0) return -1; if(j>=2*ym) return -1;
        if(k<0) return -1; if(k>=2*zm) return -1;
        return i*2*ym*2*zm+j*2*zm+k;
    }

    void Clear_Buffers()
    {
        voxels.clear();
    }

    void Read_From_Frame(const int current_frame)
    {
        Clear_Buffers();
        
        if(hierarchy!=nullptr) delete hierarchy;
        hierarchy=new Hierarchy(grid,levels);
        const T dx=grid.dX(0); const T one_over_dx=grid.one_over_dX(0); const T half_dx=(T).5*dx;
        voxels.resize(grid.Number_Of_Cells().Product());
        node_density.resize(8*voxels.size());
        std::cout<<"# cells: "<<voxels.size()<<", # nodes: "<<node_density.size()<<std::endl;
        std::cout<<"dx: "<<dx<<", one_over_dx: "<<one_over_dx<<", half dx: "<<half_dx<<std::endl;
        std::stringstream ss;ss<<directory_name<<"/"<<current_frame;
        std::istream* input1=File_Utilities::Safe_Open_Input(ss.str()+"/flags");
        std::istream* input2=File_Utilities::Safe_Open_Input(ss.str()+"/block_offsets");

        // check for density
        bool draw_density=File_Utilities::File_Exists(ss.str()+"/spgrid_density",false);
        std::istream* input3=draw_density?File_Utilities::Safe_Open_Input(ss.str()+"/spgrid_density"):nullptr;


        for(int level=0;level<levels;++level){unsigned number_of_blocks=0;
            Read_Write<unsigned>::Read(*input2,number_of_blocks);

            for(unsigned block=0;block<number_of_blocks;++block){T_INDEX base_index;
                Read_Write<T_INDEX>::Read(*input2,base_index);
                T_Range_Iterator range_iterator(T_INDEX({0,0,0}),*reinterpret_cast<T_INDEX*>(&block_size)-1);

                for(unsigned e=0;e<elements_per_block;++e){unsigned flag;
                    // 1-based
                    const T_INDEX cell_ijk=base_index+range_iterator.Index();
                    const int cell_id=id(cell_ijk(0)-1,cell_ijk(1)-1,cell_ijk(2)-1);                    
                    Read_Write<unsigned>::Read(*input1,flag);
                    TV cell_location=grid.Center(cell_ijk);
                    if(draw_density){T density;
                        Read_Write<T>::Read(*input3,density);
                        if(flag&(Cell_Type_Interior|Cell_Type_Dirichlet)){
                            voxels[cell_id]=Voxel(cell_location,density);
                        }
                    }
                    range_iterator.Next();}}}

		// openvdb::initialize();
        string output_filename=directory_name+"/nodes_data/converted_"+std::to_string(current_frame)+".txt";
        
        // openvdb::FloatGrid::Ptr mygrid = openvdb::FloatGrid::create();
        // openvdb::FloatGrid::Accessor accessor = mygrid->getAccessor();
        for(int i=0;i<xm;++i) for(int j=0;j<ym;++j) for(int k=0;k<zm;++k){
            // 0-based
            int cell_id=id(i,j,k);
            TV cell_location=voxels[cell_id].location;
            double cell_density=voxels[cell_id].density;
            // if(i==0 || j==0) std::cout<<"**********************************"<<std::endl;
            // if(i==0 || j==0) std::cout<<"cell index: "<<i<<","<<j<<","<<k<<std::endl;
            for(int ii=-1;ii<=1;ii+=2) for(int jj=-1;jj<=1;jj+=2) for(int kk=-1;kk<=1;kk+=2){
                        T interpolated_density=0.;
                        T_INDEX node_ijk({2*i+(ii+1)/2,2*j+(jj+1)/2,2*k+(kk+1)/2});
                        // if(node_ijk(0)==0) std::cout<<"node index: "<<node_ijk(0)<<","
                           //     <<node_ijk(1)<<","<<node_ijk(2)<<std::endl;
                        int node_index=node_id(node_ijk(0),node_ijk(1),node_ijk(2));
                        assert(node_index!=-1);
                        TV node_location({cell_location(0)+ii*half_dx,cell_location(1)+jj*half_dx,cell_location(2)+kk*half_dx});
                        for(T_Range_Iterator iterator(T_INDEX({0,0,0}),T_INDEX({ii,jj,kk}));iterator.Valid();iterator.Next()){
                            // 0-based
                            T_INDEX iter_cell_index=T_INDEX({i,j,k})+iterator.Index();   
                            int iter_id=id(iter_cell_index(0),iter_cell_index(1),iter_cell_index(2));
                            if(iter_id==-1) ;
                            else{TV iter_cell_location=grid.Center(iter_cell_index+T_INDEX{1,1,1});
                            T iter_density=voxels[iter_id].density;
                            // if(node_ijk(0)==0 || node_ijk(1)==0) std::cout<<"iter cell: "<<iter_cell_index(0)
                                // <<","<<iter_cell_index(1)<<","<<iter_cell_index(2)<<": "<<iter_density<<std::endl;
                            T factor=1.;
                            for(int axis=0;axis<3;++axis){
                            T delta_dis=1.-abs(iter_cell_location(axis)-node_location(axis))/dx;
                            factor*=delta_dis;}
                    interpolated_density+=factor*iter_density;}}
                    node_density[node_index]=interpolated_density;
                    // if(node_ijk(0)==0 || node_ijk(1)==0 || 
                        // node_ijk(0)==2*xm-1 || node_ijk(1)==2*ym-1) 
                        // if(interpolated_density>(T)1e-5) 
                        // std::cout<<"index: "<<node_ijk(0)<<","<<node_ijk(1)<<","
                        // <<node_ijk(2)<<", node density: "<<node_density[node_index]<<std::endl;
                // openvdb::Coord xyz(node_ijk(0),node_ijk(1),node_ijk(2));
                // accessor.setValue(xyz, float(interpolated_density));
                }
                
            // std::cout<<"**********************************"<<std::endl;
                
                }
        
        FILE* fp = fopen(output_filename.c_str(), "w");
        for(int count=0;count<node_density.size();++count)
            fprintf(fp, "%.6f\n",node_density[count]);
        fclose(fp);
        // mygrid->setName("density");
        // openvdb::io::File(output_filename).write({mygrid});
        delete input2;delete input1;
        if(input3!=nullptr) delete input3;
    }



        




};
}
#endif

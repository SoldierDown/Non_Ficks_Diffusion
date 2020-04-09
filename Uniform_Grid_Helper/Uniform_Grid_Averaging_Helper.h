//!#####################################################################
//! \file Uniform_Grid_Hierarchy_Averaging.h
//!#####################################################################
// Class Uniform_Grid_Hierarchy_Averaging
//######################################################################
#ifndef __Uniform_Grid_Hierarchy_Averaging__
#define __Uniform_Grid_Hierarchy_Averaging__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/SPGrid/Tools/SPGrid_Constrain_T_Junction_Nodes.h>
#include <nova/SPGrid/Tools/SPGrid_Downsample_Accumulate_Shared.h>
#include <nova/SPGrid/Tools/SPGrid_Masked_Average_Offset_Grid.h>
#include <nova/SPGrid/Tools/SPGrid_Masked_Normalize.h>
#include <nova/SPGrid/Tools/SPGrid_Upsample_Inject_Shared.h>
#include <nova/Tools/Log/Log.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
namespace Nova{

extern int number_of_threads;

template<class Struct_type,class T,int d>
class Uniform_Grid_Averaging_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    enum {nodes_per_cell            = Topology_Helper::number_of_nodes_per_cell};
    enum {nodes_per_face            = Topology_Helper::number_of_nodes_per_face};

  public:
    // Checked
    static void Uniform_Grid_Average_Face_Velocities_To_Cells(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                                            Channel_Vector& face_channels,Channel_Vector& cell_channels,const Vector<uint64_t,d>& other_face_offsets)
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto uniform_grid_average_face_velocities_to_cells_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&Cell_Type_Interior){
                for(int axis=0;axis<d;++axis){unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
                    auto face_data=allocator.template Get_Const_Array<Struct_type,T>(face_channels(axis));
                    auto cell_data=allocator.template Get_Array<Struct_type,T>(cell_channels(axis));
                    uint64_t other_offset=Flag_array_mask::Packed_Add(offset,other_face_offsets(axis));
                    assert(flags(offset)&face_valid_mask); assert(flags(other_face_offsets)&face_valid_mask);
                    cell_data(offset)=(T).5*(face_data(offset)+face_data(other_offset));}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,uniform_grid_average_face_velocities_to_cells_helper);
    }
    static void Uniform_Grid_Average_Face_Velocities_To_Faces(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                                                Channel_Vector& face_velocity_channels,Channel_Vector& interpolated_face_velocity_channels,
                                                                const Vector<uint64_t,d>& negative_face_offsets,uint64_t* nodes_of_cell_offsets,const int axis)
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto uniform_grid_average_face_velocities_to_faces_helper=[&](uint64_t offset)
        {
            unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&face_active_mask){
                for(int v=0;v<d;++v){
                    auto face_velocity=allocator.template Get_Const_Array<Struct_type,T>(face_velocity_channels(v));
                    auto interpolated_face_velocity=allocator.template Get_Array<Struct_type,T>(interpolated_face_velocity_channels(v));
                    if(v==axis) interpolated_face_velocity(offset)=face_velocity(offset);
                    else{uint64_t base_offset=Flag_array_mask::Packed_Add(offset,negative_face_offsets(axis)); T axis_velocity=(T)0.;
                    for(int face=0;face<nodes_per_cell;++face){ 
                        uint64_t face_offset=Flag_array_mask::Packed_Add(base_offset,nodes_of_cell_offsets[face]);
                        axis_velocity+=face_velocity(face_offset);}
                    interpolated_face_velocity(offset)=(T).25*axis_velocity;}}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,uniform_grid_average_face_velocities_to_faces_helper);
    }
};
}
#endif

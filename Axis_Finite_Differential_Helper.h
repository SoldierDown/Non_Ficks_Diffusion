//!#####################################################################
//! \file Axis_Finite_Differential_Helper.h
//!#####################################################################
// Class Axis_Finite_Differential_Helper
//######################################################################
#ifndef __Axis_Finite_Differential_Helper__
#define __Axis_Finite_Differential_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Axis_Finite_Differential_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Axis_Finite_Differential_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* cell_channel,T Struct_type::* result_channel,
                                    const T one_over_2dx,const int axis)
    {Run(allocator,blocks,cell_channel,result_channel,one_over_2dx,axis);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* cell_channel,T Struct_type::* result_channel,
            const T one_over_2dx,const int axis) const
    {
        Log::Scope scope("Axis_Finite_Differential_Helper");
        uint64_t face_neighbor_offsets[number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        uint64_t axis_neighbor_offsets[2];
        axis_neighbor_offsets[0]=face_neighbor_offsets[2*axis]; axis_neighbor_offsets[1]=face_neighbor_offsets[2*axis+1];
        auto cell_data=allocator.template Get_Const_Array<Struct_type,T>(cell_channel);
        auto result_data=allocator.template Get_Array<Struct_type,T>(result_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto axis_finite_differential_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&Cell_Type_Interior) 
            {uint64_t neighbor_offset[2];for(int nb=0;nb<2;nb++) neighbor_offset[nb]=Flag_array_mask::Packed_Add(offset,axis_neighbor_offsets[nb]);
                result_data(offset)=one_over_2dx*(cell_data(neighbor_offset[1])-cell_data(neighbor_offset[0]));}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,axis_finite_differential_helper);
    }

};
}
#endif

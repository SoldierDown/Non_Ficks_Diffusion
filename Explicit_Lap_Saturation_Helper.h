//!#####################################################################
//! \file Explicit_Lap_Saturation_Helper.h
//!#####################################################################
// Class Explicit_Lap_Saturation_Helper
//######################################################################
#ifndef __Explicit_Lap_Saturation_Helper__
#define __Explicit_Lap_Saturation_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
// #include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Explicit_Lap_Saturation_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Explicit_Lap_Saturation_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* saturation_channel,T Struct_type::* lap_saturation_channel,const T one_over_dx2)
    {Run(allocator,blocks,saturation_channel,lap_saturation_channel,one_over_dx2);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* saturation_channel,T Struct_type::* lap_saturation_channel,const T one_over_dx2) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto saturation=allocator.template Get_Const_Array<Struct_type,T>(saturation_channel); auto lap_saturation=allocator.template Get_Array<Struct_type,T>(lap_saturation_channel);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        auto explicit_lap_saturation_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){
                    for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(flags(neighbor_offset)&Cell_Type_Interior || flags(offset)&Cell_Type_Dirichlet) lap_saturation(offset)+=one_over_dx2*(saturation(neighbor_offset)-saturation(offset));}}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,explicit_lap_saturation_helper);
    }
};
}
#endif

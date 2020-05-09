//!#####################################################################
//! \file Ficks_RHS_Helper.h
//!#####################################################################
// Class Ficks_RHS_Helper
//######################################################################
#ifndef __Ficks_RHS_Helper__
#define __Ficks_RHS_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
// #include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Ficks_RHS_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Ficks_RHS_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                  T Struct_type::* saturation_channel,T Struct_type::* rhs_channel,const T a)
    {Run(allocator,blocks,saturation_channel,rhs_channel,a);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* saturation_channel,T Struct_type::* rhs_channel,const T a) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto saturation=allocator.template Get_Const_Array<Struct_type,T>(saturation_channel); auto rhs=allocator.template Get_Array<Struct_type,T>(rhs_channel);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        auto ficks_rhs_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){ rhs(offset)=saturation(offset);
                    for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(flags(neighbor_offset)&Cell_Type_Dirichlet) rhs(offset)+=a*saturation(neighbor_offset);}}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,ficks_rhs_helper);
    }

};
}
#endif

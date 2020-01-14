//!#####################################################################
//! \file Non_Ficks_RHS_Helper.h
//!#####################################################################
// Class Non_Ficks_RHS_Helper
//######################################################################
#ifndef __Non_Ficks_RHS_Helper__
#define __Non_Ficks_RHS_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Vectors/Vector.h>
#include "MPM_Flags.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Non_Ficks_RHS_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Non_Ficks_RHS_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                  T Struct_type::* saturation_channel,T Struct_type::* div_Qc_channel,T Struct_type::*& rhs_channel,unsigned Struct_type::* flags_channel,const T coeff1,const T coeff2)
    {Run(allocator,blocks,saturation_channel,div_Qc_channel,rhs_channel,flags_channel,coeff1,coeff2);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* saturation_channel,T Struct_type::* div_Qc_channel,T Struct_type::*& rhs_channel,unsigned Struct_type::* flags_channel,const T coeff1,const T coeff2) const
    {
        auto saturation=allocator.template Get_Const_Array<Struct_type,T>(saturation_channel);
        auto div_Qc=allocator.template Get_Const_Array<Struct_type,T>(div_Qc_channel);
        auto rhs=allocator.template Get_Array<Struct_type,T>(rhs_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(flags_channel);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        auto non_ficks_rhs_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Node_Saturated){
                    rhs(offset)=saturation(offset)-coeff2*div_Qc(offset);
                    for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if((flags(neighbor_offset)&Node_Active)&&(!(flags(neighbor_offset)&Node_Saturated))) 
                            rhs(offset)+=coeff1*saturation(neighbor_offset);}}}
        };

        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            non_ficks_rhs_helper(offset);}
    }

};
}
#endif

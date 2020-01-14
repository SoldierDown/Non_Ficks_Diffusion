//!#####################################################################
//! \file Diffusion_Multiply_Helper.h
//!#####################################################################
// Class Diffusion_Multiply_Helper
//######################################################################
#ifndef __Diffusion_Multiply_Helper__
#define __Diffusion_Multiply_Helper__

#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Diffusion_Multiply_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Diffusion_Multiply_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* x_channel,T Struct_type::* result_channel,
                            const bool FICKS,const T a,const T four_a_plus_one,const T coeff1)
    {Run(allocator,blocks,x_channel,result_channel,FICKS,a,four_a_plus_one,coeff1);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* x_channel,T Struct_type::* result_channel,
                const bool FICKS,const T a,const T four_a_plus_one,const T coeff1) const
    {
        if(FICKS){
        auto x=allocator.template Get_Const_Array<Struct_type,T>(x_channel);
        auto result=allocator.template Get_Array<Struct_type,T>(result_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);

        auto ficks_diffusion_multiply_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Node_Saturated){ result(offset)=four_a_plus_one*x(offset);
                    for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(flags(neighbor_offset)&Node_Active) result(offset)-=a*x(neighbor_offset); }}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,ficks_diffusion_multiply_helper);}
        else{
        auto x=allocator.template Get_Const_Array<Struct_type,T>(x_channel);
        auto result=allocator.template Get_Array<Struct_type,T>(result_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);

        auto non_ficks_diffusion_multiply_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Node_Saturated){ result(offset)=(1.+4.*coeff1)*x(offset);
                    for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(flags(neighbor_offset)&Node_Active) result(offset)-=coeff1*x(neighbor_offset); }}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,non_ficks_diffusion_multiply_helper);}
    }
};
}
#endif

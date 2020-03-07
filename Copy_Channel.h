//!#####################################################################
//! \file Copy_Channel.h
//!#####################################################################
// Class Copy_Channel
//######################################################################
#ifndef __Copy_Channel__
#define __Copy_Channel__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Base_struct_type,class Multigrid_struct_type,class T,int d>
class Copy_Channel
{
    using Base_flags_type                       = typename Base_struct_type::Flags_type;
    using Hierarchy_Multigrid                   = Grid_Hierarchy<Multigrid_struct_type,T,d>;
    using Base_allocator_type                   = SPGrid::SPGrid_Allocator<Base_struct_type,d>;
    using Base_flag_array_mask                  = typename Base_allocator_type::template Array_mask<unsigned>;
    using Multigrid_allocator_type              = SPGrid::SPGrid_Allocator<Multigrid_struct_type,d>;
    using Multigrid_flag_array_mask             = typename Multigrid_allocator_type::template Array_mask<unsigned>;

  public:
    Copy_Channel(Base_allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Hierarchy_Multigrid& hierarchy,
                 T Base_struct_type::* channel_base,T Multigrid_struct_type::* channel_multigrid,const int level,const bool copy_to=true)
    {Run(allocator,blocks,hierarchy,channel_base,channel_multigrid,level,copy_to);}

    void Run(Base_allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Hierarchy_Multigrid& hierarchy,
             T Base_struct_type::* channel_base,T Multigrid_struct_type::* channel_multigrid,const int level,const bool copy_to) const
    {
        auto data_base=allocator.template Get_Array<Base_struct_type,T>(channel_base);
        auto data_multigrid=hierarchy.Allocator(level).template Get_Array<Multigrid_struct_type,T>(channel_multigrid);

        auto copy_channel=[&](uint64_t offset)
        {
            if(copy_to) for(int e=0;e<Base_flag_array_mask::elements_per_block;++e,offset+=sizeof(Base_flags_type)){
                const uint64_t translated_offset=Multigrid_flag_array_mask::template Translate_Linear_Offset<Base_flag_array_mask>(offset);
                data_multigrid(translated_offset)=data_base(offset);}
            else for(int e=0;e<Base_flag_array_mask::elements_per_block;++e,offset+=sizeof(Base_flags_type)){
                const uint64_t translated_offset=Multigrid_flag_array_mask::template Translate_Linear_Offset<Base_flag_array_mask>(offset);
                data_base(offset)=data_multigrid(translated_offset);}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,copy_channel);
    }
};
}
#endif

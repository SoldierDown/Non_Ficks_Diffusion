//!#####################################################################
//! \file Traverse_Helper.h
//!#####################################################################
// Class Traverse_Helper
//######################################################################
#ifndef __Traverse_Helper__
#define __Traverse_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include "MPM_Flags.h"
namespace Nova{
template<class Struct_type,class T,int d>
class Traverse_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Traverse_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel,unsigned Struct_type::* flags_channel)
    {Run(allocator,blocks,channel,flags_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::*channel,unsigned Struct_type::* flags_channel) const
    {
        auto c=allocator.template Get_Const_Array<Struct_type,T>(channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(flags_channel);
        auto traverse_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Node_Saturated) Log::cout<<c(offset)<<std::endl;
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            traverse_helper(offset);}
    }
};
}
#endif

//!#####################################################################
//! \file Compare_Helper.h
//!#####################################################################
// Class Compare_Helper
//######################################################################
#ifndef __Compare_Helper__
#define __Compare_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include "MPM_Flags.h"
namespace Nova{
template<class Struct_type,class T,int d>
class Compare_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Compare_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel1,T Struct_type::* channel2,unsigned Struct_type::* flags_channel)
    {Run(allocator,blocks,channel1,channel2,flags_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::*channel1,T Struct_type::* channel2,unsigned Struct_type::* flags_channel) const
    {
        auto c1=allocator.template Get_Const_Array<Struct_type,T>(channel1);
        auto c2=allocator.template Get_Const_Array<Struct_type,T>(channel2);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(flags_channel);
        T l2_norm;
        auto compare_helper=[&](uint64_t offset, T& l2_norm)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Node_Saturated) l2_norm+=fabs(Nova_Utilities::Sqr(c1(offset)-c2(offset)));
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            compare_helper(offset,l2_norm);}
        Log::cout<<"difference: "<<std::sqrt(l2_norm)<<std::endl;
    }
};
}
#endif

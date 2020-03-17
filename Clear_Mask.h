//!#####################################################################
//! \file Clear_Mask.h
//!#####################################################################
// Class Clear_Mask
//######################################################################
#ifndef __Clear_Mask__
#define __Clear_Mask__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
 
namespace Nova{
template<class Struct_type,class T,int d>
class Clear_Mask
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Clear_Mask(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const unsigned mask)
    {Run(allocator,blocks,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const unsigned mask) const
    {
        auto flags=allocator.template Get_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto clear_mask=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                flags(offset)&=~mask;
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,clear_mask);
    }
};
}
#endif

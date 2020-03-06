//!#####################################################################
//! \file Clear_Non_Active.h
//!#####################################################################
// Class Clear_Non_Active
//######################################################################
#ifndef __Clear_Non_Active__
#define __Clear_Non_Active__

#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Clear_Non_Active
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Clear_Non_Active(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel)
    {Run(allocator,blocks,channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel) const
    {
        auto data=allocator.template Get_Array<Struct_type,T>(channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto clear_non_active=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(!(flags(offset)&Cell_Type_Interior)) data(offset)=(T)0.;
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,clear_non_active);
    }
};
}
#endif

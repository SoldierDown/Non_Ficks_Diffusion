//!#####################################################################
//! \file Saturation_Clamp_Heler.h
//!#####################################################################
// Class Saturation_Clamp_Heler
//######################################################################
#ifndef __Saturation_Clamp_Heler__
#define __Saturation_Clamp_Heler__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include "MPM_Flags.h"
namespace Nova{
template<class Struct_type,class T,int d>
class Saturation_Clamp_Heler
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Saturation_Clamp_Heler(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* new_saturation_channel,T Struct_type::* saturation_channel)
    {Run(allocator,blocks,new_saturation_channel,saturation_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* new_saturation_channel,T Struct_type::* saturation_channel) const
    {
        auto saturation=allocator.template Get_Array<Struct_type,T>(saturation_channel);
        auto new_saturation=allocator.template Get_Const_Array<Struct_type,T>(new_saturation_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto saturation_clamp_heler=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Node_Saturated) saturation(offset)=Nova_Utilities::Clamp(new_saturation(offset),(T)0.,(T)1.);
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,saturation_clamp_heler);
    }
};
}
#endif

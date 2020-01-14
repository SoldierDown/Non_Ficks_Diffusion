//!#####################################################################
//! \file Diffusion_Clear_Non_Active_Helper.h
//!#####################################################################
// Class Diffusion_Clear_Non_Active_Helper
//######################################################################
#ifndef __Diffusion_Clear_Non_Active_Helper__
#define __Diffusion_Clear_Non_Active_Helper__

#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include "../MPM_Flags.h"
namespace Nova{
template<class Struct_type,class T,int d>
class Diffusion_Clear_Non_Active_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Diffusion_Clear_Non_Active_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel)
    {Run(allocator,blocks,channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel) const
    {
        auto data=allocator.template Get_Array<Struct_type,T>(channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto diffusion_clear_non_active_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(!(flags(offset)&Node_Saturated)) data(offset)=(T)1.;
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,diffusion_clear_non_active_helper);
    }
};
}
#endif

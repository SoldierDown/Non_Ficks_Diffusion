//!#####################################################################
//! \file Multigrid_RHS_Helper.h
//!#####################################################################
// Class Multigrid_RHS_Helper
//######################################################################
#ifndef __Multigrid_RHS_Helper__
#define __Multigrid_RHS_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include "../MPM_Flags.h"


namespace Nova{
template<class Struct_type,class T,int d>
class Multigrid_RHS_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Multigrid_RHS_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* saturation_channel,T Struct_type::* rhs_channel)
    {Run(allocator,blocks,saturation_channel,rhs_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* saturation_channel,T Struct_type::* rhs_channel) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto saturation=allocator.template Get_Const_Array<Struct_type,T>(saturation_channel);
        auto rhs=allocator.template Get_Array<Struct_type,T>(rhs_channel);

        auto multigrid_rhs_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Saturated) rhs(offset)=saturation(offset);
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,multigrid_rhs_helper);        
    }

};
}
#endif

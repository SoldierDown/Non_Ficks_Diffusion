//!#####################################################################
//! \file Multiply_Helper.h
//!#####################################################################
// Class Multiply_Helper
//######################################################################
#ifndef __Multiply_Helper__
#define __Multiply_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include "MPM_Flags.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Multiply_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;

  public:
    Multiply_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector x_channel,
                         Channel_Vector& result_channel,const T coeff)
    {Run(allocator,blocks,x_channel,result_channel,coeff);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector x_channel,
             Channel_Vector& result_channel,const T coeff) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(flags_channel);
        auto multiply_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) 
              for(int v=0;v<d;v++) if(flags(offset)|Node_Saturated)
                // allocator.template Get_Array<Struct_type,T>(result_channel(v))(offset)=


        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,multiply_helper);
    }
};
}
#endif

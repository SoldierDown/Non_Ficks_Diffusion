//!#####################################################################
//! \file Multiply_Helper.h
//!#####################################################################
// Class Multiply_Helper
//######################################################################
#ifndef __Multiply_Helper__
#define __Multiply_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include "../MPM_Flags.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Multiply_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;

  public:
    Multiply_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector x_channels,Channel_Vector& result_channels,
                      const T scaled_dt_squared,const unsigned mask)
    {Run(allocator,blocks,x_channels,result_channels,scaled_dt_squared,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector x_channels,
             Channel_Vector& result_channels,const T scaled_dt_squared,const unsigned mask) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto mass=allocator.template  Get_Const_Array<Struct_type,T>(&Struct_type::ch0);
        auto multiply_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) 
              if(flags(offset)&mask) for(int v=0;v<d;v++) 
                allocator.template Get_Array<Struct_type,T>(result_channels(v))(offset)=scaled_dt_squared/mass(offset)*allocator.template Get_Array<Struct_type,T>(result_channels(v))(offset)
                                                              +allocator.template Get_Const_Array<Struct_type,T>(x_channels(v))(offset);
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,multiply_helper);
    }
};
}
#endif

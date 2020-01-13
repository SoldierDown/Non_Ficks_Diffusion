//!#####################################################################
//! \file Multiply_Helper.h
//!#####################################################################
// Class Multiply_Helper
//######################################################################
#ifndef __Multiply_Helper__
#define __Multiply_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>

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
                         Channel_Vector& result_channel,const unsigned mask)
    {Run(allocator,blocks,x_channel, result_channel,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector x_channel,
             Channel_Vector& result_channel,const unsigned mask) const
    {
//         auto data1=allocator.template Get_Const_Array<Struct_type,T>(x_channel);
//         auto data2=allocator.template Get_Array<Struct_type,T>(result_channel);
//         auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
//         double temp_result=0;

// #pragma omp parallel for reduction(+:temp_result)
//         for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
//             for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
//                 if(flags(offset)&mask) temp_result+=data1(offset)*data2(offset);}

//         result+=temp_result;
    
    }
};
}
#endif

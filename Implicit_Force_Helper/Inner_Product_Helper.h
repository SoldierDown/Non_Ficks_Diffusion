//!#####################################################################
//! \file Inner_Product_Helper.h
//!#####################################################################
// Class Inner_Product_Helper
//######################################################################
#ifndef __Inner_Product_Helper__
#define __Inner_Product_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Inner_Product_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Inner_Product_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel1,
                         T Struct_type::* channel2,double& result,const unsigned mask)
    {Run(allocator,blocks,channel1,channel2,result,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel1,
             T Struct_type::* channel2,double& result,const unsigned mask) const
    {
        auto data1=allocator.template Get_Const_Array<Struct_type,T>(channel1);
        auto data2=allocator.template Get_Const_Array<Struct_type,T>(channel2);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        double temp_result=0;

#pragma omp parallel for reduction(+:temp_result)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) temp_result+=data1(offset)*data2(offset);}

        result+=temp_result;
    }
};
}
#endif

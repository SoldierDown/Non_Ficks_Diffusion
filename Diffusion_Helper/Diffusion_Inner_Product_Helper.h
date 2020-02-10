//!#####################################################################
//! \file Diffusion_Inner_Product_Helper.h
//!#####################################################################
// Class Diffusion_Inner_Product_Helper
//######################################################################
#ifndef __Diffusion_Inner_Product_Helper__
#define __Diffusion_Inner_Product_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Diffusion_Inner_Product_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Diffusion_Inner_Product_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel0,
                         T Struct_type::* channel1,T& result,const unsigned mask)
    {Run(allocator,blocks,channel0,channel1,result,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel0,
             T Struct_type::* channel1,T& result,const unsigned mask) const
    {
        result=(T)0.;
        auto c0=allocator.template Get_Const_Array<Struct_type,T>(channel0); auto c1=allocator.template Get_Const_Array<Struct_type,T>(channel1);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        T temp_result=(T)0.;

#pragma omp parallel for reduction(+:temp_result)
        for(int b=0;b<blocks.second;b++){
            uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) temp_result+=c0(offset)*c1(offset);}
        result+=temp_result;
    }
};
}
#endif

//!#####################################################################
//! \file Diffusion_Convergence_Norm_Helper.h
//!#####################################################################
// Class Diffusion_Convergence_Norm_Helper
//######################################################################
#ifndef __Diffusion_Convergence_Norm_Helper__
#define __Diffusion_Convergence_Norm_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <algorithm>

namespace Nova{
template<class Struct_type,class T,int d>
class Diffusion_Convergence_Norm_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Diffusion_Convergence_Norm_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            T Struct_type::* channel,T& result,const unsigned mask)
    {Run(allocator,blocks,channel,result,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             T Struct_type::* channel,T& result,const unsigned mask) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto data=allocator.template Get_Const_Array<Struct_type,T>(channel); 
        result=0;
        T max_value=0;

#pragma omp parallel for reduction(max:max_value)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) max_value=std::max(max_value,Nova_Utilities::Sqr(data(offset)));}
        result=std::sqrt(std::max(result,max_value));
    }


};
}
#endif

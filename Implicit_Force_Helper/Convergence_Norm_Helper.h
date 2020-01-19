//!#####################################################################
//! \file Convergence_Norm_Helper.h
//!#####################################################################
// Class Convergence_Norm_Helper
//######################################################################
#ifndef __Convergence_Norm_Helper__
#define __Convergence_Norm_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <algorithm>

namespace Nova{
template<class Struct_type,class T,int d>
class Convergence_Norm_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator            = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;

  public:
    Convergence_Norm_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            Channel_Vector channels,T& result,const unsigned mask)
    {Run(allocator,blocks,channels,result,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Channel_Vector channels,T& result,const unsigned mask) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        double temp_result=(T)0.;
#pragma omp parallel for reduction(+:temp_result)
        for(int b=0;b<blocks.second;b++){
            uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) for(int v=0;v<d;++v) temp_result+=Nova_Utilities::Sqr(allocator.template Get_Const_Array<Struct_type,T>(channels(v))(offset));}
        result=std::sqrt(temp_result);
    }
};
}
#endif

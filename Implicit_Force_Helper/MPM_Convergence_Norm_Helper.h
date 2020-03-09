//!#####################################################################
//! \file MPM_Convergence_Norm_Helper.h
//!#####################################################################
// Class  MPM_Convergence_Norm_Helper
//######################################################################
#ifndef __MPM_Convergence_Norm_Helper__
#define __MPM_Convergence_Norm_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <algorithm>

namespace Nova{
template<class Struct_type,class T,int d>
class MPM_Convergence_Norm_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator            = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;

  public:
    MPM_Convergence_Norm_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            Channel_Vector channels,T& result,const unsigned mask)
    {Run(allocator,blocks,channels,result,mask);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Vector<T Struct_type::*,2> channels,T& result,const unsigned mask) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto c0=allocator.template Get_Const_Array<Struct_type,T>(channels(0)); auto c1=allocator.template Get_Const_Array<Struct_type,T>(channels(1));
        result=0;
        T max_value=0;

#pragma omp parallel for reduction(max:max_value)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) max_value=std::max(max_value,Nova_Utilities::Sqr(c0(offset))+Nova_Utilities::Sqr(c1(offset)));}
        result=std::sqrt(std::max(result,max_value));
    }
    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Vector<T Struct_type::*,3> channels,T& result,const unsigned mask) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto c0=allocator.template Get_Const_Array<Struct_type,T>(channels(0)); auto c1=allocator.template Get_Const_Array<Struct_type,T>(channels(1)); auto c2=allocator.template Get_Const_Array<Struct_type,T>(channels(2));
        result=0;
        T max_value=0;

#pragma omp parallel for reduction(max:max_value)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) max_value=std::max(max_value,Nova_Utilities::Sqr(c0(offset))+Nova_Utilities::Sqr(c1(offset))+Nova_Utilities::Sqr(c2(offset)));}

        result=std::sqrt(std::max(result,max_value));
    }
};
}
#endif

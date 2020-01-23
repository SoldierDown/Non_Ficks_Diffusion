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
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator            = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
  public:
    Inner_Product_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector channels1,Channel_Vector channels2,
                          double& result,const unsigned mask)
    {Run(allocator,blocks,channels1,channels2,result,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector channels1,Channel_Vector channels2,
                          double& result,const unsigned mask) const
    {
        result=(double)0.;
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        double temp_result=0;
#pragma omp parallel for reduction(+:temp_result)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask)for(int v=0;v<d;++v) temp_result+=mass(offset)*allocator.template Get_Const_Array<Struct_type,T>(channels1(v))(offset)*allocator.template Get_Const_Array<Struct_type,T>(channels2(v))(offset);}

        result+=temp_result;
    }
};
}
#endif

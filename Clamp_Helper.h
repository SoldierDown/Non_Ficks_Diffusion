//!#####################################################################
//! \file Clamp_Heler.h
//!#####################################################################
// Class Clamp_Heler
//######################################################################
#ifndef __Clamp_Heler__
#define __Clamp_Heler__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Clamp_Heler
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Clamp_Heler(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel)
    {Run(allocator,blocks,channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel) const
    {
        auto data=allocator.template Get_Array<Struct_type,T>(channel); auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto clamp_heler=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&(Cell_Type_Interior|Cell_Type_Dirichlet)) data(offset)=Nova_Utilities::Clamp(data(offset),(T)0.,(T)1.);
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,clamp_heler);
    }
};
}
#endif

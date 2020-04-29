//!#####################################################################
//! \file Flip_Helper.h
//!#####################################################################
// Class Flip_Helper
//######################################################################
#ifndef __Flip_Helper__
#define __Flip_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Flip_Helper
{
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Flip_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel)
    {Run(allocator,blocks,channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto data=allocator.template Get_Array<Struct_type,T>(channel);

        auto flip_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) 
                if(flags(offset)&Cell_Type_Interior) data(offset)=-data(offset);
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,flip_helper);
    }
};
}
#endif

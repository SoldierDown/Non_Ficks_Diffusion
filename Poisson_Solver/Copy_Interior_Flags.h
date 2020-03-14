//!#####################################################################
//! \file Copy_Interior_Flags.h
//!#####################################################################
// Class Copy_Interior_Flags
//######################################################################
#ifndef __Copy_Interior_Flags__
#define __Copy_Interior_Flags__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Copy_Interior_Flags
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Copy_Interior_Flags(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Hierarchy& hierarchy,const unsigned level)
    {Run(allocator,blocks,hierarchy,level);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Hierarchy& hierarchy,const unsigned level) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto current_flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto copy_interior_flags=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) current_flags(offset)|=Cell_Type_Interior;
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,copy_interior_flags);
    }
};
}
#endif

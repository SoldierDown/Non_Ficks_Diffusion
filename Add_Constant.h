//!#####################################################################
//! \file Add_Constant.h
//!#####################################################################
// Class Add_Constant
//######################################################################
#ifndef __Add_Constant__
#define __Add_Constant__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Add_Constant
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Add_Constant(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* cell_channel,const T constant_value)
    {Run(allocator,blocks,cell_channel,constant_value);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* cell_channel,const T constant_value) const
    {
        auto cell_data=allocator.template Get_Array<Struct_type,T>(cell_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto add_constant=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior)  cell_data(offset)+=constant_value;
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,add_constant);
    }

};
}
#endif

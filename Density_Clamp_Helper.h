//!#####################################################################
//! \file Density_Clamp_Helper.h
//!#####################################################################
// Class Density_Clamp_Helper
//######################################################################
#ifndef __Density_Clamp_Helper__
#define __Density_Clamp_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Density_Clamp_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Density_Clamp_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel)
    {Run(allocator,blocks,density_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel) const
    {
        auto density=allocator.template Get_Array<Struct_type,T>(density_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto density_clamp_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior)  if(density(offset)<(T)0.) density(offset)=(T)0.; // clamp
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,density_clamp_helper);
    }

};
}
#endif

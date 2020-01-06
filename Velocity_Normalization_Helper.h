//!#####################################################################
//! \file Velocity_Normalization_Helper.h
//!#####################################################################
// Class Velocity_Normalization_Helper
//######################################################################
#ifndef __Velocity_Normalization_Helper__
#define __Velocity_Normalization_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Velocity_Normalization_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Velocity_Normalization_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                  Channel_Vector& velocity_channels,T Struct_type::* mass_channel)
    {Run(allocator,blocks,velocity_channels,mass_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Channel_Vector& velocity_channels,T Struct_type::* mass_channel) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(mass_channel);

        auto velocity_normalization_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(mass(offset)>(T)1e-10) for(int v=0;v<d;++v)
                    allocator.template Get_Array<Struct_type,T>(velocity_channels(v))(offset)/=mass(offset);}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,velocity_normalization_helper);
    }
};
}
#endif

//!#####################################################################
//! \file Saturation_Normalization_Helper.h
//!#####################################################################
// Class Saturation_Normalization_Helper
//######################################################################
#ifndef __Saturation_Normalization_Helper__
#define __Saturation_Normalization_Helper__
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Saturation_Normalization_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Saturation_Normalization_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                  T Struct_type::*& saturation_channel,T Struct_type::* void_mass_fluid_channel ,unsigned Struct_type::* flags_channel)
    {Run(allocator,blocks,saturation_channel,void_mass_fluid_channel,flags_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             T Struct_type::*& saturation_channel,T Struct_type::* void_mass_fluid_channel ,unsigned Struct_type::* flags_channel) const
    {
        auto saturation=allocator.template Get_Array<Struct_type,T>(saturation_channel);
        auto void_mass_fluid=allocator.template Get_Const_Array<Struct_type,T>(void_mass_fluid_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(flags_channel);
        auto saturation_normalization_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Node_Saturated) {
                    if(void_mass_fluid(offset)>(T)1.e-7) saturation(offset)/=void_mass_fluid(offset);
                }
                else if(flags(offset)&Node_Active) saturation(offset)=(T)1.;
            }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,saturation_normalization_helper);
        auto clamp_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Node_Saturated) saturation(offset)=Nova_Utilities::Clamp(saturation(offset),(T)0.,(T)1.);
            
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,clamp_helper);
    }
};
}
#endif

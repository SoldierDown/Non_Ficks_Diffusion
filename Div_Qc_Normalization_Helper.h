//!#####################################################################
//! \file Div_Qc_Normalization_Helper.h
//!#####################################################################
// Class Div_Qc_Normalization_Helper
//######################################################################
#ifndef __Div_Qc_Normalization_Helper__
#define __Div_Qc_Normalization_Helper__
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Div_Qc_Normalization_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Div_Qc_Normalization_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                  T Struct_type::* div_Qc_channel,T Struct_type::* volume_channel)
    {Run(allocator,blocks,div_Qc_channel,volume_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             T Struct_type::* div_Qc_channel,T Struct_type::* volume_channel) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto div_Qc=allocator.template Get_Array<Struct_type,T>(div_Qc_channel); auto volume=allocator.template Get_Const_Array<Struct_type,T>(volume_channel);
        auto div_qc_normalization_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){
                    if(volume(offset)>(T)0.) div_Qc(offset)/=volume(offset);
                    else div_Qc(offset)=(T)0.;}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,div_qc_normalization_helper);
    }
};
}
#endif

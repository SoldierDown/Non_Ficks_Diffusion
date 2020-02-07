//!#####################################################################
//! \file Flag_Helper.h
//!#####################################################################
// Class Flag_Helper
//######################################################################
#ifndef __Flag_Helper__
#define __Flag_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>


namespace Nova{
template<class Struct_type,class T,int d>
class Flag_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Flag_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks)
    {Run(allocator,blocks);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks) const
    {
        int n_active=0;
        int n_saturated=0;
        int n_both=0;
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto active_counter=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Node_Active) n_active++;//Log::cout<<c(offset)<<std::endl;
        };
        auto saturated_counter=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Saturated) n_saturated++;//Log::cout<<c(offset)<<std::endl;
        };

        auto sna_counter=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if((flags(offset)&Cell_Saturated)&&(flags(offset)&Node_Active)) n_both++;//Log::cout<<c(offset)<<std::endl;
        };

        
        
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            active_counter(offset);
            saturated_counter(offset);
            sna_counter(offset);
            }
        Log::cout<<"Active: "<<n_active<<", Saturated: "<<n_saturated<<", Both: "<<n_both<<std::endl;
    }
};
}
#endif

//!#####################################################################
//! \file Range_Checker.h
//!#####################################################################
// Class Range_Checker
//######################################################################
#ifndef __Range_Checker__
#define __Range_Checker__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Range_Checker
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Range_Checker(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel,const T min_value,const T max_value)
    {Run(allocator,blocks,channel,min_value,max_value);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::*channel,const T min_value,const T max_value) const
    {
        auto data=allocator.template Get_Const_Array<Struct_type,T>(channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto range_checker=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&(Cell_Type_Interior)) { const T value=data(offset);
                    if(value<min_value||value>max_value) Log::cout<<"Fatal error: "<<value<<std::endl;}}
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            range_checker(offset);}
        
    }
};
}
#endif

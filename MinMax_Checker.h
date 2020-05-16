//!#####################################################################
//! \file MinMax_Checker.h
//!#####################################################################
// Class MinMax_Checker
//######################################################################
#ifndef __MinMax_Checker__
#define __MinMax_Checker__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
namespace Nova{
template<class Struct_type,class T,int d>
class MinMax_Checker
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    MinMax_Checker(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel)
    {Run(allocator,blocks,channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::*channel) const
    {
        T min_value=(T)1e10; T max_value=(T)-1e10;
        auto data=allocator.template Get_Const_Array<Struct_type,T>(channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto minmax_checker=[&](uint64_t offset,T& min_value,T& max_value)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&(Cell_Type_Interior)) { const T value=data(offset);
                    if(value<min_value) min_value=value; if(value>max_value) max_value=value;}}
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            minmax_checker(offset,min_value,max_value);}
        Log::cout<<"min: "<<min_value<<", max: "<<max_value<<std::endl;
    }
};
}
#endif

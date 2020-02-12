//!#####################################################################
//! \file Channel_Vector_Traverse_Helper.h
//!#####################################################################
// Class Channel_Vector_Traverse_Helper
//######################################################################
#ifndef __Channel_Vector_Traverse_Helper__
#define __Channel_Vector_Traverse_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Channel_Vector_Traverse_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Channel_Vector_Traverse_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Channel_Vector& channel_vector,const unsigned mask)
    {Run(allocator,blocks,channel_vector,mask);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Vector<T Struct_type::*,2>& channel_vector,const unsigned mask) const
    {
        using TV                = Vector<T,2>;
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto c0=allocator.template Get_Const_Array<Struct_type,T>(channel_vector(0)); auto c1=allocator.template Get_Const_Array<Struct_type,T>(channel_vector(1));
        T max_norm=-1;
        auto channel_vector_traverse_helper=[&](uint64_t offset,T& max_norm)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) { TV vec({c0(offset),c1(offset)});
                    if(vec.Norm()>max_norm) max_norm=vec.Norm();}
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            channel_vector_traverse_helper(offset,max_norm);}
        Log::cout<<"max norm: "<<max_norm<<std::endl;
    }
};
}
#endif

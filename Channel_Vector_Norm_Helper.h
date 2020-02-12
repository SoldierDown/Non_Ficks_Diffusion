//!#####################################################################
//! \file Channel_Vector_Norm_Helper.h
//!#####################################################################
// Class Channel_Vector_Norm_Helper
//######################################################################
#ifndef __Channel_Vector_Norm_Helper__
#define __Channel_Vector_Norm_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Channel_Vector_Norm_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Channel_Vector_Norm_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector channel_vector)
    {Run(allocator,blocks,channel_vector);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector channel_vector) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        T result=(T)0.;
        auto channel_vector_norm_helper=[&](uint64_t offset,T& result)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) for(int v=0;v<d;++v) result+=Nova_Utilities::Sqr(allocator.template Get_Const_Array<Struct_type,T>(channel_vector(v))(offset));
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            channel_vector_norm_helper(offset,result);}
        Log::cout<<"2-norm: "<<result<<std::endl;
    }
};
}
#endif

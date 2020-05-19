//!#####################################################################
//! \file Flag_Setup_Helper.h
//!#####################################################################
// Class Flag_Setup_Helper
//######################################################################
#ifndef __Flag_Setup_Helper__
#define __Flag_Setup_Helper__
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Flag_Setup_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Flag_Setup_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel)
    {Run(allocator,blocks,channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0); auto flags=allocator.template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto flag_setup_helper=[&](uint64_t offset)
        {
             for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(mass(offset)>(T)0.) {flags(offset)|=Cell_Type_Interior; flags(offset)&=~Cell_Type_Dirichlet;} 
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,flag_setup_helper);
    }
};
}
#endif

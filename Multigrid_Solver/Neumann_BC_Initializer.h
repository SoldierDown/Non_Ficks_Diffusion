//!#####################################################################
//! \file Neumann_BC_Initializer.h
//!#####################################################################
// Class Neumann_BC_Initializer
//######################################################################
#ifndef __Neumann_BC_Initializer__
#define __Neumann_BC_Initializer__

#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Neumann_BC_Initializer
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Neumann_BC_Initializer(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel)
    {Run(allocator,blocks,levelset_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel) const
    {            
        auto flags=allocator.template Get_Array<Struct_type,unsigned>(&Struct_type::flags); auto levelset=allocator.template Get_Const_Array<Struct_type,T>(levelset_channel);
        auto neumann_bc_initializer=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){
                    if(levelset(offset)<(T)0.) {
                        Log::cout<<"signed distance: "<<levelset(offset)<<std::endl;
                        flags(offset)&=~Cell_Type_Interior;}
                }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,neumann_bc_initializer);
    }

};
}
#endif

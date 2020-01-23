//!#####################################################################
//! \file Project_Helper.h
//!#####################################################################
// Class Project_Helper
//######################################################################
#ifndef __Project_Helper__
#define __Project_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Arrays/Array.h>
#include <nova/Tools/Vectors/Vector.h>

#include "../MPM_Plane_Barrier.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Project_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using TV                    = Vector<T,d>;
    using T_Barrier             = MPM_Plane_Barrier<T,d>;
  public:
    Project_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& v_channels,Array<T_Barrier> barriers)
    {Run(allocator,blocks,v_channels,barriers);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& v_channels,Array<T_Barrier> barriers) const
    {
        auto collide_nodes=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch10);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto project_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                for(int id=0;id<barriers.size();++id) if(collide_nodes(offset)==(T)1.) for(int v=0;v<d;++v) allocator.template Get_Array<Struct_type,T>(v_channels(v))(offset)=(T)0.;
                    
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,project_helper);
    }
};
}
#endif





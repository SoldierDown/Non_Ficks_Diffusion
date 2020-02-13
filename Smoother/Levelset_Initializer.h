//!#####################################################################
//! \file Levelset_Initializer.h
//!#####################################################################
// Class Levelset_Initializer
//######################################################################
#ifndef __Levelset_Initializer__
#define __Levelset_Initializer__
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include "Sphere_Levelset.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Levelset_Initializer
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;

  public:
    Levelset_Initializer(const Grid<T,d>& grid,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,Sphere_Levelset<T,d>* sphere_levelset)
    {Run(grid,allocator,blocks,levelset_channel,sphere_levelset);}

    void Run(const Grid<T,d>& grid,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,Sphere_Levelset<T,d>* sphere_levelset) const
    {            
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags); auto levelset=allocator.template Get_Array<Struct_type,T>(levelset_channel);
        auto levelset_initializer=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&(Cell_Type_Interior|Cell_Type_Dirichlet)){
                    Vector<int,d> cell_index(Flag_array_mask::LinearToCoord(offset));
                    levelset(offset)=sphere_levelset->Distance(grid.Center(cell_index));}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,levelset_initializer);
    }

};
}
#endif

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
    using T_INDEX               = Vector<int,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;

  public:
    Levelset_Initializer(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,
                        Sphere_Levelset<T,d> *levelset,const unsigned level)
    {Run(hierarchy,blocks,levelset_channel,levelset,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,
                        Sphere_Levelset<T,d> *levelset,const unsigned level) const
    {            
        const Grid<T,d>& grid=hierarchy.Lattice(level);
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(levelset_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        
        auto levelset_initializer=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            std::array<int,d> base_index_s=Flag_array_mask::LinearToCoord(offset);
            T_INDEX base_index=*reinterpret_cast<T_INDEX*>(&base_index_s);

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior) data(offset)=levelset->Distance(grid.Center(index));
                range_iterator.Next();}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,levelset_initializer);
    }

};
}
#endif

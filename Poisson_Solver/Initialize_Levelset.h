//!#####################################################################
//! \file Initialize_Levelset.h
//!#####################################################################
// Class Initialize_Levelset
//######################################################################
#ifndef __Initialize_Levelset__
#define __Initialize_Levelset__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include "Analytic_Levelset.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Initialize_Levelset
{
    using T_INDEX                   = Vector<int,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Initialize_Levelset(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,
                        Analytic_Levelset<T,d> *levelset,const unsigned level)
    {Run(hierarchy,blocks,levelset_channel,levelset,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,
             Analytic_Levelset<T,d> *levelset,const unsigned level) const
    {
        const Grid<T,d>& grid=hierarchy.Lattice(level);
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(levelset_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto initialize_levelset=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            std::array<int,d> base_index_s=Flag_array_mask::LinearToCoord(offset);
            T_INDEX base_index=*reinterpret_cast<T_INDEX*>(&base_index_s);

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior) data(offset)=levelset->Signed_Distance(grid.Center(index));
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,initialize_levelset);
    }
};
}
#endif

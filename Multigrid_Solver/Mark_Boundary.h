//!#####################################################################
//! \file Mark_Boundary.h
//!#####################################################################
// Class Mark_Boundary
//######################################################################
#ifndef __Mark_Boundary__
#define __Mark_Boundary__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Mark_Boundary
{
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Mark_Boundary(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
                  const Array<uint64_t>& neighbor_offsets,const int level,const unsigned mask)
    {Run(hierarchy,blocks,neighbor_offsets,level,mask);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
             const Array<uint64_t>& neighbor_offsets,const int level,const unsigned mask) const
    {
        auto flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto mark_boundary=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){bool mark_cell=false;
                    for(int k=0;k<(int)neighbor_offsets.size();++k){
                        const uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,neighbor_offsets(k));
                        mark_cell=!hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Interior);
                        if(!mark_cell) mark_cell=hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Dirichlet);
                        if(mark_cell) break;}
                    if(mark_cell) flags(offset)|=mask;}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,mark_boundary);
    }
};
}
#endif

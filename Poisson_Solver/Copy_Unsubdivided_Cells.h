//!#####################################################################
//! \file Copy_Unsubdivided_Cells.h
//!#####################################################################
// Class Copy_Unsubdivided_Cells
//######################################################################
#ifndef __Copy_Unsubdivided_Cells__
#define __Copy_Unsubdivided_Cells__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Copy_Unsubdivided_Cells
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

    enum {number_of_nodes_per_cell  = Topology_Helper::number_of_nodes_per_cell};

  public:
    Copy_Unsubdivided_Cells(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Hierarchy& hierarchy,
                            const uint64_t* nodes_of_cell_offsets,const unsigned level)
    {Run(allocator,blocks,hierarchy,nodes_of_cell_offsets,level);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Hierarchy& hierarchy,
             const uint64_t* nodes_of_cell_offsets,const unsigned level) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto copy_unsubdivided_cells=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&Cell_Type_Interior){
                bool no_children_interior=true;uint64_t base_offset=Flag_array_mask::UpsampleOffset(offset);
                for(int node=0;node<number_of_nodes_per_cell;++node){uint64_t node_offset=Flag_array_mask::Packed_Add(base_offset,nodes_of_cell_offsets[node]);
                    if(hierarchy.template Set<unsigned>(level-1,&Struct_type::flags).Is_Set(node_offset,Cell_Type_Interior)){
                        no_children_interior=false;break;}}
                if(no_children_interior) hierarchy.Activate_Cell(level,offset,Cell_Type_Interior);}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,copy_unsubdivided_cells);
    }
};
}
#endif

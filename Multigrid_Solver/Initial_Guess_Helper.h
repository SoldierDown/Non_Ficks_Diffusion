//!#####################################################################
//! \file Initial_Guess_Helper.h
//!#####################################################################
// Class Initial_Guess_Helper
//######################################################################
#ifndef __Initial_Guess_Helper__
#define __Initial_Guess_Helper__

#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Initial_Guess_Helper
{
    using TV                    = Vector<T,d>;
    using T_INDEX               = Vector<int,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    using Hierarchy_Base        = Grid_Hierarchy<Struct_type,T,d>;
    // using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Initial_Guess_Helper(Hierarchy_Base& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* x_channel,const bool random_value)
    {Run(hierarchy,allocator,blocks,x_channel,random_value);}

    void Run(Hierarchy_Base& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* x_channel,const bool random_value) const
    {            
        Random_Numbers<T> random;
        random.Set_Seed(0);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags); auto x=allocator.template Get_Array<Struct_type,T>(x_channel);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        const Grid<T,d>& grid=hierarchy.Lattice(0);
        TV size_of_domain(grid.domain.max_corner-grid.domain.min_corner);
        auto initial_guess_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){
                    T_INDEX index(Flag_array_mask::LinearToCoord(offset)); TV X=grid.Center(index);
                    TV ratio=(X-grid.domain.min_corner)/size_of_domain;
                    if(random_value) x(offset)=random.Get_Uniform_Number(-1,1);
                    else x(offset)=sin(ratio(0)*two_pi)*sin(ratio(1)*two_pi);}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,initial_guess_helper);
    }

};
}
#endif

//!#####################################################################
//! \file Boundary_Check.h
//!#####################################################################
// Class Boundary_Check
//######################################################################
#ifndef __Boundary_Check__
#define __Boundary_Check__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Boundary_Check
{
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Boundary_Check(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks)
    {   Log::cout<<"Boundary Check"<<std::endl;
        Run(hierarchy,allocator,blocks);}

    void Run(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto boundary_check=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                Vector<int,d> index(Flag_array_mask::LinearToCoord(offset));
                if(flags(offset)&Cell_Type_Dirichlet){
                    Log::cout<<"dirichlet index: "<<index<<std::endl;
                }
                if(flags(offset)&Cell_Type_Interior){
                    Log::cout<<"interior index: "<<index<<std::endl;
                }}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,boundary_check);
    }

};
}
#endif

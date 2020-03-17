//!#####################################################################
//! \file Divergence_Helper.h
//!#####################################################################
// Class Divergence_Helper
//######################################################################
#ifndef __Divergence_Helper__
#define __Divergence_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Divergence_Helper
{
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Divergence_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_velocity_channels,
                      T Struct_type::* divergence_channel,const Vector<uint64_t,d>& other_face_offsets,const int level)
    {Run(hierarchy,blocks,face_velocity_channels,divergence_channel,other_face_offsets,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_velocity_channels,
             T Struct_type::* divergence_channel,const Vector<uint64_t,d>& other_face_offsets,const int level) const
    {
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto divergence=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(divergence_channel);

        double one_over_dx=hierarchy.Lattice(0).one_over_dX(0);
        auto divergence_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset) & (Cell_Type_Interior|Cell_Type_Ghost)){T result=(T)0.;
                for(int axis=0;axis<d;++axis){unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
                    auto face_velocity=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(face_velocity_channels(axis));
                    if(flags(offset)&face_valid_mask) result += one_over_dx*face_velocity(offset);
                    uint64_t other_offset=Flag_array_mask::Packed_Add(offset,other_face_offsets(axis));
                    if(flags(other_offset)&face_valid_mask) result -= one_over_dx*face_velocity(other_offset);}
                divergence(offset)=result;}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,divergence_helper);
    }
};
}
#endif

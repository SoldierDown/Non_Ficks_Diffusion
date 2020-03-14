//!#####################################################################
//! \file Compute_Time_Step.h
//!#####################################################################
// Class Compute_Time_Step
//######################################################################
#ifndef __Compute_Time_Step__
#define __Compute_Time_Step__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Compute_Time_Step
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Compute_Time_Step(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_velocity_channels,
                      Vector<uint64_t,d>& other_face_offsets,const int level,T& dt)
    {Run(hierarchy,blocks,face_velocity_channels,other_face_offsets,level,dt);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_velocity_channels,
             Vector<uint64_t,d>& other_face_offsets,const int level,T& dt) const
    {
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        T max_value=0;

#pragma omp parallel for reduction(max:max_value)
        for(unsigned b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){T local_V_norm=(T)0.;
                for(int axis=0;axis<d;++axis){T axis_side_max=(T)0.;
                    unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
                    uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,other_face_offsets(axis));
                    auto face_velocity=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(face_velocity_channels(axis));

                    if(flags(offset)&face_valid_mask) axis_side_max=std::max(axis_side_max,std::fabs(face_velocity(offset)));
                    if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,face_valid_mask))
                        axis_side_max=std::max(axis_side_max,std::fabs(face_velocity(neighbor_offset)));

                    local_V_norm+=hierarchy.Lattice(level).one_over_dX(axis)*axis_side_max;}
                max_value=std::max(max_value,local_V_norm);}}

        dt=std::max(dt,max_value);
    }
};
}
#endif

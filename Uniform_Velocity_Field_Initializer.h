//!#####################################################################
//! \file Uniform_Velocity_Field_Initializer.h
//!#####################################################################
// Class Uniform_Velocity_Field_Initializer
//######################################################################
#ifndef __Uniform_Velocity_Field_Initializer__
#define __Uniform_Velocity_Field_Initializer__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Uniform_Velocity_Field_Initializer
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Uniform_Velocity_Field_Initializer(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_velocity_channels,const T bv,const int level)
    {
        Run(hierarchy,blocks,face_velocity_channels,bv,level);
    }

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_velocity_channels,const T bv,const int level) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto uniform_velocity_field_initializer=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            TV background_velocity=TV::Axis_Vector(1)*bv;
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                for(int axis=0;axis<d;++axis){const unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis); const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
                    if(!(flags(offset)&face_active_mask)) hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_velocity_channels(axis))(offset)=(T)0.;
                    if(flags(offset)&face_valid_mask)  hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_velocity_channels(axis))(offset)=background_velocity(axis);
                    }
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,uniform_velocity_field_initializer);
    }
};
}
#endif

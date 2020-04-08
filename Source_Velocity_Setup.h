//!#####################################################################
//! \file Source_Velocity_Setup.h
//!#####################################################################
// Class Source_Velocity_Setup
//######################################################################
#ifndef __Source_Velocity_Setup__
#define __Source_Velocity_Setup__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Source_Velocity_Setup
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Implicit_Object_Array     = Array<Implicit_Object<T,d>*>;

  public:
    Source_Velocity_Setup(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Implicit_Object_Array& sources,Array<TV> source_velocity,
                               Channel_Vector& face_velocity_channels,const int level)
    {assert(sources.size()==source_velocity.size());Run(hierarchy,blocks,sources,source_velocity,face_velocity_channels,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Implicit_Object_Array& sources,Array<TV> source_velocity,
             Channel_Vector& face_velocity_channels,const int level) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto source_velocity_setup=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){const T_INDEX index=base_index+range_iterator.Index();
                for(int axis=0;axis<d;++axis){TV X=hierarchy.Lattice(level).Face(axis,index);
                    const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
                    const unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);

                    if(flags(offset)&face_valid_mask){bool inside=false;
                        for(size_t i=0;i<sources.size();++i) if(sources(i)->Inside(X)){
                            hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_velocity_channels(axis))(offset)=source_velocity(i)(axis);inside=true;break;}
                        if(!inside && !(flags(offset)&face_active_mask))  hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_velocity_channels(axis))(offset)=(T)0.;}}
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,source_velocity_setup);
    }
};
}
#endif

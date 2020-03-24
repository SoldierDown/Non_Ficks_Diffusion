//!#####################################################################
//! \file Boundary_Value_Initializer.h
//!#####################################################################
// Class Boundary_Value_Initializer
//######################################################################
#ifndef __Boundary_Value_Initializer__
#define __Boundary_Value_Initializer__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Boundary_Value_Initializer
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

    Array<TV> source_velocity;

  public:
    Boundary_Value_Initializer(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Implicit_Object_Array& sources,
                               Channel_Vector& face_velocity_channels,T Struct_type::* pressure_channel,const int level)
    {
        source_velocity.Append(TV::Axis_Vector(1));
        Run(hierarchy,blocks,sources,face_velocity_channels,pressure_channel,level);
    }

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Implicit_Object_Array& sources,
             Channel_Vector& face_velocity_channels,T Struct_type::* pressure_channel,const int level) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto pressure=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(pressure_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto boundary_value_initializer=[&](uint64_t offset)
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
                        if(!inside && !(flags(offset)&face_active_mask))  hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_velocity_channels(axis))(offset)=(T)0.;
                            
                            }}
                if(flags(offset)&Cell_Type_Dirichlet) pressure(offset)=(T)0.;
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,boundary_value_initializer);
    }
};
}
#endif

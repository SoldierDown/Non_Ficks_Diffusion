//!#####################################################################
//! \file Uniform_Grid_Face_Vector_Advection_Helper.h
//!#####################################################################
// Class Uniform_Grid_Face_Vector_Advection_Helper
//######################################################################
#ifndef __Uniform_Grid_Face_Vector_Advection_Helper__
#define __Uniform_Grid_Face_Vector_Advection_Helper__

#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include "Uniform_Grid_Backtrace.h"
#include "Uniform_Grid_Interpolation_Helper.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Uniform_Grid_Face_Vector_Advection_Helper
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Hierarchy_Interpolation   = Grid_Hierarchy_Interpolation<Struct_type,T,d>;
    enum {number_of_nodes_per_face  = Topology_Helper::number_of_nodes_per_face};

  public:
    Uniform_Grid_Face_Vector_Advection_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
                                                Channel_Vector& face_vector_channels,Channel_Vector& face_velocity_channels,T Struct_type::* result_channel,
                                                const int level,const T dt,const unsigned mask,const int axis)
    {Run(hierarchy,blocks,face_vector_channels,face_velocity_channels,result_channel,level,dt,mask,axis);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
            Channel_Vector& face_vector_channels,Channel_Vector& interpolated_face_velocity_channels,T Struct_type::* result_channel,
            const int level,const T dt,const unsigned mask,const int axis) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto result=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(result_channel);

        auto face_vector_advection_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&mask){TV velocity;
                    for(int axis=0;axis<d;++axis) velocity(axis)=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(interpolated_face_velocity_channels(axis))(offset);
                    // backtrace
                    TV dX=-velocity*dt;
                    uint64_t new_offset=offset;int new_level=level;
                    TV weights=Uniform_Grid_Backtrace_Helper<Struct_type,T,d>::Uniform_Grid_Face_Backtrace(hierarchy,new_level,index,new_offset,dX,axis);
                    T value=Uniform_Grid_Interpolation_Helper<Struct_type,T,d>::Uniform_Grid_Face_Interpolation_Helper(hierarchy,axis,new_level,new_offset,weights,face_vector_channels);
                    result(offset)=value;}
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,face_vector_advection_helper);
    }
};
}
#endif

//!#####################################################################
//! \file Grid_Hierarchy_Projection.h
//!#####################################################################
// Class Grid_Hierarchy_Projection
//######################################################################
#ifndef __Grid_Hierarchy_Projection__
#define __Grid_Hierarchy_Projection__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Projection/Ghost_Value_Accumulate.h>
#include <nova/Dynamics/Hierarchy/Projection/Ghost_Vector_Value_Accumulate.h>
#include <nova/Dynamics/Hierarchy/Projection/Ghost_Value_Propagate.h>
#include <nova/Dynamics/Hierarchy/Projection/Ghost_Vector_Value_Propagate.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include "Clear_Non_Active.h"
#include "Divergence_Helper.h"
#include "Divergence_Star_Helper.h"
#include "Gradient_Helper.h"
#include "Interior_Laplace_Helper.h"
#include "Laplace_Gradient_Helper.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Grid_Hierarchy_Projection
{
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Grid_Hierarchy_Projection() {}
    ~Grid_Hierarchy_Projection() {}

    static void Propagate_Ghost_Values(Hierarchy& hierarchy,T Struct_type::* v_channel)
    {
        for(int level=hierarchy.Levels()-2;level>=0;--level)
            Ghost_Value_Propagate<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),v_channel,v_channel,level);
    }

    static void Propagate_Ghost_Vector_Values(Hierarchy& hierarchy,Channel_Vector& v_channel)
    {
        for(int level=hierarchy.Levels()-2;level>=0;--level)
            Ghost_Vector_Value_Propagate<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),v_channel,v_channel,level);
    }

    static void Accumulate_Ghost_Values(Hierarchy& hierarchy,T Struct_type::* u_channel)
    {
        for(int level=0;level<hierarchy.Levels()-1;++level)
            Ghost_Value_Accumulate<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level+1),u_channel,u_channel,level+1);
    }

    static void Accumulate_Ghost_Vector_Values(Hierarchy& hierarchy,Channel_Vector& u_channel)
    {
        for(int level=0;level<hierarchy.Levels()-1;++level)
            Ghost_Vector_Value_Accumulate<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level+1),u_channel,u_channel,level+1);
    }

    static void Compute_Laplacian(Hierarchy& hierarchy,T Struct_type::* u_channel,T Struct_type::* Lu_channel)
    {
        const int levels=hierarchy.Levels();
        // compute laplace
        for(int level=0;level<levels;++level)
            Interior_Laplace_Helper<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),u_channel,Lu_channel,level);
    }

    static void Compute_Divergence(Hierarchy& hierarchy,Channel_Vector& face_velocity_channels,T Struct_type::* divergence_channel)
    {
        Vector<uint64_t,d> other_face_offsets;
        for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);

        for(int level=0;level<hierarchy.Levels();++level)
            Divergence_Helper<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),face_velocity_channels,
                                               divergence_channel,other_face_offsets,level);
        
    }

    static void Compute_Divergence_Star(Hierarchy& hierarchy,Channel_Vector& face_velocity_channels,T Struct_type::* divergence_channel,Vector<T,d> zeta)
    {
        Vector<uint64_t,d> other_face_offsets;
        for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);

        for(int level=0;level<hierarchy.Levels();++level)
            Divergence_Star_Helper<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),face_velocity_channels,
                                               divergence_channel,zeta,other_face_offsets,level);
        
    }


    static void Compute_Gradient(Hierarchy& hierarchy,Channel_Vector& gradient_channels,T Struct_type::* pressure_channel)
    {
        const int levels=hierarchy.Levels();

        Propagate_Ghost_Values(hierarchy,pressure_channel);

        // clear gradient channels
        for(int axis=0;axis<d;++axis) for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),gradient_channels(axis));

        // compute gradients
        for(int level=0;level<levels;++level)
            Gradient_Helper<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),gradient_channels,pressure_channel,level);

        // accumulate gradients
        for(int axis=0;axis<d;++axis) Accumulate_Ghost_Values(hierarchy,gradient_channels(axis));

        // propagate gradients
        for(int axis=0;axis<d;++axis) Propagate_Ghost_Values(hierarchy,gradient_channels(axis));
    }
};
}
#endif

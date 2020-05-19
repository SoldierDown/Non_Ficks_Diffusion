//!#####################################################################
//! \file Uniform_Grid_Advection_Helper.h
//!#####################################################################
// Class Uniform_Grid_Advection_Helper
//######################################################################
#ifndef __Uniform_Grid_Advection_Helper__
#define __Uniform_Grid_Advection_Helper__

#include "Uniform_Grid_Density_Advection_Helper.h"
#include "Uniform_Grid_Face_Velocity_Advection_Helper.h"
#include "Uniform_Grid_Face_Vector_Advection_Helper.h"
#include "Uniform_Grid_Averaging_Helper.h"
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/SPGrid/Tools/SPGrid_Copy.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Uniform_Grid_Advection_Helper
{
    using TV                        = Vector<T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

    enum {number_of_nodes_per_cell  = Topology_Helper::number_of_nodes_per_cell};
    enum {number_of_nodes_per_face  = Topology_Helper::number_of_nodes_per_face};

  public:
    static void Uniform_Grid_Advect_Density(Hierarchy& hierarchy,Channel_Vector& cell_velocity_channels,T Struct_type::* density_channel,
                                            T Struct_type::* temp_channel,const T dt)
    {
        uint64_t nodes_of_cell_offsets[number_of_nodes_per_cell];
        Topology_Helper::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
        const int levels=hierarchy.Levels();

        // clear temporary channel
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);

        // advect
        for(int level=0;level<levels;++level)
            Uniform_Grid_Density_Advection_Helper<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),cell_velocity_channels,density_channel,temp_channel,
                                                      nodes_of_cell_offsets,level,dt,(unsigned)Cell_Type_Interior);

        // copy result
        for(int level=0;level<levels;++level)
            SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel,
                                                 density_channel,(unsigned)Cell_Type_Interior);
    }

    static void Uniform_Grid_Advect_Face_Velocities(Hierarchy& hierarchy,Channel_Vector& face_velocity_channels,Channel_Vector& face_velocity_backup_channels,
                                                    Channel_Vector& interpolated_face_velocity_channels,T Struct_type::* temp_channel,const T dt)
    {
        Vector<uint64_t,d> other_face_offsets; for(int v=0;v<d;++v) other_face_offsets(v)=Topology_Helper::Axis_Vector_Offset(v);
        Vector<uint64_t,d> negative_face_offsets; for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
        uint64_t nodes_of_cell_offsets[number_of_nodes_per_cell]; Topology_Helper::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
        const int levels=hierarchy.Levels();

        // advect velocities for each axis
        for(int axis=0;axis<d;++axis){unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
            uint64_t nodes_of_face_offsets[number_of_nodes_per_face];
            Topology_Helper::Nodes_Of_Face_Offsets(nodes_of_face_offsets,axis);

            // clear
            for(int level=0;level<levels;++level){for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),interpolated_face_velocity_channels(v));
				SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);}

            // interpolate velocity at face
            for(int level=0;level<levels;++level)
                Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Faces(hierarchy,hierarchy.Allocator(level),hierarchy.Blocks(level),
                                            face_velocity_backup_channels,interpolated_face_velocity_channels,negative_face_offsets,nodes_of_cell_offsets,axis);

            // advect
            for(int level=0;level<levels;++level)
                Uniform_Grid_Face_Velocity_Advection_Helper<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),face_velocity_backup_channels,interpolated_face_velocity_channels,temp_channel,
                                                                other_face_offsets,level,dt,face_active_mask,axis);

            // copy result
            for(int level=0;level<levels;++level)
                SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel,face_velocity_channels(axis),face_active_mask);}
    }

    static void Uniform_Grid_Advect_Face_Vector(Hierarchy& hierarchy,Channel_Vector& face_vector_channels,Channel_Vector& face_velocity_channels,
                                    				Channel_Vector& interpolated_face_velocity_channels,T Struct_type::* temp_channel,const T dt)
    {
        const int levels=hierarchy.Levels();
        Vector<uint64_t,d> other_face_offsets; for(int v=0;v<d;++v) other_face_offsets(v)=Topology_Helper::Axis_Vector_Offset(v);
        Vector<uint64_t,d> negative_face_offsets; for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
        uint64_t nodes_of_cell_offsets[number_of_nodes_per_cell]; Topology_Helper::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
        // advect vector along each axis
        for(int axis=0;axis<d;++axis){unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);

            // clear
            for(int level=0;level<levels;++level){ for(int v=0;v<d;++v) 
				SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),interpolated_face_velocity_channels(v));
				SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);}

            // interpolate velocity at face
            for(int level=0;level<levels;++level)
                Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Faces(hierarchy,hierarchy.Allocator(level),hierarchy.Blocks(level),
                                            face_velocity_channels,interpolated_face_velocity_channels,negative_face_offsets,nodes_of_cell_offsets,axis);

            // advect
            for(int level=0;level<levels;++level)
                Uniform_Grid_Face_Vector_Advection_Helper<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),face_vector_channels,interpolated_face_velocity_channels,temp_channel,
                                                                level,dt,face_active_mask,axis);

            // copy result
            for(int level=0;level<levels;++level)
                SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel,face_vector_channels(axis),face_active_mask);}
    }

};
}
#endif

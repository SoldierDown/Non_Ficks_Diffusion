//!#####################################################################
//! \file Uniform_Grid_Interpolation_Helper.h
//!#####################################################################
// Class Uniform_Grid_Interpolation_Helper
//######################################################################
#ifndef __Uniform_Grid_Interpolation_Helper__
#define __Uniform_Grid_Interpolation_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include "Uniform_Grid_Linear_Interpolation.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Uniform_Grid_Interpolation_Helper
{
    using TV                        = Vector<T,d>;
    using TV2                       = Vector<T,d-1>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Hierarchy_Lookup          = Grid_Hierarchy_Lookup<Struct_type,T,d>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

    enum {number_of_nodes_per_cell  = Topology_Helper::number_of_nodes_per_cell};
    enum {number_of_nodes_per_face  = Topology_Helper::number_of_nodes_per_face};

  public:
    static T Uniform_Grid_Cell_Interpolation_Helper(Hierarchy& hierarchy,uint64_t* nodes_of_cell_offsets,
													const int level,const uint64_t offset,const TV weights,T Struct_type::* cell_channel)
    {
		T density_cell_array[number_of_nodes_per_cell];
        auto cell_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(cell_channel);

        for(int cell=0;cell<number_of_nodes_per_cell;++cell){
            uint64_t cell_offset=Flag_array_mask::Packed_Add(offset,nodes_of_cell_offsets[cell]);
            T cell_value=cell_data(cell_offset);
            density_cell_array[cell]=cell_value;}

        // weights are reversed to account for "reverse" ordering of nodes_of_cell_offsets array
        T multilinear_from_cells=Uniform_Grid_Linear_Interpolation<T,T>::Linear(density_cell_array,weights.Reversed());
        return multilinear_from_cells;
    }
        
    // static T Face_Interpolation_Helper(Hierarchy& hierarchy,uint64_t* nodes_of_face_offsets,const int axis,const int level,const uint64_t offset,
    //                                    const TV2& weights,Channel_Vector& face_velocity_channels,Channel_Vector& node_velocity_channels)
    // {
    //     unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
    //     int face_level=level;uint64_t face_offset=offset;TV2 face_weights=weights;
    //     const bool face_found=Hierarchy_Lookup::Face_Lookup(hierarchy,face_offset,face_level,face_weights,face_valid_mask,axis);

    //     T velocity_node_array[number_of_nodes_per_face],interpolated_face_value=(T)0.;
    //     for(int node=0;node<number_of_nodes_per_face;++node){uint64_t node_offset=Flag_array_mask::Packed_Add(face_offset,nodes_of_face_offsets[node]);
    //         const T node_value=hierarchy.Allocator(face_level).template Get_Const_Array<Struct_type,T>(node_velocity_channels(axis))(node_offset);
    //         interpolated_face_value+=node_value;
    //         velocity_node_array[node]=node_value;}
    //     interpolated_face_value/=(T)number_of_nodes_per_face;

    //     // weights are reversed to account for "reverse" ordering of nodes_of_face_offsets array
    //     const T multilinear_from_nodes=Nova::Linear_Interpolation<T,T>::Linear(velocity_node_array,face_weights.Reversed());
    //     const T phi=(T)2.*std::min(face_weights.Min(),(TV2(1)-face_weights).Min());
    //     T actual_face_value=hierarchy.Allocator(face_level).template Get_Const_Array<Struct_type,T>(face_velocity_channels(axis))(face_offset);

    //     return multilinear_from_nodes + phi*(actual_face_value-interpolated_face_value);
    // }
};
}
#endif

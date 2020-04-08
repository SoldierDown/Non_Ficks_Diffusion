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
	using T_INDEX 					= Vector<int,d>;
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
            density_cell_array[cell]=cell_data(cell_offset);}

        // weights are reversed to account for "reverse" ordering of nodes_of_cell_offsets array
        T multilinear_from_cells=Uniform_Grid_Linear_Interpolation<T,T>::Linear(density_cell_array,weights.Reversed());
        return multilinear_from_cells;
    }
        
    static T Uniform_Grid_Face_Interpolation_Helper(Hierarchy& hierarchy,const int axis,const int level,const uint64_t offset,const TV& weights,Channel_Vector& face_vector_channels)
    {
		T face_array[number_of_nodes_per_cell];
        auto face_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(face_vector_channels(axis));
		Vector<int,d> base_index(Flag_array_mask::LinearToCoord(offset));
		int counter=0;
		for(Range_Iterator<d> iterator(T_INDEX(),T_INDEX(1));iterator.Valid();iterator.Next()){
			const T_INDEX& index=base_index+iterator.Index();
			uint64_t face_offset=Flag_array_mask::Linear_Offset(index._data);
			face_array[counter++]=hierarchy.Channel(level,face_vector_channels(axis))(index._data);}

        // weights are reversed to account for "reverse" ordering of nodes_of_cell_offsets array
        T multilinear_from_faces=Uniform_Grid_Linear_Interpolation<T,T>::Linear(face_array,weights.Reversed());
        return multilinear_from_faces;
    }
};
}
#endif

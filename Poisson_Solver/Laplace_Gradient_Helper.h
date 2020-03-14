//!#####################################################################
//! \file Laplace_Gradient_Helper.h
//!#####################################################################
// Class Laplace_Gradient_Helper
//######################################################################
#ifndef __Laplace_Gradient_Helper__
#define __Laplace_Gradient_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Laplace_Gradient_Helper
{
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Laplace_Gradient_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
                            Channel_Vector& gradient_channels,T Struct_type::* Lu_channel,const int level)
    {Run(hierarchy,blocks,gradient_channels,Lu_channel,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
             Channel_Vector& gradient_channels,T Struct_type::* Lu_channel,const int level) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto Lu_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(Lu_channel);
        auto flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);

        uint64_t face_neighbor_offsets[number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);

        T face_areas[d];
        for(int axis=0;axis<d;++axis){face_areas[axis]=(T)1.;
            for(int other_axis=0;other_axis<d;++other_axis) if(other_axis!=axis) face_areas[axis]*=hierarchy.Lattice(level).dX[other_axis];}

        auto laplace_gradient_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) for(int face=0;face<number_of_faces_per_cell;++face){int axis=(face/2);
                    uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                    auto gradient=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(gradient_channels(axis));

                    if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Ghost)){
                        Lu_data(offset)+=face_areas[axis]*gradient(neighbor_offset);
                        Lu_data(neighbor_offset)-=face_areas[axis]*gradient(neighbor_offset);}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,laplace_gradient_helper);
    }
};
}
#endif

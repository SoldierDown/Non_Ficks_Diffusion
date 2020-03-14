//!#####################################################################
//! \file Interior_Laplace_Helper.h
//!#####################################################################
// Class Interior_Laplace_Helper
//######################################################################
#ifndef __Interior_Laplace_Helper__
#define __Interior_Laplace_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Interior_Laplace_Helper
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Interior_Laplace_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& gradient_channels,
                            T Struct_type::* u_channel,T Struct_type::* Lu_channel,const int level)
    {Run(hierarchy,blocks,gradient_channels,u_channel,Lu_channel,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& gradient_channels,
             T Struct_type::* u_channel,T Struct_type::* Lu_channel,const int level) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto Lu_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(Lu_channel);
        auto data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(u_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        uint64_t face_neighbor_offsets[number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);

        T face_areas[d];
        for(int axis=0;axis<d;++axis){face_areas[axis]=(T)1.;
            for(int other_axis=0;other_axis<d;++other_axis) if(other_axis!=axis) face_areas[axis]*=hierarchy.Lattice(level).dX[other_axis];}

        double scaling_factor=hierarchy.Lattice(0).one_over_dX.Product();

        auto interior_laplace_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            std::array<int,d> base_index_s=Flag_array_mask::LinearToCoord(offset);
            T_INDEX base_index=*reinterpret_cast<T_INDEX*>(&base_index_s);

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior){T laplacian=(T)0.;const TV X=hierarchy.Lattice(level).Center(index);
                    for(int face=0;face<number_of_faces_per_cell;++face){int axis=(face/2),side=(face%2);
                        uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);

                        if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Ghost)){
                            auto gradient=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(gradient_channels(axis));
                            int coarse_level=level+1;T coarse_face_area=(T)1.;
                            // compute coarse face area
                            for(int other_axis=0;other_axis<d;++other_axis) if(other_axis!=axis) coarse_face_area*=hierarchy.Lattice(coarse_level).dX[other_axis];
                            T distance=(T).5*(hierarchy.Lattice(level).dX(axis)+hierarchy.Lattice(coarse_level).dX(axis)),averaging_factor=face_areas[axis]/coarse_face_area;

                            gradient(neighbor_offset)+=scaling_factor*averaging_factor*(data(offset)-data(neighbor_offset))/distance;}
                        else if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Interior|Cell_Type_Dirichlet)){
                            laplacian+=scaling_factor*face_areas[axis]*hierarchy.Lattice(level).one_over_dX[axis]*(data(offset)-data(neighbor_offset));}}
                    Lu_data(offset)=laplacian;}
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,interior_laplace_helper);
    }
};
}
#endif

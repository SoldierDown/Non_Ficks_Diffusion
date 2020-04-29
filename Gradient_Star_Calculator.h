//!#####################################################################
//! \file Gradient_Star_Calculator.h
//!#####################################################################
// Class Gradient_Star_Calculator
//######################################################################
#ifndef __Gradient_Star_Calculator__
#define __Gradient_Star_Calculator__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Gradient_Star_Calculator
{
    using TV                    = Vector<T,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    // Checked
    Gradient_Star_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* cell_channel,Channel_Vector& gradient_channels,
                        const TV one_over_2dx,const TV zata)
    {Run(allocator,blocks,cell_channel,gradient_channels,one_over_2dx,zata);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* cell_channel,Channel_Vector& gradient_channels,
            const TV one_over_2dx,const TV zeta) const
    {
        uint64_t face_neighbor_offsets[number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        auto cell_data=allocator.template Get_Const_Array<Struct_type,T>(cell_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto gradient_star_calculator=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&Cell_Type_Interior) for(int axis=0;axis<d;++axis){   
                auto result_data=allocator.template Get_Array<Struct_type,T>(gradient_channels(axis)); const T coeff=zeta(axis)*one_over_2dx(axis);
                uint64_t neighbor_offset[2];for(int nb=0;nb<2;++nb) neighbor_offset[nb]=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[2*axis+nb]);
                result_data(offset)=coeff*(cell_data(neighbor_offset[1])-cell_data(neighbor_offset[0]));}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,gradient_star_calculator);
    }

};
}
#endif

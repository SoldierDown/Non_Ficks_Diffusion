//!#####################################################################
//! \file Lap_Star_Calculator.h
//!#####################################################################
// Class Lap_Star_Calculator
//######################################################################
#ifndef __Lap_Star_Calculator__
#define __Lap_Star_Calculator__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
// #include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Lap_Star_Calculator
{
    using TV                    = Vector<T,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Lap_Star_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* channel,T Struct_type::* lap_channel,const TV zeta,const T one_over_dx2)
    {Run(allocator,blocks,channel,lap_channel,zeta,one_over_dx2);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* channel,T Struct_type::* lap_channel,const TV zeta,const T one_over_dx2) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto data=allocator.template Get_Const_Array<Struct_type,T>(channel); auto lap=allocator.template Get_Array<Struct_type,T>(lap_channel);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        auto lap_star_calculator=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){
                    for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        const T coeff=one_over_dx2*Nova_Utilities::Sqr(zeta(face/2));
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]); 
                        if(flags(neighbor_offset)&(Cell_Type_Interior|Cell_Type_Dirichlet)) lap(offset)+=coeff*(data(neighbor_offset)-data(offset));}}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,lap_star_calculator);


    }
};
}
#endif

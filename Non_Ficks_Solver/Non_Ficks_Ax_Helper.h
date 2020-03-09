//!#####################################################################
//! \file Non_Ficks_Ax_Helper.h
//!#####################################################################
// Class Laplace_Helper
//######################################################################
#ifndef __Non_Ficks_Ax_Helper__
#define __Non_Ficks_Ax_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Utilities.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Non_Ficks_Ax_Helper
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Non_Ficks_Ax_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
                   T Struct_type::* u_channel,T Struct_type::* Lu_channel,const T coeff1,const int level)
    {Run(hierarchy,blocks,u_channel,Lu_channel,coeff1,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
             T Struct_type::* u_channel,T Struct_type::* Lu_channel,const T coeff1,const int level) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto Lu_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(Lu_channel);
        auto data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(u_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        uint64_t face_neighbor_offsets[number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);

        T scaling_factor=coeff1*Nova_Utilities::Sqr<T>(hierarchy.Lattice(level).one_over_dX[0]);

        auto non_ficks_ax_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){T Ax_entry=data(offset);
                    for(int face=0;face<number_of_faces_per_cell;++face){
                        uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Interior))
                            Ax_entry+=scaling_factor*(data(offset)-data(neighbor_offset));
                        else if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Dirichlet))
                            Ax_entry+=scaling_factor*data(offset);}
                    Lu_data(offset)=Ax_entry;}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,non_ficks_ax_helper);
    }
};
}
#endif

//!#####################################################################
//! \file Epsilon_Squared_Grad_S_Helper.h
//!#####################################################################
// Class Epsilon_Squared_Grad_S_Helper
//######################################################################
#ifndef __Epsilon_Squared_Grad_S_Helper__
#define __Epsilon_Squared_Grad_S_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Epsilon_Squared_Grad_S_Helper
{
    using TV                        = Vector<T,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Epsilon_Squared_Grad_S_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_qs_channels,
                   T Struct_type::* epsilon_channel,T Struct_type::* density_backup_channel,const T minus_dt_over_tau_1,const TV one_over_dx,const int level)
    {Run(allocator,blocks,face_qs_channels,epsilon_channel,density_backup_channel,minus_dt_over_tau_1,one_over_dx,level);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_qs_channels,
            T Struct_type::* epsilon_channel,T Struct_type::* density_backup_channel,const T minus_dt_over_tau_1,const TV one_over_dx,const int level) const
    {
        auto density_backup_data=allocator.template Get_Const_Array<Struct_type,T>(density_backup_channel);
        auto epsilon_data=allocator.template Get_Const_Array<Struct_type,T>(epsilon_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        Vector<uint64_t,d> negative_face_offsets;
        for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
        auto epsilon_squared_grad_s_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) for(int axis=0;axis<d;++axis){
                auto face_qp_data=allocator.template Get_Array<Struct_type,T>(face_qs_channels(axis));
                if(flags(offset)&Topology_Helper::Face_Active_Mask(axis)){uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,negative_face_offsets(axis));
                    const T interpolated_epsilon_value=(T).5*(epsilon_data(offset)+epsilon_data(neighbor_offset));
                    face_qp_data(offset)+=minus_dt_over_tau_1*Nova_Utilities::Sqr<T>(interpolated_epsilon_value)*(density_backup_data(offset)-density_backup_data(neighbor_offset))*one_over_dx(axis);}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,epsilon_squared_grad_s_helper);
    }
};
}
#endif

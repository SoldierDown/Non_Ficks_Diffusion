//!#####################################################################
//! \file DG_Ficks_T_RHS_Helper.h
//!#####################################################################
// Class DG_Ficks_T_RHS_Helper
//######################################################################
#ifndef __DG_Ficks_T_RHS_Helper__
#define __DG_Ficks_T_RHS_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
// #include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class DG_Ficks_T_RHS_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    DG_Ficks_T_RHS_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                        T Struct_type::* density_channel,T Struct_type::* density_backup_channel,
                        T Struct_type::* T_channel,T Struct_type::* rhs_channel,
                        const T dt_over_dx2,const T SR_dt)
    {Run(allocator,blocks,density_channel,density_backup_channel,T_channel,rhs_channel,dt_over_dx2,SR_dt);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            T Struct_type::* density_channel,T Struct_type::* density_backup_channel,
            T Struct_type::* T_channel,T Struct_type::* rhs_channel,
            const T dt_over_dx2,const T SR_dt) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto rhs_data=allocator.template Get_Array<Struct_type,T>(rhs_channel);  
        auto density_data=allocator.template Get_Const_Array<Struct_type,T>(density_channel); 
        auto density_backup_data=allocator.template Get_Const_Array<Struct_type,T>(density_backup_channel); 
        auto T_data=allocator.template Get_Const_Array<Struct_type,T>(T_channel); 
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        auto dg_ficks_t_rhs_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){ 
                    rhs_data(offset)=T_data(offset)+SR_dt+density_data(offset)-density_backup_data(offset);
                    for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(flags(neighbor_offset)&Cell_Type_Dirichlet) rhs_data(offset)+=dt_over_dx2*T_data(neighbor_offset);}}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,dg_ficks_t_rhs_helper);
    }

};
}
#endif

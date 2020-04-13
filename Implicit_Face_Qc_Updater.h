//!#####################################################################
//! \file Implicit_Face_Qc_Updater.h
//!#####################################################################
// Class Implicit_Face_Qc_Updater
//######################################################################
#ifndef __Implicit_Face_Qc_Updater__
#define __Implicit_Face_Qc_Updater__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Implicit_Face_Qc_Updater
{
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Implicit_Face_Qc_Updater(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_qc_channels,Channel_Vector& face_qc_backup_channels,
                   T Struct_type::* density_channel,const T coeff1,const T coeff2,const int level)
    {Run(hierarchy,blocks,face_qc_channels,face_qc_backup_channels,density_channel,coeff1,coeff2,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_qc_channels,Channel_Vector& face_qc_backup_channels,
                   T Struct_type::* density_channel,const T coeff1,const T coeff2,const int level) const
    {
        auto density=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(density_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        Vector<uint64_t,d> negative_face_offsets;
        for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
        T one_over_dx=hierarchy.Lattice(0).one_over_dX(0);
        auto implicit_face_qc_updater=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) for(int axis=0;axis<d;++axis){
                auto face_qc_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_qc_channels(axis));
                auto face_qc_backup_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_qc_backup_channels(axis));
                if(flags(offset)&Topology_Helper::Face_Active_Mask(axis)){uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,negative_face_offsets(axis));
                    face_qc_data(offset)=coeff1*face_qc_backup_data(offset)+coeff2*(density(offset)-density(neighbor_offset))*one_over_dx;}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,implicit_face_qc_updater);
    }
};
}
#endif

//!#####################################################################
//! \file Face_Vector_Copy.h
//!#####################################################################
// Class Face_Vector_Copy
//######################################################################
#ifndef __Face_Vector_Copy__
#define __Face_Vector_Copy__

#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Face_Vector_Copy
{
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Face_Vector_Copy(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_vector_channels,
                   Channel_Vector& face_vector_backup_channels)
    {Run(allocator,blocks,face_vector_channels,face_vector_backup_channels);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_vector_channels,
                   Channel_Vector& face_vector_backup_channels) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto face_vector_copy=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) for(int axis=0;axis<d;++axis){
                auto face_vector_data=allocator.template Get_Array<Struct_type,T>(face_vector_channels(axis));
                auto face_vector_backup_data=allocator.template Get_Array<Struct_type,T>(face_vector_backup_channels(axis));
                if(flags(offset)&Topology_Helper::Face_Active_Mask(axis)) face_vector_backup_data(offset)=face_vector_data(offset);}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,face_vector_copy);
    }
};
}
#endif

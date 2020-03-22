//!#####################################################################
//! \file Non_Ficks_Smoke_Density_Explicit_Update_Helper.h
//!#####################################################################
// Class Non_Ficks_Smoke_Density_Explicit_Update_Helper
//######################################################################
#ifndef __Non_Ficks_Smoke_Density_Explicit_Update_Helper__
#define __Non_Ficks_Smoke_Density_Explicit_Update_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
// #include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Non_Ficks_Smoke_Density_Explicit_Update_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Non_Ficks_Smoke_Density_Explicit_Update_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                  T Struct_type::* density_channel,T Struct_type::* lap_density_channel,T Struct_type::* div_qc_channel,const T coeff1,const T coeff2)
    {Run(allocator,blocks,density_channel,lap_density_channel,div_qc_channel,coeff1,coeff2);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* density_channel,T Struct_type::* lap_density_channel,T Struct_type::* div_qc_channel,const T coeff1,const T coeff2) const
    {
        auto density=allocator.template Get_Array<Struct_type,T>(density_channel);
        auto lap_density=allocator.template Get_Const_Array<Struct_type,T>(lap_density_channel);
        auto div_qc=allocator.template Get_Array<Struct_type,T>(div_qc_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto non_ficks_smoke_density_update_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) density(offset)+=coeff1*lap_density(offset)+coeff2*div_qc(offset);
                
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,non_ficks_smoke_density_update_helper);
    }

};
}
#endif

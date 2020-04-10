//!#####################################################################
//! \file Epsilon_Squared_Lap_S_Helper.h
//!#####################################################################
// Class Epsilon_Squared_Lap_S_Helper
//######################################################################
#ifndef __Epsilon_Squared_Lap_S_Helper__
#define __Epsilon_Squared_Lap_S_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Epsilon_Squared_Lap_S_Helper
{
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Epsilon_Squared_Lap_S_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel,
                   T Struct_type::* epsilon_channel,T Struct_type::* lap_density_backup_channel,const T dt_over_tau_s,const int level)
    {Run(allocator,blocks,density_channel,epsilon_channel,lap_density_backup_channel,dt_over_tau_s,level);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel,
                   T Struct_type::* epsilon_channel,T Struct_type::* lap_density_backup_channel,const T dt_over_tau_s,const int level) const
    {
        auto density_data=allocator.template Get_Array<Struct_type,T>(density_channel);
        auto epsilon_data=allocator.template Get_Const_Array<Struct_type,T>(epsilon_channel);
        auto lap_density_backup_data=allocator.template Get_Const_Array<Struct_type,T>(lap_density_backup_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto epsilon_squared_lap_s_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior)
                    density_data(offset)+=dt_over_tau_s*Nova_Utilities::Sqr(epsilon_data(offset))*lap_density_backup_data(offset);}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,epsilon_squared_lap_s_helper);
    }
};
}
#endif

//!#####################################################################
//! \file Density_Plus_Poly_Term_Helper.h
//!#####################################################################
// Class Density_Plus_Poly_Term_Helper
//######################################################################
#ifndef __Density_Plus_Poly_Term_Helper__
#define __Density_Plus_Poly_Term_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Tools/Utilities/Constants.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Density_Plus_Poly_Term_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Density_Plus_Poly_Term_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                    T Struct_type::* density_channel,T Struct_type::* density_backup_channel,T Struct_type::* m_channel,const T dt_over_taus)
    {Run(allocator,blocks,density_channel,density_backup_channel,m_channel,dt_over_taus);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,                                    
            T Struct_type::* density_channel,T Struct_type::* density_backup_channel,T Struct_type::* m_channel,const T dt_over_taus) const
    {
        auto density_data=allocator.template Get_Array<Struct_type,T>(density_channel);
        auto density_backup_data=allocator.template Get_Const_Array<Struct_type,T>(density_backup_channel);
        auto m_data=allocator.template Get_Const_Array<Struct_type,T>(m_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto density_plus_poly_term_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){ const T cell_density=density_backup_data(offset);
                    density_data(offset)+=dt_over_taus*cell_density*((T)1.-cell_density)*(cell_density+m_data(offset)-(T).5);}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,density_plus_poly_term_helper);
    }

};
}
#endif

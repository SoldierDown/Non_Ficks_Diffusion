//!#####################################################################
//! \file DG_M_Calculator.h
//!#####################################################################
// Class DG_M_Calculator
//######################################################################
#ifndef __DG_M_Calculator__
#define __DG_M_Calculator__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Tools/Utilities/Constants.h>

namespace Nova{
template<class Struct_type,class T,int d>
class DG_M_Calculator
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    DG_M_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                    T Struct_type::* m_channel,T Struct_type::* sigma_channel,T Struct_type::* T_backup_channel,
                    const T m_alpha,const T gamma,const T Teq)
    {Run(allocator,blocks,m_channel,sigma_channel,T_backup_channel,m_alpha,gamma,Teq);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            T Struct_type::* m_channel,T Struct_type::* sigma_channel,T Struct_type::* T_backup_channel,
            const T m_alpha,const T gamma,const T Teq) const
    {
        auto m_data=allocator.template Get_Array<Struct_type,T>(m_channel); 
        auto sigma_data=allocator.template Get_Const_Array<Struct_type,T>(sigma_channel);
        auto T_backup_data=allocator.template Get_Const_Array<Struct_type,T>(T_backup_channel); 
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto dg_m_calculator=[&](uint64_t offset)
        {
            const T m_alpha_over_pi=m_alpha/pi;
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior)  { const T sig=sigma_data(offset);
                    m_data(offset)=-m_alpha_over_pi*atan(gamma*sig*T_backup_data(offset));}
                    // m_data(offset)=m_alpha_over_pi*atan(gamma*sig*(Teq-T_backup_data(offset)));}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,dg_m_calculator);
    }

};
}
#endif

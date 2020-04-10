//!#####################################################################
//! \file M_Calculator.h
//!#####################################################################
// Class M_Calculator
//######################################################################
#ifndef __M_Calculator__
#define __M_Calculator__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Tools/Utilities/Constants.h>

namespace Nova{
template<class Struct_type,class T,int d>
class M_Calculator
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    M_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* m_channel,T Struct_type::* T_channel,const T m_alpha)
    {Run(allocator,blocks,m_channel,T_channel,m_alpha);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* m_channel,T Struct_type::* T_channel,const T m_alpha) const
    {
        auto m_data=allocator.template Get_Array<Struct_type,T>(m_channel);
        auto T_data=allocator.template Get_Const_Array<Struct_type,T>(T_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto m_calculator=[&](uint64_t offset)
        {
            const T m_alpha_over_pi=m_alpha/pi; const T gamma=(T)10.; const T Teq=(T)1.;
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior)  m_data(offset)=m_alpha_over_pi*atan(gamma*(Teq-T_data(offset)));
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,m_calculator);
    }

};
}
#endif

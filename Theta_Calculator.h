//!#####################################################################
//! \file Theta_Calculator.h
//!#####################################################################
// Class Theta_Calculator
//######################################################################
#ifndef __Theta_Calculator__
#define __Theta_Calculator__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Theta_Calculator
{
    using TV                        = Vector<T,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Theta_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& v_channels,T Struct_type::* theta_channel)
    {Run(allocator,blocks,v_channels,theta_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& v_channels,T Struct_type::* theta_channel) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto theta_data=allocator.template Get_Array<Struct_type,T>(theta_channel);
        auto v0_data=allocator.template Get_Const_Array<Struct_type,T>(v_channels(0)); auto v1_data=allocator.template Get_Const_Array<Struct_type,T>(v_channels(1));
        auto theta_calculator=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){ const T v0=v0_data(offset); const T v1=v1_data(offset); const T half_pi=(T).5*pi;
                    if(v0==(T)0.){ if(v1>(T)0.) theta_data(offset)=half_pi; else if(v1<(T)0.) theta_data(offset)=-half_pi;}
                    else{ const T atan_value=atan(v1/v0);
                        if(v0>(T)0.){ if(v1>(T)0.) theta_data(offset)=atan_value; else if(v1<(T)0.) theta_data(offset)=two_pi+atan_value;}
                        else theta_data(offset)=pi+atan_value;}   
                }
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,theta_calculator);
    }
};

}
#endif

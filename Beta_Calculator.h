//!#####################################################################
//! \file Beta_Calculator.h
//!#####################################################################
// Class Beta_Calculator
//######################################################################
#ifndef __Beta_Calculator__
#define __Beta_Calculator__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Beta_Calculator
{
    using TV                        = Vector<T,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Beta_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& dSdX_channels,T Struct_type::* beta_channel)
    {Run(allocator,blocks,dSdX_channels,beta_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& dSdX_channels,T Struct_type::* beta_channel) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto beta_data=allocator.template Get_Array<Struct_type,T>(beta_channel);
        auto dSdx_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(0)); 
        auto dSdy_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(1));
        auto dSdz_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(2)); 
        auto beta_calculator=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){ 
                    const T minus_dSdz=dSdz_data(offset);
                    const T partial_sqrt=sqrt(Nova_Utilities::Sqr(dSdx_data(offset))+Nova_Utilities::Sqr(dSdy_data(offset))); const T half_pi=(T).5*pi;
                    if(minus_dSdz==(T)0.) beta_data(offset)=half_pi;
                    else{ const T atan_value=atan(partial_sqrt/minus_dSdz);
                        if(minus_dSdz>(T)0.) beta_data(offset)=atan_value; 
                        else beta_data(offset)=pi+atan_value;} 
                }
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,beta_calculator);
    }
};

}
#endif

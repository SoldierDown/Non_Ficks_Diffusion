//!#####################################################################
//! \file Av_Calculator.h
//!#####################################################################
// Class Av_Calculator
//######################################################################
#ifndef __Av_Calculator__
#define __Av_Calculator__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Av_Calculator
{
    using TV                        = Vector<T,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Av_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* Av_channel,T Struct_type::* theta_channel,T Struct_type::* beta_channel,
                const Vector<T,2> epsilon,const T delta,const int omega,const int zeta)
    {Run(allocator,blocks,Av_channel,theta_channel,beta_channel,epsilon,delta,omega,zeta);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            T Struct_type::* Av_channel,T Struct_type::* theta_channel,T Struct_type::* beta_channel,
            const Vector<T,2> epsilon,const T delta,const int omega,const int zeta) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags); 
        auto Av_data=allocator.template Get_Array<Struct_type,T>(Av_channel);
        auto theta_data=allocator.template Get_Const_Array<Struct_type,T>(theta_channel);
        const T eps_xy=epsilon(0);
        auto av_calculator=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) Av_data(offset)=delta*((T)1.+eps_xy*cos(omega*theta_data(offset)));
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,av_calculator);
    }
    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            T Struct_type::* Av_channel,T Struct_type::* theta_channel,T Struct_type::* beta_channel,
            const Vector<T,2> epsilon,const T delta,const int omega,const int zeta) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags); 
        auto Av_data=allocator.template Get_Array<Struct_type,T>(Av_channel);
        auto theta_data=allocator.template Get_Const_Array<Struct_type,T>(theta_channel); 
        auto beta_data=allocator.template Get_Const_Array<Struct_type,T>(beta_channel);
        const T eps_xy=epsilon(0); const T eps_z=epsilon(1);
        auto av_calculator=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) Av_data(offset)=delta*((T)1.+eps_xy*cos(omega*theta_data(offset))+eps_z*cos(zeta*beta_data(offset)));
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,av_calculator);
    }

};

}
#endif

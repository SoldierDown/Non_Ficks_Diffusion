//!#####################################################################
//! \file Epsilon_Epsilon_Prime_dSdX_Helper.h
//!#####################################################################
// Class Epsilon_Epsilon_Prime_dSdX_Helper
//######################################################################
#ifndef __Epsilon_Epsilon_Prime_dSdX_Helper__
#define __Epsilon_Epsilon_Prime_dSdX_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Tools/Utilities/Constants.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Epsilon_Epsilon_Prime_dSdX_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Epsilon_Epsilon_Prime_dSdX_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            Channel_Vector& dSdX_channels,T Struct_type::* psi_channel,T Struct_type::* epsilon_channel,
                            const int omega,const T delta)
    {Run(allocator,blocks,dSdX_channels,psi_channel,epsilon_channel,omega,delta);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Channel_Vector& dSdX_channels,T Struct_type::* psi_channel,T Struct_type::* epsilon_channel,
            const int omega,const T delta) const
    {
        auto dSdx_data=allocator.template Get_Array<Struct_type,T>(dSdX_channels(0));
        auto dSdy_data=allocator.template Get_Array<Struct_type,T>(dSdX_channels(1));
        auto psi_data=allocator.template Get_Const_Array<Struct_type,T>(psi_channel);
        auto epsilon_data=allocator.template Get_Array<Struct_type,T>(epsilon_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto epsilon_epsilon_prime_dSdX_helper=[&](uint64_t offset)
        {
            const T coeff1=(T).01; const T coeff2=(T).01*delta; const T coeff3=(T)-.01*delta*omega;
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&Cell_Type_Interior){
                const T psi=psi_data(offset); const T epsilon=coeff1+coeff2*cos(omega*psi); const T epsilon_prime=coeff3*sin(omega*psi);
                for(int axis=0;axis<d;++axis) allocator.template Get_Array<Struct_type,T>(dSdX_channels(axis))(offset)*=epsilon*epsilon_prime;
                epsilon_data(offset)=epsilon;}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,epsilon_epsilon_prime_dSdX_helper);
    }


};
}
#endif
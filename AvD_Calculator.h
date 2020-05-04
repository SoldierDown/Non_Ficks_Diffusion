//!#####################################################################
//! \file AvD_Calculator.h
//!#####################################################################
// Various classes for calculating AvD.
//######################################################################
#ifndef __AvD_Calculator__
#define __AvD_Calculator__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
namespace Nova{
template<class Struct_type,class T,int d>
class Compute_AvD0
{
    using TV                        = Vector<T,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Compute_AvD0(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                Channel_Vector& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD0_channel,T Struct_type::* theta_channel,
                const Vector<T,2> epsilon,const T delta,const int omega)
    {Run(allocator,blocks,dSdX_channels,Av_channel,AvD0_channel,theta_channel,epsilon,delta,omega);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,2>& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD0_channel,T Struct_type::* theta_channel,
            const Vector<T,2> epsilon,const T delta,const int omega) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto AvD0_data=allocator.template Get_Array<Struct_type,T>(AvD0_channel); auto Av_data=allocator.template Get_Const_Array<Struct_type,T>(Av_channel);
        auto dSdy_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(1));
        auto theta_data=allocator.template Get_Const_Array<Struct_type,T>(theta_channel);
        const T eps_xy=epsilon(0);
        auto compute_avd0=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) {const T Av_prime=-delta*eps_xy*omega*sin(omega*theta_data(offset));
                    AvD0_data(offset)=-Av_data(offset)*Av_prime*dSdy_data(offset);}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_avd0);
    }
    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,3>& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD0_channel,T Struct_type::* theta_channel,
            const Vector<T,2> epsilon,const T delta,const int omega) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto AvD0_data=allocator.template Get_Array<Struct_type,T>(AvD0_channel); auto Av_data=allocator.template Get_Const_Array<Struct_type,T>(Av_channel);
        auto theta_data=allocator.template Get_Const_Array<Struct_type,T>(theta_channel);
        auto dSdx_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(0)); 
        auto dSdy_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(1));
        auto dSdz_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(2));
        const T eps_xy=epsilon(0); const T eps_z=epsilon(1);
        auto compute_avd0=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){ 
                    const T dSdx=dSdx_data(offset); const T dSdy=dSdy_data(offset); const T dSdz=dSdz_data(offset); 
                    const T dSdx2_plus_dSdy2=Nova_Utilities::Sqr(dSdx)+Nova_Utilities::Sqr(dSdy); 
                    const T dSdX_squared=dSdx2_plus_dSdy2+Nova_Utilities::Sqr(dSdz);
                    const T Av=Av_data(offset); const T Av_prime=-delta*eps_xy*omega*sin(omega*theta_data(offset));
                    const T FP=Av_prime*dSdy/dSdx2_plus_dSdy2;
                    const T SP=(T)4.*eps_z*delta*dSdx*Nova_Utilities::Sqr(dSdz)/Nova_Utilities::Sqr(dSdX_squared);
                    AvD0_data(offset)=-dSdX_squared*Av*(FP+SP);}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_avd0);
    }
};

template<class Struct_type,class T,int d>
class Compute_AvD1
{
    using TV                        = Vector<T,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Compute_AvD1(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                Channel_Vector& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD1_channel,T Struct_type::* theta_channel,
                const Vector<T,2> epsilon,const T delta,const int omega)
    {Run(allocator,blocks,dSdX_channels,Av_channel,AvD1_channel,theta_channel,epsilon,delta,omega);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,2>& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD1_channel,T Struct_type::* theta_channel,
            const Vector<T,2> epsilon,const T delta,const int omega) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto AvD1_data=allocator.template Get_Array<Struct_type,T>(AvD1_channel); auto Av_data=allocator.template Get_Const_Array<Struct_type,T>(Av_channel);
        auto theta_data=allocator.template Get_Const_Array<Struct_type,T>(theta_channel);
        auto dSdx_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(0));
        const T eps_xy=epsilon(0);
        auto compute_avd1=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) { const T Av_prime=-delta*eps_xy*omega*sin(omega*theta_data(offset));
                    AvD1_data(offset)=Av_data(offset)*Av_prime*dSdx_data(offset);}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_avd1);
    }
    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,3>& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD1_channel,T Struct_type::* theta_channel,
            const Vector<T,2> epsilon,const T delta,const int omega) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto AvD1_data=allocator.template Get_Array<Struct_type,T>(AvD1_channel); auto Av_data=allocator.template Get_Const_Array<Struct_type,T>(Av_channel);
        auto theta_data=allocator.template Get_Const_Array<Struct_type,T>(theta_channel);
        auto dSdx_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(0)); 
        auto dSdy_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(1));
        auto dSdz_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(2));
        const T eps_xy=epsilon(0); const T eps_z=epsilon(1);
        auto compute_avd1=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){ 
                    const T dSdx=dSdx_data(offset); const T dSdy=dSdy_data(offset); const T dSdz=dSdz_data(offset);
                    const T dSdx2_plus_dSdy2=Nova_Utilities::Sqr(dSdx)+Nova_Utilities::Sqr(dSdy); 
                    const T dSdX_squared=dSdx2_plus_dSdy2+Nova_Utilities::Sqr(dSdz);
                    const T Av=Av_data(offset); const T Av_prime=-delta*eps_xy*omega*sin(omega*theta_data(offset));
                    const T FP=-Av_prime*dSdx/dSdx2_plus_dSdy2;
                    const T SP=(T)4.*eps_z*delta*dSdy*delta*Nova_Utilities::Sqr(dSdz)/Nova_Utilities::Sqr(dSdX_squared);
                    AvD1_data(offset)=-dSdX_squared*Av*(FP+SP);}
                                 
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_avd1);
    }
};


template<class Struct_type,class T,int d>
class Compute_AvD2
{
    using TV                        = Vector<T,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Compute_AvD2(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                Channel_Vector& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD2_channel,
                const Vector<T,2> epsilon,const T delta,const int omega)
    {Run(allocator,blocks,dSdX_channels,Av_channel,AvD2_channel,epsilon,delta,omega);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,2>& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD2_channel,
            const Vector<T,2> epsilon,const T delta,const int omega) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto AvD2_data=allocator.template Get_Array<Struct_type,T>(AvD2_channel);
        auto compute_avd2=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) AvD2_data(offset)=(T)0.;
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_avd2);
    }

    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,3>& dSdX_channels,T Struct_type::* Av_channel,T Struct_type::* AvD2_channel,
            const Vector<T,2> epsilon,const T delta,const int omega) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto AvD2_data=allocator.template Get_Array<Struct_type,T>(AvD2_channel);
        auto Av_data=allocator.template Get_Const_Array<Struct_type,T>(Av_channel);
        auto dSdx_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(0)); 
        auto dSdy_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(1));
        auto dSdz_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(2));
        const T eps_z=epsilon(1);
        auto compute_avd2=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){ 
                    const T dSdx=dSdx_data(offset); const T dSdy=dSdy_data(offset); const T dSdz=dSdz_data(offset);
                    const T dSdx2_plus_dSdy2=Nova_Utilities::Sqr(dSdx)+Nova_Utilities::Sqr(dSdy); 
                    const T dSdX_squared=dSdx2_plus_dSdy2+Nova_Utilities::Sqr(dSdz);
                    const T Av=Av_data(offset);
                    const T D2=-(T)4.*eps_z*delta*dSdz*dSdx2_plus_dSdy2/Nova_Utilities::Sqr(dSdX_squared);
                    AvD2_data(offset)=-dSdX_squared*Av*D2;}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_avd2);
    }
};


}
#endif

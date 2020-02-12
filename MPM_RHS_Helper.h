//!#####################################################################
//! \file MPM_RHS_Helper.h
//!#####################################################################
// Class MPM_RHS_Helper
//######################################################################
#ifndef __MPM_RHS_Helper__
#define __MPM_RHS_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>


namespace Nova{
template<class Struct_type,class T,int d>
class MPM_RHS_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    MPM_RHS_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Channel_Vector velocity_star_channels,Channel_Vector& rhs_channels)
    {Run(allocator,blocks,velocity_star_channels,rhs_channels);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Vector<T Struct_type::*,2>& velocity_star_channels,Vector<T Struct_type::*,2>& rhs_channels) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto vs0=allocator.template Get_Const_Array<Struct_type,T>(velocity_star_channels(0)); auto vs1=allocator.template Get_Const_Array<Struct_type,T>(velocity_star_channels(1));
        auto rhs0=allocator.template Get_Array<Struct_type,T>(rhs_channels(0)); auto rhs1=allocator.template Get_Array<Struct_type,T>(rhs_channels(1));
        auto mpm_rhs_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){rhs0(offset)=vs0(offset); rhs1(offset)=vs1(offset);}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,mpm_rhs_helper);        
    }

    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Vector<T Struct_type::*,3>& velocity_star_channels,Vector<T Struct_type::*,3>& rhs_channels) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto vs0=allocator.template Get_Const_Array<Struct_type,T>(velocity_star_channels(0)); auto vs1=allocator.template Get_Const_Array<Struct_type,T>(velocity_star_channels(1)); auto vs2=allocator.template Get_Const_Array<Struct_type,T>(velocity_star_channels(2));
        auto rhs0=allocator.template Get_Array<Struct_type,T>(rhs_channels(0)); auto rhs1=allocator.template Get_Array<Struct_type,T>(rhs_channels(1)); auto rhs2=allocator.template Get_Array<Struct_type,T>(rhs_channels(2));
        auto mpm_rhs_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){rhs0(offset)=vs0(offset); rhs1(offset)=vs1(offset); rhs2(offset)=vs2(offset);}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,mpm_rhs_helper);        
    }

};
}
#endif

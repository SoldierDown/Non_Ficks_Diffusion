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
#include "MPM_Flags.h"


namespace Nova{
template<class Struct_type,class T,int d>
class MPM_RHS_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    MPM_RHS_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector rhs_channels,Channel_Vector velocity_star_channels)
    {Run(allocator,blocks,rhs_channels,velocity_star_channels);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector rhs_channels,Channel_Vector& velocity_star_channels) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto mpm_rhs_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Node_Saturated) for(int v=0;v<d;++v){
                    allocator.template Get_Array<Struct_type,T>(rhs_channels(v))(offset)=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(v))(offset);
                                                                                                        
                // Log::cout<<"vi: "<<allocator.template Get_Array<Struct_type,T>(velocity_channels(v))(offset)
                // <<", dt/mass: "<< dt/mass(offset) <<", fi: "<< allocator.template Get_Array<Struct_type,T>(f_channels(v))(offset)
                // <<", v*i: "<<allocator.template Get_Array<Struct_type,T>(velocity_star_channels(v))(offset)<<std::endl;
            }}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,mpm_rhs_helper);        
    }

};
}
#endif

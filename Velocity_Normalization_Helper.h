//!#####################################################################
//! \file Velocity_Normalization_Helper.h
//!#####################################################################
// Class Velocity_Normalization_Helper
//######################################################################
#ifndef __Velocity_Normalization_Helper__
#define __Velocity_Normalization_Helper__
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include "MPM_Flags.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Velocity_Normalization_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Velocity_Normalization_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& velocity_channels)
    {Run(allocator,blocks,velocity_channels);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Vector<T Struct_type::*,2>& velocity_channels) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0); auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto v0=allocator.template Get_Array<Struct_type,T>(velocity_channels(0)); auto v1=allocator.template Get_Array<Struct_type,T>(velocity_channels(1));

        auto velocity_normalization_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Saturated) { const T mass_inverse=1/mass(offset);
                    // Log::cout<<"p: "<<v0(offset)<<","<<v1(offset)<<", m: "<<mass(offset)<<", mass_inverse: "<<mass_inverse<<", v: "<<v0(offset)*mass_inverse<<","<<v1(offset)*mass_inverse<<std::endl;
                    v0(offset)*=mass_inverse; v1(offset)*=mass_inverse;}
            }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,velocity_normalization_helper);
    }

    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Vector<T Struct_type::*,3>& velocity_channels) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0); auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto v0=allocator.template Get_Array<Struct_type,T>(velocity_channels(0)); auto v1=allocator.template Get_Array<Struct_type,T>(velocity_channels(1)); auto v2=allocator.template Get_Array<Struct_type,T>(velocity_channels(2));

        auto velocity_normalization_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Saturated) { const T mass_inverse=1/mass(offset);
                    v0(offset)*=mass_inverse; v1(offset)*=mass_inverse; v2(offset)*=mass_inverse;}
            }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,velocity_normalization_helper);
    }

    void Min_Max(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& velocity_channels) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto min_max_helper=[&](uint64_t offset)
        {
            T min_mass=FLT_MAX, max_mass=-FLT_MAX;
            T min_v=FLT_MAX,    max_v=-FLT_MAX;
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Saturated) {
                    T m=mass(offset);
                    if(m>max_mass) max_mass=m;
                    if(m<min_mass) min_mass=m;
                    for(int v=0;v<d;++v){
                    T velocity=abs(allocator.template Get_Array<Struct_type,T>(velocity_channels(v))(offset));
                    if(velocity<min_v) min_v=velocity;
                    if(velocity>max_v) max_v=velocity;
                    Vector<int,d> index(Flag_array_mask::LinearToCoord(offset));}}}
            Log::cout<<"min v: "<<min_v<<", max v: "<<max_v<<std::endl;
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,min_max_helper);
    }
};
}
#endif

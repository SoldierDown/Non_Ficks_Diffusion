//!#####################################################################
//! \file Explicit_Force_Helper.h
//!#####################################################################
// Class Explicit_Force_Helper
//######################################################################
#ifndef __Explicit_Force_Helper__
#define __Explicit_Force_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>


namespace Nova{
template<class Struct_type,class T,int d>
class Explicit_Force_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Explicit_Force_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                  Channel_Vector f_channels,Channel_Vector velocity_channels,Channel_Vector& velocity_star_channels,const T dt)
    {Run(allocator,blocks,f_channels,velocity_channels,velocity_star_channels,dt);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Vector<T Struct_type::*,2> f_channels, Vector<T Struct_type::*,2> velocity_channels, Vector<T Struct_type::*,2>& velocity_star_channels,const T dt) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto v0=allocator.template Get_Const_Array<Struct_type,T>(velocity_channels(0)); auto v1=allocator.template Get_Const_Array<Struct_type,T>(velocity_channels(1)); 
        auto vs0=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(0)); auto vs1=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(1));
        auto f0=allocator.template Get_Const_Array<Struct_type,T>(f_channels(0)); auto f1=allocator.template Get_Const_Array<Struct_type,T>(f_channels(1));
        auto explicit_velocity_update_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){vs0(offset)=v0(offset)+dt/mass(offset)*f0(offset); vs1(offset)=v1(offset)+dt/mass(offset)*f1(offset);}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,explicit_velocity_update_helper);        
    }

    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Vector<T Struct_type::*,3> f_channels, Vector<T Struct_type::*,3> velocity_channels, Vector<T Struct_type::*,3>& velocity_star_channels,const T dt) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0);  auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto v0=allocator.template Get_Const_Array<Struct_type,T>(velocity_channels(0)); auto v1=allocator.template Get_Const_Array<Struct_type,T>(velocity_channels(1)); auto v2=allocator.template Get_Const_Array<Struct_type,T>(velocity_channels(2));
        auto vs0=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(0)); auto vs1=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(1)); auto vs2=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(2));
        auto f0=allocator.template Get_Const_Array<Struct_type,T>(f_channels(0)); auto f1=allocator.template Get_Const_Array<Struct_type,T>(f_channels(1)); auto f2=allocator.template Get_Const_Array<Struct_type,T>(f_channels(2));
        auto explicit_velocity_update_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){vs0(offset)=v0(offset)+dt/mass(offset)*f0(offset); vs1(offset)=v1(offset)+dt/mass(offset)*f1(offset); vs2(offset)=v2(offset)+dt/mass(offset)*f2(offset);}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,explicit_velocity_update_helper);        
    }

};
}
#endif

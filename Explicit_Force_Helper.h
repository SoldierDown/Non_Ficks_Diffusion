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
                                  Channel_Vector& f_channels,Channel_Vector velocity_channels,Channel_Vector& velocity_star_channels,
                                  T Struct_type::* mass_channel,T Struct_type::* valid_nodes_channel,const T dt)
    {Run(allocator,blocks,f_channels,velocity_channels,velocity_star_channels,mass_channel,valid_nodes_channel,dt);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Channel_Vector& f_channels,Channel_Vector velocity_channels,Channel_Vector& velocity_star_channels,
             T Struct_type::* mass_channel,T Struct_type::* valid_nodes_channel,const T dt) const
    {
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(mass_channel);
        auto valid_nodes=allocator.template Get_Const_Array<Struct_type,T>(valid_nodes_channel);
        auto explicit_velocity_update_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(valid_nodes(offset)>(T)0.) 
                for(int v=0;v<d;++v){
                    allocator.template Get_Array<Struct_type,T>(velocity_star_channels(v))(offset)=allocator.template Get_Array<Struct_type,T>(velocity_channels(v))(offset)
                                                                                                        +dt/mass(offset)*allocator.template Get_Array<Struct_type,T>(f_channels(v))(offset);
                    // Log::cout<<"v*i: "<<allocator.template Get_Array<Struct_type,T>(velocity_star_channels(v))(offset)
                    //     <<", vi: "<<allocator.template Get_Array<Struct_type,T>(velocity_channels(v))(offset)<<std::endl;
                    // T velocity=allocator.template Get_Array<Struct_type,T>(f_channels(v))(offset);
                    // Vector<int,d> index(Flag_array_mask::LinearToCoord(offset));
            }}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,explicit_velocity_update_helper);        
    }

};
}
#endif

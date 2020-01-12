//!#####################################################################
//! \file Grid_Based_Collision_Helper.h
//!#####################################################################
// Class Grid_Based_Collision_Helper
//######################################################################
#ifndef __Grid_Based_Collision_Helper__
#define __Grid_Based_Collision_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Arrays/Array.h>
#include <nova/Tools/Vectors/Vector.h>

#include "MPM_Plane_Barrier.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Grid_Based_Collision_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using TV                    = Vector<T,d>;
    using T_Barrier             = MPM_Plane_Barrier<T,d>;
  public:
    Grid_Based_Collision_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& velocity_star_channels,
                T Struct_type::* collide_nodes_channel,unsigned Struct_type::* flags_channel,Array<T_Barrier> barriers)
    {Run(allocator,blocks,velocity_star_channels,collide_nodes_channel,flags_channel,barriers);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& velocity_star_channels,
                T Struct_type::* collide_nodes_channel,unsigned Struct_type::* flags_channel,Array<T_Barrier> barriers) const
    {
        auto collide_nodes=allocator.template Get_Const_Array<Struct_type,T>(collide_nodes_channel);
        auto valid_nodes=allocator.template Get_Const_Array<Struct_type,unsigned>(flags_channel);
        auto grid_based_collision_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                for(int id=0;id<barriers.size();++id){
                    TV normal_vector=barriers(id).normal; T mu=barriers(id).mu;
                        if((valid_nodes(offset)&Node_Saturated)
                    &&(collide_nodes(offset)==((T)id+(T)1.)||collide_nodes(offset)==(T)3.)){
                    TV vel=TV();
                    for(int v=0;v<d;++v) vel(v)=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(v))(offset);
                    T projection=vel.Dot_Product(normal_vector);
                    if(projection<(T)0.){
                        vel-=normal_vector*projection;
                        if(-projection*mu<vel.Norm()) vel+=vel.Normalized()*projection*mu;
                        else vel=TV();}
                    for (int v=0;v<d;++v) allocator.template Get_Array<Struct_type,T>(velocity_star_channels(v))(offset)=vel(v);}}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,grid_based_collision_helper);
    }
};
}
#endif






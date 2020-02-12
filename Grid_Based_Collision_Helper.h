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
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using TV                    = Vector<T,d>;
    using T_Barrier             = MPM_Plane_Barrier<T,d>;
  public:
    Grid_Based_Collision_Helper(const Grid<T,d>& grid,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& velocity_star_channels,Array<T_Barrier> barriers)
    {Run(grid,allocator,blocks,velocity_star_channels,barriers);}

    void Run(const Grid<T,2>& grid,SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,2>& velocity_star_channels, Array<MPM_Plane_Barrier<T,2> > barriers) const
    {
        using TV=Vector<T,2>; using T_INDEX=Vector<int,2>;
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto vs0=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(0)); auto vs1=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(1));
        auto grid_based_collision_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){
                    T_INDEX index(Flag_array_mask::LinearToCoord(offset)); TV cell_location=grid.Center(index);
                for(int id=0;id<barriers.size();++id){
                    const TV& normal_vector=barriers(id).normal; const T& mu=barriers(id).mu;
                    const T result=(cell_location-barriers(id).surface).Dot_Product(barriers(id).normal);
                    if(result<(T)0.){
                        TV vel=TV({vs0(offset),vs1(offset)});
                        T projection=vel.Dot_Product(normal_vector);
                        if(projection<(T)0.){
                            vel-=normal_vector*projection;
                            if(-projection*mu<vel.Norm()) vel+=vel.Normalized()*projection*mu;
                            else vel=TV();}
                        vs0(offset)=vel(0); vs1(offset)=vel(1);
                    }}}
                }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,grid_based_collision_helper);
    }
    void Run(const Grid<T,3>& grid,SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,3>& velocity_star_channels, Array<MPM_Plane_Barrier<T,3> > barriers) const
    {
        using TV=Vector<T,3>; using T_INDEX=Vector<int,3>;
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto vs0=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(0));  auto vs1=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(1)); auto vs2=allocator.template Get_Array<Struct_type,T>(velocity_star_channels(2));
        auto grid_based_collision_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){
                    T_INDEX index(Flag_array_mask::LinearToCoord(offset)); TV cell_location=grid.Center(index);
                for(int id=0;id<barriers.size();++id){
                    const TV& normal_vector=barriers(id).normal; const T& mu=barriers(id).mu;
                    const T result=(cell_location-barriers(id).surface).Dot_Product(barriers(id).normal);
                    if(result<(T)0.){
                        TV vel=TV({vs0(offset),vs1(offset),vs2(offset)});
                        T projection=vel.Dot_Product(normal_vector);
                        if(projection<(T)0.){
                            vel-=normal_vector*projection;
                            if(-projection*mu<vel.Norm()) vel+=vel.Normalized()*projection*mu;
                            else vel=TV();}
                        vs0(offset)=vel(0); vs1(offset)=vel(1); vs2(offset)=vel(2);
                    }}}
                }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,grid_based_collision_helper);
    }
};
}
#endif






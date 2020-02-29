//!#####################################################################
//! \file Levelset_Initializer.h
//!#####################################################################
// Class Levelset_Initializer
//######################################################################
#ifndef __Levelset_Initializer__
#define __Levelset_Initializer__
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include "Sphere_Levelset.h"
namespace Nova{
template<class Struct_type,class T,int d>
class Levelset_Initializer
{
    using TV                    = Vector<T,d>;    
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Levelset_Initializer(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,Sphere_Levelset<T,d>& levelset)
    {Run(hierarchy,allocator,blocks,levelset_channel,levelset);}

    void Run(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* levelset_channel,Sphere_Levelset<T,d>& levelset) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto data=allocator.template Get_Array<Struct_type,T>(levelset_channel);
        const Grid<T,d>& grid=hierarchy.Lattice(0);
        auto levelset_initializer=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior) {
                    Vector<int,d> index(Flag_array_mask::LinearToCoord(offset));
                    data(offset)=levelset.Distance(grid.Center(index));}}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,levelset_initializer);
    }



    // void Min_Max(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& velocity_channels) const
    // {
    //     auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0);
    //     auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
    //     auto min_max_helper=[&](uint64_t offset)
    //     {
    //         T min_mass=FLT_MAX, max_mass=-FLT_MAX;
    //         T min_v=FLT_MAX,    max_v=-FLT_MAX;
    //         for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
    //             if(flags(offset)&Cell_Type_Interior) {
    //                 T m=mass(offset);
    //                 if(m>max_mass) max_mass=m;
    //                 if(m<min_mass) min_mass=m;
    //                 for(int v=0;v<d;++v){
    //                 T velocity=abs(allocator.template Get_Array<Struct_type,T>(velocity_channels(v))(offset));
    //                 if(velocity<min_v) min_v=velocity;
    //                 if(velocity>max_v) max_v=velocity;
    //                 Vector<int,d> index(Flag_array_mask::LinearToCoord(offset));}}}
    //         Log::cout<<"min v: "<<min_v<<", max v: "<<max_v<<std::endl;
    //     };
    //     SPGrid_Computations::Run_Parallel_Blocks(blocks,min_max_helper);
    // }
};
}
#endif

//!#####################################################################
//! \file Grid_Saturation_Initialization_Helper.h
//!#####################################################################
// Class Grid_Saturation_Initialization_Helper
//######################################################################
#ifndef __Grid_Saturation_Initialization_Helper__
#define __Grid_Saturation_Initialization_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Grid_Saturation_Initialization_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    // using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Grid_Saturation_Initialization_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                            T Struct_type::* saturation_channel,T Struct_type::* void_mass_fluid_channel,const T volume)
    {Run(allocator,blocks,saturation_channel,void_mass_fluid_channel,volume);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* saturation_channel,T Struct_type::* void_mass_fluid_channel,const T volume) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto saturation=allocator.template Get_Array<Struct_type,T>(saturation_channel); 
        auto void_fluid_mass=allocator.template Get_Array<Struct_type,T>(void_mass_fluid_channel);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        auto grid_saturation_initialization_helper=[&](uint64_t offset)
        {
            const T fluid_density=(T)1.; const T unit_mass_fluid=fluid_density*volume;
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Dirichlet){ void_fluid_mass(offset)=unit_mass_fluid; saturation(offset)=unit_mass_fluid;}
                else if(flags(offset)&Cell_Type_Interior) void_fluid_mass(offset)=unit_mass_fluid;}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,grid_saturation_initialization_helper);
    }

};
}
#endif

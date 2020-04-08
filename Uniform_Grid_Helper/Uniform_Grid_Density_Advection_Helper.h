//!#####################################################################
//! \file Uniform_Grid_Density_Advection_Helper.h
//!#####################################################################
// Class Uniform_Grid_Density_Advection_Helper
//######################################################################
#ifndef __Uniform_Grid_Density_Advection_Helper__
#define __Uniform_Grid_Density_Advection_Helper__

#include "Uniform_Grid_Backtrace.h"
#include <nova/Dynamics/Hierarchy/Interpolation/Grid_Hierarchy_Interpolation.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include "Uniform_Grid_Interpolation_Helper.h"
namespace Nova{
template<class Struct_type,class T,int d>
class Uniform_Grid_Density_Advection_Helper
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Hierarchy_Interpolation   = Grid_Hierarchy_Interpolation<Struct_type,T,d>;

    enum {number_of_nodes_per_cell  = Topology_Helper::number_of_nodes_per_cell};

  public:
    // Checked
    Uniform_Grid_Density_Advection_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
                                    Channel_Vector& cell_velocity_channels,T Struct_type::* cell_channel,T Struct_type::* result_channel,
                                    uint64_t* nodes_of_cell_offsets,const int level,const T dt,const unsigned mask)
    {Run(hierarchy,blocks,cell_velocity_channels,cell_channel,result_channel,nodes_of_cell_offsets,level,dt,mask);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
            Channel_Vector& cell_velocity_channels,T Struct_type::* cell_channel,T Struct_type::* result_channel,
            uint64_t* nodes_of_cell_offsets,const int level,const T dt,const unsigned mask) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto result=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(result_channel);

        auto density_advection_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&mask){TV velocity;
                    for(int axis=0;axis<d;++axis) velocity(axis)=hierarchy.Channel(level,cell_velocity_channels(axis))(offset);
                    // backtrace
                    TV dX=-velocity*dt;
                    uint64_t new_offset=offset;int new_level=level;
                    TV weights=Uniform_Grid_Backtrace_Helper<Struct_type,T,d>::Uniform_Grid_Cell_Backtrace(hierarchy,new_level,index,new_offset,dX);
                    // Log::cout<<"velocity: "<<velocity<<", weights: "<<weights<<std::endl;
                    T value=Uniform_Grid_Interpolation_Helper<Struct_type,T,d>::Uniform_Grid_Cell_Interpolation_Helper(hierarchy,nodes_of_cell_offsets,new_level,new_offset,weights,cell_channel);
                    // Log::cout<<"velocity: "<<velocity<<", dt: "<<dt<<", dX:"<<dX<<", weights: "<<weights<<", value: "<<value<<std::endl;
                    result(offset)=value;}
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,density_advection_helper);
    }
};
}
#endif

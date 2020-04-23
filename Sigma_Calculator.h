//!#####################################################################
//! \file SF_M_Calculator.h
//!#####################################################################
// Class SF_M_Calculator
//######################################################################
#ifndef __SF_M_Calculator__
#define __SF_M_Calculator__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Tools/Utilities/Constants.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Sigma_Calculator
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Sigma_Calculator(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& vT_channels,T Struct_type::* sigma_channel,
                    const T delta)
    {Run(allocator,blocks,vT_channels,sigma_channel,delta);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& vT_channels,T Struct_type::* sigma_channel,
                    const T delta) const
    {
        auto sigma_data=allocator.template Get_Array<Struct_type,T>(sigma_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto sigma_calculator=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior)  {T a=(T)0.; T b=(T)0.; for(int axis=0;axis<d;++axis){ 
                const T axis_cell_vT_value=allocator.template Get_Const_Array<Struct_type,T>(vT_channels(axis))(offset);
                a+=pow(axis_cell_vT_value,4); b+=pow(axis_cell_vT_value,2);}
                sigma_data(offset)=(T)1.-delta*((T)1.-a/b);}

        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,sigma_calculator);
    }

};
}
#endif

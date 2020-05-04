//!#####################################################################
//! \file Add_Random_Term.h
//!#####################################################################
// Class Add_Random_Term
//######################################################################
#ifndef __Add_Random_Term__
#define __Add_Random_Term__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Tools/Utilities/Constants.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Add_Random_Term
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Add_Random_Term(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                                    T Struct_type::* density_channel,T Struct_type::* density_backup_channel,
                                    const T random_factor,const T random_value)
    {Run(allocator,blocks,density_channel,density_backup_channel,random_factor,random_value);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,                                    
            T Struct_type::* density_channel,T Struct_type::* density_backup_channel,
            const T random_factor,const T random_value) const
    {
        auto density_data=allocator.template Get_Array<Struct_type,T>(density_channel);
        auto density_backup_data=allocator.template Get_Const_Array<Struct_type,T>(density_backup_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto Add_Random_Term=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){ 
                    const T density_backup_value=density_backup_data(offset);
                    if(density_backup_value>(T)0. && density_backup_value<=(T).5) 
                        density_data(offset)+=random_factor*density_backup_value*((T)1.-density_backup_value)*(random_value-(T).5);}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,Add_Random_Term);
    }

};
}
#endif

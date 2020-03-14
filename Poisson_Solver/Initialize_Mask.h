//!#####################################################################
//! \file Initialize_Mask.h
//!#####################################################################
// Class Initialize_Mask
//######################################################################
#ifndef __Initialize_Mask__
#define __Initialize_Mask__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Base_struct_type,class Multigrid_struct_type,class T,int d>
class Initialize_Mask
{
    using Base_flags_type                       = typename Base_struct_type::Flags_type;
    using Hierarchy_Multigrid                   = Grid_Hierarchy<Multigrid_struct_type,T,d>;
    using Multigrid_Hierarchy                   = Array<Hierarchy_Multigrid*>;
    using Base_allocator_type                   = SPGrid::SPGrid_Allocator<Base_struct_type,d>;
    using Base_flag_array_mask                  = typename Base_allocator_type::template Array_mask<unsigned>;
    using Multigrid_allocator_type              = SPGrid::SPGrid_Allocator<Multigrid_struct_type,d>;
    using Multigrid_flag_array_mask             = typename Multigrid_allocator_type::template Array_mask<unsigned>;

  public:
    Initialize_Mask(Base_allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                    Multigrid_Hierarchy& multigrid_hierarchy,const int spgrid_level,const unsigned mask)
    {Run(allocator,blocks,multigrid_hierarchy,spgrid_level,mask);}

    void Run(Base_allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Multigrid_Hierarchy& multigrid_hierarchy,const int spgrid_level,const unsigned mask) const
    {
        const int mg_levels=(int)multigrid_hierarchy.size();
        auto flags=allocator.template Get_Const_Array<Base_struct_type,unsigned>(&Base_struct_type::flags);

        auto initialize_mask=[&](uint64_t offset)
        {
            for(int e=0;e<Base_flag_array_mask::elements_per_block;++e,offset+=sizeof(Base_flags_type))
                if(flags(offset)&mask){uint64_t current_offset=offset;
                    uint64_t translated_offset=Multigrid_flag_array_mask::template Translate_Linear_Offset<Base_flag_array_mask>(offset);
                    for(int mg_level=0;mg_level<mg_levels;++mg_level){
                        if(mg_level<=spgrid_level){
                            multigrid_hierarchy(mg_level)->Channel(spgrid_level,&Multigrid_struct_type::flags)(translated_offset)&=~Cell_Type_Interior;
                            multigrid_hierarchy(mg_level)->Channel(spgrid_level,&Multigrid_struct_type::flags)(translated_offset)|=mask;}
                        else{current_offset=Base_flag_array_mask::DownsampleOffset(current_offset);
                            translated_offset=Multigrid_flag_array_mask::template Translate_Linear_Offset<Base_flag_array_mask>(current_offset);
                            multigrid_hierarchy(mg_level)->Channel(mg_level,&Multigrid_struct_type::flags)(translated_offset)&=~Cell_Type_Interior;
                            multigrid_hierarchy(mg_level)->Channel(mg_level,&Multigrid_struct_type::flags)(translated_offset)|=mask;}}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,initialize_mask);
    }
};
}
#endif

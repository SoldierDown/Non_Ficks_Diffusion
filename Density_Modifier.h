//!#####################################################################
//! \file Density_Modifier.h
//!#####################################################################
// Class Density_Modifier
//######################################################################
#ifndef __Density_Modifier__
#define __Density_Modifier__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Density_Modifier
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Implicit_Object_Array     = Array<Implicit_Object<T,d>*>;

  public:
    Density_Modifier(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel,
                     Implicit_Object_Array& sources,const int level)
    {Run(hierarchy,blocks,density_channel,sources,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel,
             Implicit_Object_Array& sources,const int level) const
    {
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto density=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(density_channel);

        auto density_modifier=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior){const TV X=hierarchy.Lattice(level).Center(index);
                    for(size_t i=0;i<sources.size();++i) if(sources(i)->Inside(X)) density(offset)=(T)1.;}
                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,density_modifier);
    }

};
}
#endif

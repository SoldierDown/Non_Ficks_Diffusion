//!#####################################################################
//! \file Multiply_Inverse_Diagonal.h
//!#####################################################################
// Class Multiply_Inverse_Diagonal
//######################################################################
#ifndef __Multiply_Inverse_Diagonal__
#define __Multiply_Inverse_Diagonal__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Utilities.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Multiply_Inverse_Diagonal
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Multiply_Inverse_Diagonal(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* source_channel,
                              T Struct_type::* destination_channel,const T coeff1,const unsigned mask,const int level)
    {Run(hierarchy,blocks,source_channel,destination_channel,coeff1,mask,level);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* source_channel,
             T Struct_type::* destination_channel,const T coeff1,const unsigned mask,const int level) const
    {
        auto data_source=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(source_channel);
        auto data_destination=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(destination_channel);
        auto flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);

        const T laplace_scale_uniform   = Nova_Utilities::Sqr<T>(hierarchy.Lattice(level).one_over_dX[0]);
        const T diagonal                = (d==2)?((T)1.+(T)4.*coeff1*laplace_scale_uniform):((T)1.+(T)6.*coeff1*laplace_scale_uniform);
        const T one_over_diagonal       = (T)1./diagonal;

        auto multiply_inverse_diagonal=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) data_destination(offset)=data_source(offset)*one_over_diagonal;
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,multiply_inverse_diagonal);
    }
};
}
#endif

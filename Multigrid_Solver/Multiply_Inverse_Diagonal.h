//!#####################################################################
//! \file Multiply_Inverse_Diagonal.h
//!#####################################################################
// Class Multiply_Inverse_Diagonal
//######################################################################
#ifndef __Multiply_Inverse_Diagonal__
#define __Multiply_Inverse_Diagonal__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

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
                              T Struct_type::* destination_channel,const unsigned mask,const int level,
                                bool FICKS,const T dt,const T diff_coeff,const T Fc,const T tau)
    {Run(hierarchy,blocks,source_channel,destination_channel,mask,level,FICKS,dt,diff_coeff,Fc,tau);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* source_channel,
             T Struct_type::* destination_channel,const unsigned mask,const int level,
                bool FICKS,const T dt,const T diff_coeff,const T Fc,const T tau)
                              
    {
        auto flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto data_source=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(source_channel);
        auto data_destination=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(destination_channel);
        
        T face_areas[d];
        for(int axis=0;axis<d;++axis){face_areas[axis]=(T)1.;
            for(int other_axis=0;other_axis<d;++other_axis) if(other_axis!=axis) face_areas[axis]*=hierarchy.Lattice(level).dX[other_axis];}
        
        double scaling_factor=hierarchy.Lattice(0).one_over_dX.Product();
        const T laplace_scale_uniform   = scaling_factor*face_areas[0]*hierarchy.Lattice(level).one_over_dX[0];
        const T a=dt*diff_coeff*laplace_scale_uniform; const T twod_a_plus_one=2.*d*a+1.;
        const T coeff1=dt*diff_coeff*(Fc*tau+dt)*laplace_scale_uniform/(dt+tau);
        const T diagonal                = (FICKS)?twod_a_plus_one:((T)1.+(T)2.*d*coeff1);
        const T one_over_diagonal       = (T)1./diagonal;

        auto multiply_inverse_diagonal=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&mask) data_destination(offset)=data_source(offset)*one_over_diagonal;
            }
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,multiply_inverse_diagonal);
    }
};
}
#endif

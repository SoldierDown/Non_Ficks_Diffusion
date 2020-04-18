//!#####################################################################
//! \file Psi_Evaluation_Helper.h
//!#####################################################################
// Class Psi_Evaluation_Helper
//######################################################################
#ifndef __Psi_Evaluation_Helper__
#define __Psi_Evaluation_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>
#include <nova/Tools/Utilities/Constants.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Psi_Evaluation_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Psi_Evaluation_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            Channel_Vector& dSdX_channels,T Struct_type::* psi_channel)
    {Run(allocator,blocks,dSdX_channels,psi_channel);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,2>& dSdX_channels,T Struct_type::* psi_channel) const
    {
        auto dSdx_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(0));
        auto dSdy_data=allocator.template Get_Const_Array<Struct_type,T>(dSdX_channels(1));
        auto psi_data=allocator.template Get_Array<Struct_type,T>(psi_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto psi_evaluation_helper_2d=[&](uint64_t offset)
        {
            const T half_pi=(T).5*pi;
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&Cell_Type_Interior){
                if(dSdx_data(offset)==(T)0.){
                    if(dSdy_data(offset)>(T)0.) psi_data(offset)=half_pi;
                    else if(dSdy_data(offset)<(T)0.) psi_data(offset)=-half_pi;}
                else{ const T atan_value=atan(dSdy_data(offset)/dSdx_data(offset));
                    if(dSdx_data(offset)>(T)0.){
                        if(dSdy_data(offset)>(T)0.) psi_data(offset)=atan_value;
                        else if(dSdy_data(offset)<(T)0.) psi_data(offset)=two_pi+atan_value;}
                    else psi_data(offset)=pi+atan_value;}} 
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,psi_evaluation_helper_2d);
    }
    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Vector<T Struct_type::*,3>& dSdX_channels,T Struct_type::* psi_channel) const
    {

    }

};
}
#endif

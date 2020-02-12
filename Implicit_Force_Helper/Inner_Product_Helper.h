//!#####################################################################
//! \file Inner_Product_Helper.h
//!#####################################################################
// Class Inner_Product_Helper
//######################################################################
#ifndef __Inner_Product_Helper__
#define __Inner_Product_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Inner_Product_Helper
{
    using TV                        = Vector<T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator            = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
  public:
    Inner_Product_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Channel_Vector& channels1,const Channel_Vector& channels2,
                          T& result,const unsigned mask)
    {Run(allocator,blocks,channels1,channels2,result,mask);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Vector<T Struct_type::*,2>& channels0,const Vector<T Struct_type::*,2>& channels1,
                          T& result,const unsigned mask) const
    {
        result=(T)0.;
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto c00=allocator.template Get_Const_Array<Struct_type,T>(channels0(0)); auto c01=allocator.template Get_Const_Array<Struct_type,T>(channels0(1));
        auto c10=allocator.template Get_Const_Array<Struct_type,T>(channels1(0)); auto c11=allocator.template Get_Const_Array<Struct_type,T>(channels1(1));
        T temp_result=(T)0.;
#pragma omp parallel for reduction(+:temp_result)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) {
                  TV vec1({c00(offset),c01(offset)}),vec2({c10(offset)*mass(offset),c11(offset)*mass(offset)});
                  temp_result+=vec1.Dot_Product(vec2);}}
        result+=temp_result;
    }

    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,const Vector<T Struct_type::*,3>& channels0,const Vector<T Struct_type::*,3>& channels1,
                          T& result,const unsigned mask) const
    {
        result=(T)0.;
        auto mass=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch0);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto c00=allocator.template Get_Const_Array<Struct_type,T>(channels0(0));
        auto c01=allocator.template Get_Const_Array<Struct_type,T>(channels0(1));
        auto c02=allocator.template Get_Const_Array<Struct_type,T>(channels0(2));
        auto c10=allocator.template Get_Const_Array<Struct_type,T>(channels1(0));
        auto c11=allocator.template Get_Const_Array<Struct_type,T>(channels1(1));
        auto c12=allocator.template Get_Const_Array<Struct_type,T>(channels1(2));
        T temp_result=(T)0.;
#pragma omp parallel for reduction(+:temp_result)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) { 
                  // TV tmp_vec; for(int v=0;v<d;++v) tmp_vec(v)=allocator.template Get_Const_Array<Struct_type,T>(channels1(v))(offset);
                  // Log::cout<<"mass: "<<node_mass<<", tmp_vec: "<<tmp_vec<<std::endl;
                  TV vec1({c00(offset),c01(offset),c02(offset)}), vec2({c10(offset)*mass(offset),c11(offset)*mass(offset),c12(offset)*mass(offset)});
                  temp_result+=vec1.Dot_Product(vec2);
                  // for(int v=0;v<d;++v) temp_result+=node_mass*allocator.template Get_Const_Array<Struct_type,T>(channels1(v))(offset)*allocator.template Get_Const_Array<Struct_type,T>(channels2(v))(offset);
          }}

        result+=temp_result;
    }
};
}
#endif

//!#####################################################################
//! \file Convergence_Norm_Helper.h
//!#####################################################################
// Class Convergence_Norm_Helper
//######################################################################
#ifndef __Convergence_Norm_Helper__
#define __Convergence_Norm_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <algorithm>

namespace Nova{
template<class Struct_type,class T,int d>
class Convergence_Norm_Helper
{
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator            = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;

  public:
    Convergence_Norm_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            Channel_Vector channels,T& result,const unsigned mask,const int threads)
    {Run(allocator,blocks,channels,result,mask,threads);}

//     void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
//              Vector<T Struct_type::*,2> channels,T& result,const unsigned mask) const
//     {
//         auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
//         auto c0=allocator.template Get_Const_Array<Struct_type,T>(channels(0));
//         auto c1=allocator.template Get_Const_Array<Struct_type,T>(channels(1));
//         T max_value=0;

// #pragma omp parallel for reduction(max:max_value)
//         for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
//             for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
//                 if(flags(offset)&mask) max_value=std::max(max_value,Nova_Utilities::Sqr(c0(offset))+Nova_Utilities::Sqr(c1(offset)));}

//         result=std::sqrt(std::max(result,max_value));
//     }
//     void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
//              Vector<T Struct_type::*,3> channels,T& result,const unsigned mask) const
//     {
//         auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
//         auto c0=allocator.template Get_Const_Array<Struct_type,T>(channels(0));
//         auto c1=allocator.template Get_Const_Array<Struct_type,T>(channels(1));
//         auto c2=allocator.template Get_Const_Array<Struct_type,T>(channels(2));
//         T max_value=0;

// #pragma omp parallel for reduction(max:max_value)
//         for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
//             for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
//                 if(flags(offset)&mask) max_value=std::max(max_value,Nova_Utilities::Sqr(c0(offset))+Nova_Utilities::Sqr(c1(offset))+Nova_Utilities::Sqr(c2(offset)));}

//         result=std::sqrt(std::max(result,max_value));
//     }
    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
             Channel_Vector channels,T& result,const unsigned mask,const int threads) const
    {
        // auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        // result=(T)0.;

        // auto convergence_norm_helper=[&](uint64_t offset,T& result){
        //   for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
        //       for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
        //           if(flags(offset)&mask) for(int v=0;v<d;++v) result+=Nova_Utilities::Sqr(allocator.template Get_Const_Array<Struct_type,T>(channels(v))(offset));}};
        
        // for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
        //     uint64_t offset=iterator.Offset();
        //     convergence_norm_helper(offset,result);}

        // result=std::sqrt(result);
//         Array<T> result_per_thread(threads);
// #pragma omp parallel for
//         for(unsigned i=0;i<simulated_particles.size();++i){const int tid=omp_get_thread_num();
//           T& r=result_per_thread(tid); const int id=simulated_particles(i);
//           r=std::max(r,particles(id).V.Norm_Squared());}
//         T result=(T)0.;
//         for(int tid=0;tid<threads;++tid) result=std::max(result,result_per_thread(tid));
//         std::cout<<"max v: "<<std::sqrt(result)<<std::endl;
//         return std::sqrt(result);

        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        T max_value=0;

#pragma omp parallel for reduction(max:max_value)
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) for(int v=0;v<d;++v) max_value=std::max(max_value,std::fabs(allocator.template Get_Const_Array<Struct_type,T>(channels(v))(offset)));}
        result=std::max(result,max_value);

    }
};
}
#endif

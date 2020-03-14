//!#####################################################################
//! \file Multigrid_Refinement.h
//!#####################################################################
// Class Multigrid_Refinement
//######################################################################
#ifndef __Multigrid_Refinement__
#define __Multigrid_Refinement__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Projection/Ghost_Value_Accumulate.h>
#include <nova/Dynamics/Hierarchy/Projection/Ghost_Value_Propagate.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/SPGrid/Tools/SPGrid_Multiple_Allocator_Masked_Plus_Equals_Helper.h>
#include "Prolongation_Stencil_Helper.h"
#include "Restriction_Stencil_Helper.h"

namespace Nova{

extern int number_of_threads;

template<class Struct_type,class T,int d>
class Multigrid_Refinement
{
    using TV                                = Vector<T,d>;
    using T_INDEX                           = Vector<int,d>;
    using Hierarchy                         = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type                    = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                   = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper                   = Grid_Topology_Helper<Flag_array_mask>;
    enum {nodes_per_cell                    = Topology_Helper::number_of_nodes_per_cell};
    enum {restriction_stencil_denominator   = (d==2)?16:64};

  public:
    Multigrid_Refinement() {}
    ~Multigrid_Refinement() {}

    static void Restrict(Hierarchy& fine_hierarchy,Hierarchy& coarse_hierarchy,T Struct_type::* fine_data_channel,
                         T Struct_type::* coarse_data_channel,const Vector<int,2>& finest_active_level)
    {
        // grabbing levels
        const int fine_level=finest_active_level(0),coarse_level=finest_active_level(1),levels=fine_hierarchy.Levels();
        assert(levels==coarse_hierarchy.Levels());
        assert(fine_level>=0 && coarse_level<levels);

        // clear coarse data
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(coarse_hierarchy.Allocator(level),coarse_hierarchy.Blocks(level),coarse_data_channel);

        // restrict from fine to coarse for finest level
        {
            auto coarse_data=coarse_hierarchy.Allocator(coarse_level).template Get_Array<Struct_type,T>(coarse_data_channel);
            auto fine_data=fine_hierarchy.Allocator(fine_level).template Get_Const_Array<Struct_type,T>(fine_data_channel);
            auto coarse_flags=coarse_hierarchy.Allocator(coarse_level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
            const T scale=(T)1./(T)restriction_stencil_denominator;

            Restriction_Stencil_Helper<Struct_type,T,d> restriction_stencil_helper((T*)coarse_data.Get_Data_Ptr(),(T*)fine_data.Get_Data_Ptr(),
                (unsigned*)coarse_flags.Get_Data_Ptr(),coarse_hierarchy.Blocks(coarse_level).first,coarse_hierarchy.Blocks(coarse_level).second,
                scale,(unsigned)(Cell_Type_Interior|Cell_Type_Ghost));
            if(number_of_threads) restriction_stencil_helper.Run_Parallel(number_of_threads);
            else restriction_stencil_helper.Run();
        }

        // accumulate within coarse
        for(int level=coarse_level;level<levels-1;++level)
            Ghost_Value_Accumulate<Struct_type,T,d>(coarse_hierarchy,coarse_hierarchy.Blocks(level+1),
                                                    coarse_data_channel,coarse_data_channel,level+1);

        // copy fine to coarse for all higher levels
        for(int level=levels-1;level>=coarse_level;--level){
            SPGrid::Multiple_Allocator_Masked_Plus_Equals_Helper<Struct_type,T,unsigned,d> helper(
                coarse_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),coarse_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),
                fine_hierarchy.Blocks(level),coarse_data_channel,fine_data_channel,coarse_data_channel,&Struct_type::flags,(unsigned)Cell_Type_Interior);
            if(number_of_threads) helper.Run_Parallel(number_of_threads);
            else helper.Run();}
    }

    static void Prolongate(Hierarchy& fine_hierarchy,Hierarchy& coarse_hierarchy,T Struct_type::* fine_data_channel,
                           T Struct_type::* coarse_data_channel,const Vector<int,2>& finest_active_level)
    {
        // auxiliary stuff
        uint64_t nodes_of_cell_offsets[nodes_per_cell];
        Topology_Helper::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
        Vector<uint64_t,d> parity_masks,negative_axis_vector_offsets;
        for(int v=0;v<d;++v){parity_masks(v)=Topology_Helper::Axis_Vector_Offset(v);
            negative_axis_vector_offsets(v)=Topology_Helper::Negative_Axis_Vector_Offset(v);}

        // grabbing levels
        const int fine_level=finest_active_level(0),coarse_level=finest_active_level(1),levels=fine_hierarchy.Levels();
        assert(levels==coarse_hierarchy.Levels());
        assert(fine_level>=0 && coarse_level<levels);

        // clear fine data
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(fine_hierarchy.Allocator(level),fine_hierarchy.Blocks(level),fine_data_channel);

        // copy coarse to fine for all higher levels
        for(int level=levels-1;level>=coarse_level;--level){
            SPGrid::Multiple_Allocator_Masked_Plus_Equals_Helper<Struct_type,T,unsigned,d> helper(
                fine_hierarchy.Allocator(level),coarse_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),
                fine_hierarchy.Blocks(level),fine_data_channel,coarse_data_channel,fine_data_channel,&Struct_type::flags,(unsigned)Cell_Type_Interior);
            if(number_of_threads) helper.Run_Parallel(number_of_threads);
            else helper.Run();}

        // propagate within coarse
        for(int level=levels-2;level>coarse_level;--level)
            Ghost_Value_Propagate<Struct_type,T,d>(coarse_hierarchy,coarse_hierarchy.Blocks(level),coarse_data_channel,coarse_data_channel,level);

        // prolongate from coarse to fine for finest level
        {
            auto coarse_data=coarse_hierarchy.Allocator(coarse_level).template Get_Array<Struct_type,T>(coarse_data_channel);
            auto fine_data=fine_hierarchy.Allocator(fine_level).template Get_Const_Array<Struct_type,T>(fine_data_channel);
            auto fine_flags=fine_hierarchy.Allocator(fine_level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

            Prolongation_Stencil_Helper<Struct_type,T,d> prolongation_stencil_helper(
                (T*)fine_data.Get_Data_Ptr(),(T*)coarse_data.Get_Data_Ptr(),(unsigned*)fine_flags.Get_Data_Ptr(),
                fine_hierarchy.Blocks(fine_level).first,fine_hierarchy.Blocks(fine_level).second,(unsigned)Cell_Type_Interior);
            if(number_of_threads) prolongation_stencil_helper.Run_Parallel(number_of_threads);
            else prolongation_stencil_helper.Run();
        }
    }
};
}
#endif

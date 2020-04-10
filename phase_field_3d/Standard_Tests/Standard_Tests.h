//!#####################################################################
//! \file Standard_Tests.h
//!#####################################################################
// Class Standard_Tests
//######################################################################
#ifndef __Standard_Tests__
#define __Standard_Tests__

#include <nova/Geometry/Implicit_Objects/Box_Implicit_Object.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include "../../PF_Data.h"
#include "../../PF_Example.h"
#include "../../Rasterizers/Adaptive_Sphere_Rasterizer.h"
#include "../../Rasterizers/Randomized_Rasterizer.h"
#include "../../Sphere_Implicit_Object.h"

namespace Nova{
template<class T,int d>
class Standard_Tests: public PF_Example<T,d>
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Struct_type               = PC_Data<T>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Base                      = PF_Example<T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    using Base::output_directory; using Base::test_number;using Base::counts;using Base::levels;using Base::domain_walls;using Base::hierarchy;using Base::rasterizer;
    using Base::cfl;    using Base::sources;    using Base::density_channel;
    using Base::omega;  using Base::FICKS;
    using Base::explicit_diffusion;
    /****************************
     * example explanation:
     *
     * 1. Simple inflow/outflow.
     ****************************/

    Standard_Tests()
        :Base()
    {}

//######################################################################
    void Parse_Options() override
    {
        Base::Parse_Options();
        output_directory=(explicit_diffusion?(FICKS?"Phase_Change_F_":"Phase_Change_NF_"):(FICKS?"Implicit_Phase_Change_F_":"Implicit_Phase_Change_NF_"))+std::to_string(d)+"d_"+std::to_string(omega)+"branches_Resolution_"+std::to_string(counts(0))+"x"+std::to_string(counts(1));
        output_directory="3Qt";
        for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=false;
        TV min_corner,max_corner=TV(7.68);
        hierarchy=new Hierarchy(counts,Range<T,d>(min_corner,max_corner),levels);
    }
//######################################################################
    void Initialize_Rasterizer() override
    {
        rasterizer=new Randomized_Rasterizer<Struct_type,T,d>(*hierarchy);
    }
//######################################################################
    void Initialize_State() override
    {
        // clear density channel
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);

        for(int level=0;level<levels;++level){auto blocks=hierarchy->Blocks(level);
            auto block_size=hierarchy->Allocator(level).Block_Size();
            auto data=hierarchy->Allocator(level).template Get_Array<Struct_type,T>(density_channel);
            auto flags=hierarchy->Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

            for(unsigned block=0;block<blocks.second;++block){uint64_t offset=blocks.first[block];
                Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
                T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));

                for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                    const T_INDEX index=base_index+range_iterator.Index();
                    if(flags(offset)&Cell_Type_Interior && sources(0)->Inside(hierarchy->Lattice(level).Center(index))) data(offset)=(T)1.;
                    range_iterator.Next();}}}
    }
//######################################################################
    void Initialize_Sources() override
    {
        const T radius=0.03*4;
        const TV center=TV(3.84);
        Implicit_Object<T,d>* obj=new Sphere_Implicit_Object<T,d>(center,radius);
        sources.Append(obj);
    }
//######################################################################
};
}
#endif

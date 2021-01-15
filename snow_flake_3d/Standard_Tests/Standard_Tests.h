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
#include "../../SF_Data.h"
#include "../../SF_Example.h"
#include "../../Rasterizers/Adaptive_Sphere_Rasterizer.h"
#include "../../Rasterizers/Randomized_Rasterizer.h"
#include "../../Sphere_Implicit_Object.h"

namespace Nova{
template<class T,int d>
class Standard_Tests: public SF_Example<T,d>
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Struct_type               = SF_Data<T>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Base                      = SF_Example<T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    using Base::output_directory; using Base::test_number;using Base::counts;using Base::levels;using Base::domain_walls;using Base::hierarchy;using Base::rasterizer;
    using Base::cfl;    using Base::density_sources; using Base::velocity_sources;    using Base::density_channel; using Base::T_channel;
    using Base::omega;  using Base::FICKS; using Base::const_density_value; using Base::K; using Base::cell_width;
    using Base::T0;
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
        output_directory=(explicit_diffusion?(FICKS?"Snow_Flake_F_":"Snow_Flake_NF_"):(FICKS?"Implicit_Snow_Flake_F_":"Implicit_Snow_Flake_NF_"))+std::to_string(d)+"d_"
                        +std::to_string(omega)+"branches_Resolution_"+std::to_string(counts(0))+"x"+std::to_string(counts(1))+"x"+std::to_string(counts(2))+"_case"+std::to_string(test_number);
        for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=false;
        TV min_corner,max_corner=TV({cell_width*counts(0),cell_width*counts(1),cell_width*counts(2)});
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
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(0),hierarchy->Blocks(0),density_channel);

        auto blocks=hierarchy->Blocks(0);
        auto block_size=hierarchy->Allocator(0).Block_Size();
        auto density_data=hierarchy->Allocator(0).template Get_Array<Struct_type,T>(density_channel);
        auto T_data=hierarchy->Allocator(0).template Get_Array<Struct_type,T>(T_channel);
        auto flags=hierarchy->Allocator(0).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        for(unsigned block=0;block<blocks.second;++block){uint64_t offset=blocks.first[block];
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior && density_sources(0)->Inside(hierarchy->Lattice(0).Center(index))) density_data(offset)=const_density_value;
                T_data(offset)=-(T)1./K;
                range_iterator.Next();}}
    }
//######################################################################
    void Initialize_Sources() override
    {
        const T radius=(T)2.*cell_width;
        const TV center=TV({(T).5*cell_width*counts(0),(T)0.,(T).5*cell_width*counts(2)});
        Implicit_Object<T,d>* obj=new Sphere_Implicit_Object<T,d>(center,radius);
        density_sources.Append(obj);
    }
//######################################################################
};
}
#endif
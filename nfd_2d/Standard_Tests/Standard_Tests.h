//!#####################################################################
//! \file Standard_Tests.h
//!#####################################################################
// Class Standard_Tests
//######################################################################
#ifndef __Standard_Tests__
#define __Standard_Tests__

#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/Tools/Utilities/Range.h>
#include "../../MPM_Example.h"

namespace Nova{
template<class T,int d>
class Standard_Tests: public MPM_Example<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = MPM_Example<T,d>;
    using T_Particle        = MPM_Particle<T,d>;
    using Struct_type       = MPM_Data<T>;
    using Hierarchy         = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::output_directory;using Base::test_number;using Base::particles;using Base::parse_args;using Base::counts;using Base::domain;

    /****************************
     * example explanation:
     *
     * 1. Block falling on a mound.
     ****************************/

    Standard_Tests()
        :Base()
    {}

//######################################################################
    void Parse_Options() override
    {
        Base::Parse_Options();
        output_directory="nfd_2d_"+std::to_string(test_number);

        domain.min_corner=TV();domain.max_corner=TV(1);
    }
//######################################################################
    void Initialize_Particles() override
    {
        Log::Scope scope("Initialize_Particles");

        const T mass_density=(T)2.;
        const int number_of_particles=200;
        const Range<T,d> block(TV({.3,.7}),TV({.7,.9}));
        const T block_area=block.Size();
        const T area_per_particle=block_area/number_of_particles;

        Random_Numbers<T> random;
        random.Set_Seed(0);

        particles.resize(number_of_particles);
        for(int i=0;i<number_of_particles;++i){
            particles(i).X=random.Get_Uniform_Vector(block);
            particles(i).mass=mass_density*area_per_particle;}
    }
//######################################################################
};
}
#endif

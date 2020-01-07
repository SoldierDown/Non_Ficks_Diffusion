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

        //Base::simulated_particles.Clear();
        particles.Clear();

        const T mass_density=(T)2.;
        const int number_of_particles=2000;
        const Range<T,d> block(TV({.4,.3}),TV({.6,.5}));
        // const T block_area=block.Size();
        const T block_area=block.Area();
        const T area_per_particle=block_area/number_of_particles;
        // std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
        const T E=(T)40.,nu=(T).2;

        Random_Numbers<T> random;
        random.Set_Seed(0);
        for(int i=0;i<number_of_particles;++i){
            T_Particle p;
            p.X=random.Get_Uniform_Vector(block);
            p.V(1)=(T)-1.;
            p.mass=mass_density*area_per_particle;
            p.constitutive_model.Compute_Lame_Parameters(E,nu);
            p.constitutive_model.plastic=false;
            p.constitutive_model.stretching_yield=(T)1.005;
            p.constitutive_model.compression_yield=(T)0.985;
            p.constitutive_model.hardening_factor=(T)7.;
            particles.Append(p);
            // Base::simulated_particles.Append(i);
        }
        // Log::cout<<"simulated size: "<<Base::simulated_particles.size()<<", particles: "<<particles.size()<<std::endl;        
    }
//######################################################################
};
}
#endif

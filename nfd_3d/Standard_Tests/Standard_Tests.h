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
    using T_Barrier         = MPM_Plane_Barrier<T,d>;
    using Base              = MPM_Example<T,d>;
    using T_Particle        = MPM_Particle<T,d>;
    using Struct_type       = MPM_Data<T>;
    using Hierarchy         = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::output_directory;using Base::test_number;using Base::particles;using Base::parse_args;using Base::counts;using Base::domain;
    using Base::barriers; using Base::FICKS;
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
        output_directory=std::to_string(d)+"d_"+(FICKS?"F_":"NF_");

        domain.min_corner=TV();domain.max_corner=TV(1);
    }
//######################################################################
    void Initialize_Particles(int test_case) override
    {
        Log::Scope scope("Initialize_Particles");
        switch (test_case)
        {
        case 1:{
            Random_Numbers<T> random;
            random.Set_Seed(0);
            // T_Barrier ground(0.,TV({0.,1.,0.}),TV({0.,.1,0.}));
            // barriers.Append(ground);
            // T_Barrier ceiling(0.,TV({0.,-1.,0.}),TV({0.,.9,0.}));
            // Base::barriers.Append(ceiling);
            // T_Barrier left_wall(0.,TV({1.,0.,0.}),TV({.1,0.,0.}));
            // barriers.Append(left_wall);
            // T_Barrier right_wall(0.,TV({-1.,0.,0.}),TV({.9,0.,0.}));
            // barriers.Append(right_wall);
            // T_Barrier front_wall(0.,TV({0.,0.,1.}),TV({0.,0.,.1}));
            // barriers.Append(front_wall);
            // T_Barrier back_wall(0.,TV({0.,0.,-1.}),TV({0.,0.,.9}));
            // barriers.Append(back_wall);
            {
                const T solid_density=(T)10.;
                const T fluid_density=(T)1.;
                const int number_of_particles=40000;
                const Range<T,d> block(TV({.45,.45,.49}),TV({.55,.55,.51}));
                const T block_area=block.Area();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                const T E=(T)100.,nu=(T).4;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Uniform_Vector(block);
                    p.V=TV();
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.eta=(T).01;
                    p.constitutive_model.plastic=false;
                    p.saturation=(T)0.;
                    p.volume_fraction_0=(T).7;
                    p.mass_solid=solid_density*area_per_particle*((T)1.-p.volume_fraction_0);
                    p.mass_fluid=fluid_density*p.saturation*area_per_particle*p.volume_fraction_0;
                    p.mass=p.mass_solid+p.mass_fluid;
                    particles.Append(p);
                }  
            }
        }
    break;

    default:
        break;
    }}
};
//######################################################################
}

#endif

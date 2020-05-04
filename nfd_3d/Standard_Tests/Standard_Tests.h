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
    using Base::barriers;using Base::FICKS;using Base::E;using Base::nu;using Base::eta;using Base::Fc;using Base::tau;using Base::diff_coeff;
    using Base::mg_levels; using Base::np; using Base::nbw; using Base::fluid_source;
    using Base::min_dt; using Base::max_dt;
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
        output_directory=std::to_string(d)+"d_"+(FICKS?"F_":"NF_")+"E_"+std::to_string(E)+"_nu_"+std::to_string(nu)+"_diff_"+std::to_string(diff_coeff)
                                +"_eta_"+std::to_string(eta)+"_Fc_"+std::to_string(Fc)+"_tau_"+std::to_string(tau)+"_Resolution_"+std::to_string(counts(0))+"_"+std::to_string(mg_levels)+"_levels_"+std::to_string(np)+"_particles_dt_"
                                +std::to_string(min_dt)+"_"+std::to_string(max_dt)+"_"+std::to_string(nbw)+"cw";

        domain.min_corner=TV();domain.max_corner=TV(11);
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
                const int number_of_particles=np;//40000*std::pow(counts(0)/(T)256.,3);
                const Range<T,d> block(TV({5.4,5.4,5.49}),TV({5.6,5.6,5.51}));
                const T block_area=block.Area();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Uniform_Vector(block);
                    p.V=TV();
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.eta=eta;
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
        // example 19: fluid
        case 19:{
            Random_Numbers<T> random;
            random.Set_Seed(0);
            // T_Barrier ceiling(0.,TV({0.,-1.,0.}),TV({0.,.9,0.}));
            // Base::barriers.Append(ceiling);
            T_Barrier ground(0.,TV({0.,1.,0.}),TV({0.,1.,0.}));
            barriers.Append(ground);
            T_Barrier left_wall(0.,TV({1.,0.,0.}),TV({1,0.,0.}));
            Base::barriers.Append(left_wall);
            T_Barrier right_wall(0.,TV({-1.,0.,0.}),TV({9,0.,0.}));
            Base::barriers.Append(right_wall);
            T_Barrier front_wall(0.,TV({0.,0.,1.}),TV({0.,0.,1.}));
            Base::barriers.Append(front_wall);
            T_Barrier back_wall(0.,TV({0.,0.,-1.}),TV({0.,0.,9.}));
            Base::barriers.Append(back_wall);
            
            {
                const T mass_density=(T)2.;
                const int number_of_particles=200000;
                const Sphere<T,d> ball(TV({5.5,5.5,5.5}),1.5);
                const T block_area=ball.Size();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                const T E=(T)30.,nu=(T).4;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Vector_In_Sphere(ball)+ball.center;
                    p.V(0)=(T)0.;
                    p.V(1)=(T)-20.;
                    p.mass=mass_density*area_per_particle;
                    p.mass_solid=p.mass;
                    p.mass_fluid=(T)0.;
                    p.eos=true;
                    
                    p.bulk_modulus=(T)1.;
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.plastic=false;

                    p.constitutive_model.eta=(T)0.;
                    p.saturation=(T)0.;
                    p.volume_fraction_0=(T)1.;

                    particles.Append(p);
                }  
            }
            
        }break;
        // example 23: hydrogel falling
        case 23:{
            Random_Numbers<T> random;
            random.Set_Seed(0);
            
            fluid_source.min_corner=TV({4.75,4.75,4.75});
            fluid_source.max_corner=TV({5.25,5.25,5.25});

            T_Barrier ground(0.,TV({0.,1.,0.}),TV({0.,1.,0.}));
            barriers.Append(ground);
            T_Barrier left_wall(0.,TV({1.,0.,0.}),TV({1,0.,0.}));
            Base::barriers.Append(left_wall);
            T_Barrier right_wall(0.,TV({-1.,0.,0.}),TV({10,0.,0.}));
            Base::barriers.Append(right_wall);
            T_Barrier front_wall(0.,TV({0.,0.,1.}),TV({0.,0.,1.}));
            Base::barriers.Append(front_wall);
            T_Barrier back_wall(0.,TV({0.,0.,-1.}),TV({0.,0.,10.}));
            Base::barriers.Append(back_wall);          

            // hydrogel  
            {
                const T solid_density=(T)10.;
                const T fluid_density=(T)1.;
                const int number_of_particles=80000;
                const Range<T,d> block(TV({4.5,5.3,4.5}),TV({5.5,6.3,5.5}));
                const T block_area=block.Area();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.eos=false;
                    p.X=random.Get_Uniform_Vector(block);
                    p.V=TV({0.,-10.,0});
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.eta=(T)eta;
                    p.constitutive_model.plastic=false;
                    p.saturation=(T)0.;
                    p.volume_fraction_0=(T).7;
                    p.mass_solid=solid_density*area_per_particle*((T)1.-p.volume_fraction_0);
                    p.mass_fluid=fluid_density*p.saturation*area_per_particle*p.volume_fraction_0;
                    p.mass=p.mass_solid+p.mass_fluid;
                    particles.Append(p);
                }  
            }

            // fluid
            {
                const T solid_density=(T)10.;
                const T fluid_density=(T)1.;
                const int number_of_particles=400000;
                const Range<T,d> block(TV({1.,1.,1.}),TV({10.,5.,10.}));
                const T block_area=block.Area();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.valid=true; 
                    p.X=random.Get_Uniform_Vector(block);
                    p.V=TV();
                    p.mass=fluid_density*area_per_particle;
                    p.mass_fluid=p.mass;
                    p.mass_solid=(T)0.;
                    p.volume=(T)0.;
                    p.scp=Matrix<T,3>();
                    p.eos_scp=Matrix<T,3>();
                    // EOS fluid particle
                    p.eos=true;
                    p.density=1;
                    p.bulk_modulus=(T)1.;
                    p.gamma=(T)7;
                    p.saturation=(T)1.;
                    p.volume_fraction_0=(T)1.;
                    particles.Append(p);
                }  
            }
            
        }break;
    default:
        break;
    }}
};
//######################################################################
}

#endif

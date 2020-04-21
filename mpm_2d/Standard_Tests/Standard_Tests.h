//!#####################################################################
//! \file Standard_Tests.h
//!#####################################################################
// Class Standard_Tests
//######################################################################
#ifndef __Standard_Tests__
#define __Standard_Tests__

#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/Tools/Utilities/Range.h>
#include <nova/Geometry/Basic_Geometry/Sphere.h>
#include "../../MPM_Example.h"


namespace Nova{
template<class T,int d>
class Standard_Tests: public MPM_Example<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = MPM_Example<T,d>;
    using T_Particle        = MPM_Particle<T,d>;
    using T_Barrier         = MPM_Plane_Barrier<T,d>;
    using Struct_type       = MPM_Data<T>;
    using Hierarchy         = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::output_directory;using Base::test_number;using Base::particles;using Base::parse_args;using Base::counts;using Base::domain;
    using Base::barriers;
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
        output_directory="mpm_2d_"+std::to_string(test_number);

        domain.min_corner=TV();domain.max_corner=TV(1);
    }
//######################################################################
    void Initialize_Particles(int test_case) override
    {
        Log::Scope scope("Initialize_Particles");
        particles.Clear();
        Log::cout<<"Pi: "<<pi<<std::endl;
        switch (test_case)
        {
        // Test: one particle
        case 1:{
            Random_Numbers<T> random;
            random.Set_Seed(0);
            {
                const T mass_density=(T)2.;
                const int number_of_particles=1;
                const Sphere<T,d> ball(TV({.2,.5}),.1);
                const T block_area=ball.Size();
                const T area_per_particle=block_area/number_of_particles;
                Log::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                const T E=(T)40.,nu=(T).4;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Vector_In_Sphere(ball)+ball.center;
                    Log::cout<<"position: "<<p.X<<std::endl;
                    p.V=TV({.5,0.});
                    p.mass=mass_density*area_per_particle;
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.plastic=true;
                    p.constitutive_model.stretching_yield=(T)1.01;
                    p.constitutive_model.compression_yield=(T)0.95;
                    p.constitutive_model.hardening_factor=(T)10.;
                    particles.Append(p);
                }
            }
            

        }break;


        case 15:{
            const T mass_density=(T)2.;
            const int number_of_particles=2000;
            const Range<T,d> block(TV({.4,.4}),TV({.6,.6}));
            const T block_area=block.Area();
            const T area_per_particle=block_area/number_of_particles;
            // std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
            const T E=(T)40.,nu=(T).2;

            Random_Numbers<T> random;
            random.Set_Seed(0);
            for(int i=0;i<number_of_particles;++i){
                T_Particle p;
                p.X=random.Get_Uniform_Vector(block);
                p.V=TV({-2.,.1});
                p.mass=mass_density*area_per_particle;
                p.constitutive_model.Compute_Lame_Parameters(E,nu);
                p.constitutive_model.plastic=true;
                p.constitutive_model.stretching_yield=(T)1.01;
                p.constitutive_model.compression_yield=(T)0.95;
                p.constitutive_model.hardening_factor=(T)10.;
                particles.Append(p);
            }  

        }break;
        // snow ball hitting the ground
        case 16:{
            Random_Numbers<T> random;
            random.Set_Seed(0);
            
            {
                const T mass_density=(T)2.;
                const int number_of_particles=2000;
                const Sphere<T,d> ball(TV({.5,.5}),.1);
                const T block_area=ball.Size();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                const T E=(T)50.,nu=(T).45;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Vector_In_Sphere(ball)+ball.center;
                    p.V(1)=(T)-1.5;
                    p.mass=mass_density*area_per_particle;
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.plastic=true;
                    p.constitutive_model.stretching_yield=(T)1.005;
                    p.constitutive_model.compression_yield=(T)0.975;
                    p.constitutive_model.hardening_factor=(T)10.;
                    particles.Append(p);
                }  
            }
            {
                const T mass_density=(T)1.;
                const int number_of_particles=1000;
                const Range<T,d> block(TV({0.,.1}),TV({1.,.125}));
                const T block_area=block.Area();
                const T area_per_particle=block_area/number_of_particles;
                const T E=(T)5., nu=(T).4;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Uniform_Vector(block);
                    p.mass=mass_density*area_per_particle;
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.plastic=true;
                    p.constitutive_model.stretching_yield=(T)1.001;
                    p.constitutive_model.compression_yield=(T)0.975;
                    p.constitutive_model.hardening_factor=(T)5.;
                    particles.Append(p);
                }              
            }

        }break;
        // example 17: snow ball hitting the wall
        case 17:{
            Random_Numbers<T> random;
            random.Set_Seed(0);
            T_Barrier wall(0.,TV({1.,0.}),TV({.1,0.}));
            barriers.Append(wall);
            T_Barrier ground(0.,TV({0.,1.}),TV({0.,.1}));
            barriers.Append(ground);
            {
                const T mass_density=(T)2.;
                const int number_of_particles=2000;
                const Sphere<T,d> ball(TV({.5,.5}),.1);
                const T block_area=ball.Size();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                const T E=(T)30.,nu=(T).4;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Vector_In_Sphere(ball)+ball.center;
                    p.V(0)=(T)-2.;
                    p.V(1)=(T).1;
                    p.mass=mass_density*area_per_particle;
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.plastic=true;
                    p.constitutive_model.stretching_yield=(T)1.01;
                    p.constitutive_model.compression_yield=(T)0.95;
                    p.constitutive_model.hardening_factor=(T)10.;
                    particles.Append(p);
                }  
            }
        } break;
        // example 18: two snow balls collide
        case 18:{
            Random_Numbers<T> random;
            random.Set_Seed(0);
            T_Barrier ground(0.,TV({0.,1.}),TV({0.,.1}));
            barriers.Append(ground);
            // T_Barrier ceiling(0.,TV({0.,-1.}),TV({0.,.9}));
            // Base::barriers.Append(ceiling);
            // T_Barrier left_wall(0.,TV({1.,0.}),TV({.1,0.}));
            // barriers.Append(left_wall);
            // T_Barrier right_wall(0.,TV({-1.,0.}),TV({.9,0.}));
            // barriers.Append(right_wall);
            
            {
                const T mass_density=(T)2.;
                const int number_of_particles=2000;
                const Sphere<T,d> ball(TV({.8,.5}),.1);
                const T block_area=ball.Size();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                const T E=(T)30.,nu=(T).4;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Vector_In_Sphere(ball)+ball.center;
                    p.V(0)=(T)-2.;
                    p.V(1)=(T).1;
                    p.mass=mass_density*area_per_particle;
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.plastic=true;
                    p.constitutive_model.stretching_yield=(T)1.01;
                    p.constitutive_model.compression_yield=(T)0.95;
                    p.constitutive_model.hardening_factor=(T)10.;
                    particles.Append(p);
                }  
            }

            {
                const T mass_density=(T)2.;
                const int number_of_particles=2000;
                const Sphere<T,d> ball(TV({.2,.5}),.1);
                const T block_area=ball.Size();
                const T area_per_particle=block_area/number_of_particles;
                std::cout<<"block area: "<<block_area<<", area per particle:"<<area_per_particle<<std::endl;
                const T E=(T)30.,nu=(T).4;
                for(int i=0;i<number_of_particles;++i){
                    T_Particle p;
                    p.X=random.Get_Vector_In_Sphere(ball)+ball.center;
                    p.V(0)=(T)2.;
                    p.V(1)=(T).1;
                    p.mass=mass_density*area_per_particle;
                    p.constitutive_model.Compute_Lame_Parameters(E,nu);
                    p.constitutive_model.plastic=true;
                    p.constitutive_model.stretching_yield=(T)1.01;
                    p.constitutive_model.compression_yield=(T)0.95;
                    p.constitutive_model.hardening_factor=(T)10.;
                    particles.Append(p);
                }  
            }
            
        }break;
        default:
            break;
        }

    }
//######################################################################
};
}
#endif

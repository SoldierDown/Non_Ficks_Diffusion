//!#####################################################################
//! \file MPM_Example.h
//!#####################################################################
// Class MPM_Example
//######################################################################
#ifndef __MPM_Example__
#define __MPM_Example__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Tools/Utilities/Constants.h>
#include <nova/Tools/Utilities/Example.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include <nova/Tools/Utilities/Utilities.h>
#include <nova/Tools/Arrays/Array.h>
#include "MPM_Data.h"
#include "Poisson_Data.h"
#include "MPM_Particle.h"
#include "MPM_Plane_Barrier.h"
#include "./Tools/Matrix_MXN.h"
#include "./Tools/Interval.h"
#include "./Tools/Influence_Iterator.h"
#include "./Tools/Cropped_Influence_Iterator.h"

namespace Nova{
template<class T>
static T N2(const T x)
{
    if(fabs(x)<(T).5) return (T).75-Nova_Utilities::Sqr(x);
    else if(fabs(x)<(T)1.5) return (T).5*Nova_Utilities::Sqr((T)1.5-fabs(x));
    else return (T)0.;
}
template<class T,int d>
inline T N2(const Vector<T,d>& X,const Vector<T,d> one_over_dX)
{
    T value=(T)1.;
    for(int axis=0;axis<d;++axis) value*=N2(X(axis)*one_over_dX(axis));
    return value;
}
template<class T>
static T dN2(const T x)
{
    T sign=x>=(T)0.?(T)1.:(T)-1.;
    if(fabs(x)<(T).5) return -(T)2.*x;
    else if(fabs(x)<(T)1.5) return x-(T)1.5*sign;
    else return (T)0.;  
}
template<class T,int d>
inline T dN2(const Vector<T,d> X,const int axis,const Vector<T,d> one_over_dX)
{
    T value=(T)1.;
    for(int v=0;v<d;++v) 
        if(v==axis) value*=one_over_dX(v)*dN2(X(v)*one_over_dX(v));
        else value*=N2(X(v)*one_over_dX(v));
    return value;
}
template<class T,int d>
inline Vector<T,d> dN2(const Vector<T,d> X,const Vector<T,d> one_over_dX)
{
    Vector<T,d> value=Vector<T,d>();
    for(int axis=0;axis<d;++axis) value(axis)=dN2(X,axis,one_over_dX);
    return value;
}
template<class T,int d>
static T N3(const T x)
{
    if(fabs(x)<(T)1.) return (T).5*Nova_Utilities::Cube(fabs(x))-Nova_Utilities::Sqr(x)+(T)two_thirds;
    else if(fabs(x)<(T)2.) return (T)-one_sixth*Nova_Utilities::Cube(fabs(x))+Nova_Utilities::Sqr(x)-(T)2.*fabs(x)+(T)four_thirds;
    else return (T)0.;
}
template<class T,int d>
inline T N3(const Vector<T,d> X,const Vector<T,d> one_over_dX)
{
    T value=(T)1.;
    for(int axis=0;axis<d;++axis) value*=N3(X(axis)*one_over_dX(axis));
    return value;
}
template<class T>
static T dN3(const T x)
{
    T sign=x>=(T)0.?(T)1.:(T)-1.;
    if(fabs(x)<(T)1.) return (T)1.5*Nova_Utilities::Sqr(x)*sign-(T)2.*x;
    else if(fabs(x)<(T)2.) return (T)-.5*Nova_Utilities::Sqr(x)*sign+(T)2.*x-(T)2.*sign;
    else return (T)0.;  
}
template<class T,int d>
inline T dN3(const Vector<T,d> X,const int axis,const Vector<T,d> one_over_dX)
{
    T value=(T)1.;
    for(int v=0;v<d;++v) 
        if(v==axis) value*=one_over_dX(v)*dN3(X(v)*one_over_dX(v));
        else value*=N3(X(v)*one_over_dX(v));
    return value;
}
template<class T,int d>
inline Vector<T,d> dN3(const Vector<T,d> X,const Vector<T,d> one_over_dX)
{
    Vector<T,d> value=Vector<T,d>();
    for(int axis=0;axis<d;++axis) value(axis)=dN3(X,axis,one_over_dX);
    return value;
}

template<class T> void
Vector_To_Flag(Vector<int,2> current_cell)
{
    using Struct_type               = MPM_Data<T>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,2>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    std::cout<<Flag_array_mask::Linear_Offset(current_cell(0),current_cell(1))<<std::endl;
}

template<class T> void
Vector_To_Flag(Vector<int,3> current_cell)
{
    using Struct_type               = MPM_Data<T>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,3>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    std::cout<<Flag_array_mask::Linear_Offset(current_cell(0),current_cell(1),current_cell(2))<<std::endl;
}

template <class T,int d> class MPM_Example;

template<class T>
class MPM_Example<T,2>: public Example<T,2>
{
    using T_Mat                         = Matrix<T,2>;
    using TV                            = Vector<T,2>;
    using Base                          = Example<T,2>;
    using T_INDEX                       = Vector<int,2>;
    using T_Particle                    = MPM_Particle<T,2>;
    using T_Barrier                     = MPM_Plane_Barrier<T,2>;
    using MPM_struct_type               = MPM_Data<T>;
    using Diff_struct_type              = Poisson_Data<T>;
    using MPM_flags_type                = typename MPM_struct_type::Flags_type;
    using Diff_flags_type               = typename Diff_struct_type::Flags_type;
    using MPM_allocator_type            = SPGrid::SPGrid_Allocator<MPM_struct_type,2>;
    using Diff_allocator_type           = SPGrid::SPGrid_Allocator<Diff_struct_type,2>;
    using MPM_flag_array_mask           = typename MPM_allocator_type::template Array_mask<unsigned>;
    using Diff_flag_array_mask          = typename Diff_allocator_type::template Array_mask<unsigned>;
    using MPM_Hierarchy                 = Grid_Hierarchy<MPM_struct_type,T,2>;
    using Diff_Hierarchy                = Grid_Hierarchy<Diff_struct_type,T,2>;
    using MPM_Channel_Vector            = Vector<T MPM_struct_type::*,2>;
    using Diff_Channel_Vector           = Vector<T Diff_struct_type::*,2>;
    using T_Influence_Iterator          = Influence_Iterator<T,2,T_INDEX>;
    using T_Cropped_Influence_Iterator  = Cropped_Influence_Iterator<T,2,T_INDEX>;
  public:
    using Base::frame_title;using Base::output_directory;using Base::parse_args;using Base::first_frame;

    T explicit_force_rt=0.;
    int explicit_force_cnt=0;
    T apply_force_rt=0.;
    int apply_force_cnt=0;
    T ras_rt=0.;
    int ras_cnt=0;
    T ras_vx_rt=0.;
    int ras_vx_cnt=0;
    T update_x_v_rt=0.;
    int update_x_v_cnt=0;
    T flip;
    T cfl;
    T solver_tolerance;
    int solver_iterations;
    int levels,threads;
    T_INDEX counts;
    Range<T,2> domain;
    Range<T,2> bbox;
    Array<T_Particle> particles;
    Array<T_Barrier> barriers;
    Array<int> simulated_particles;
    Array<int> invalid_particles;
    Array<int> valid_grid_indices;
    Array<Array<int> > valid_grid_indices_thread;
    Array<Interval<int> > x_intervals;
    Matrix_MxN<Array<int> > particle_bins;
    TV gravity;

    MPM_Hierarchy *mpm_hierarchy;
    Diff_Hierarchy *diff_hierarchy;

    unsigned MPM_struct_type::* mpm_flags_channel;
    T MPM_struct_type::* mass_channel;
    MPM_Channel_Vector velocity_channels;
    MPM_Channel_Vector velocity_star_channels;
    MPM_Channel_Vector f_channels;



    // Hydrogel variables
    T E;
    T nu;
    T eta;
    T diff_coeff;
    T tau;
    T Fc;
    bool FICKS;
    bool explicit_diffusion;

    // Hydrogel channels
    unsigned Diff_struct_type::* diff_flags_channel;
    T Diff_struct_type::* saturation_channel;
    T Diff_struct_type::* lap_saturation_channel;
    T Diff_struct_type::* void_mass_fluid_channel;
    T Diff_struct_type::* volume_channel;
    T Diff_struct_type::* div_Qc_channel;

    T Diff_struct_type::* diff_rhs_channel;
    T Diff_struct_type::* diff_q_channel;
    T Diff_struct_type::* diff_s_channel;
    T Diff_struct_type::* diff_r_channel;
    T Diff_struct_type::* diff_z_channel;

    // Krylov solver channels
    MPM_Channel_Vector mpm_rhs_channels;
    MPM_Channel_Vector mpm_q_channels;
    MPM_Channel_Vector mpm_s_channels;
    MPM_Channel_Vector mpm_r_channels;
    MPM_Channel_Vector mpm_z_channels;

    MPM_Example();

    ~MPM_Example();


//######################################################################
    virtual void Initialize_Particles(int test_case)=0;
//######################################################################
    void Initialize_SPGrid();
    void Initialize();
    void Reset_Grid_Based_Variables();
    void Reset_Solver_Channels();
    void Populate_Simulated_Particles();
    void Update_Particle_Weights();
    void Group_Particles();
    void Rasterize();
    void Grid_Based_Collision(const bool detect_collision);
    void Update_Constitutive_Model_State();
    void Update_Particle_Velocities_And_Positions(const T dt);
    void Estimate_Particle_Volumes();
    void Ficks_Diffusion(T dt);
    void Non_Ficks_Diffusion(T dt);
    void Apply_Force(const T dt);
    void Apply_Explicit_Force(const T dt);
    void Grid_Based_Collison(const bool detect_collision);
    T    Max_Particle_Velocity() const;
    void Limit_Dt(T& dt,const T time) override;
    void Test();
    void Register_Options() override;
    void Parse_Options() override;
    void Read_Output_Files(const int frame);
    void Write_Output_Files(const int frame) const override;
    void Rasterize_Voxels();
  protected:
    void Compute_Bounding_Box(Range<T,2>& bbox);
//######################################################################
};

template<class T>
class MPM_Example<T,3>: public Example<T,3>
{
    using T_Mat                         = Matrix<T,3>;
    using TV                            = Vector<T,3>;
    using Base                          = Example<T,3>;
    using T_INDEX                       = Vector<int,3>;
    using T_Particle                    = MPM_Particle<T,3>;
    using T_Barrier                     = MPM_Plane_Barrier<T,3>;
    using MPM_struct_type               = MPM_Data<T>;
    using Diff_struct_type              = Poisson_Data<T>;
    using MPM_flags_type                = typename MPM_struct_type::Flags_type;
    using Diff_flags_type               = typename Diff_struct_type::Flags_type;
    using MPM_allocator_type            = SPGrid::SPGrid_Allocator<MPM_struct_type,3>;
    using Diff_allocator_type           = SPGrid::SPGrid_Allocator<Diff_struct_type,3>;
    using MPM_flag_array_mask           = typename MPM_allocator_type::template Array_mask<unsigned>;
    using Diff_flag_array_mask          = typename Diff_allocator_type::template Array_mask<unsigned>;
    using MPM_Hierarchy                 = Grid_Hierarchy<MPM_struct_type,T,3>;
    using Diff_Hierarchy                = Grid_Hierarchy<Diff_struct_type,T,3>;
    using MPM_Channel_Vector            = Vector<T MPM_struct_type::*,3>;
    using Diff_Channel_Vector           = Vector<T Diff_struct_type::*,3>;
    using T_Influence_Iterator          = Influence_Iterator<T,3,T_INDEX>;
    using T_Cropped_Influence_Iterator  = Cropped_Influence_Iterator<T,3,T_INDEX>;
  public:
    using Base::frame_title;using Base::output_directory;using Base::parse_args;using Base::first_frame;

    T explicit_force_rt=0.;
    int explicit_force_cnt=0;
    T apply_force_rt=0.;
    int apply_force_cnt=0;
    T ras_rt=0.;
    int ras_cnt=0;
    T ras_vx_rt=0.;
    int ras_vx_cnt=0;
    T update_x_v_rt=0.;
    int update_x_v_cnt=0;
    T flip;
    T cfl;
    T solver_tolerance;
    int solver_iterations;
    int levels,threads;
    T_INDEX counts;
    Range<T,3> domain;
    Range<T,3> bbox;
    Array<T_Particle> particles;
    Array<T_Barrier> barriers;
    Array<int> simulated_particles;
    Array<int> invalid_particles;
    Array<int> valid_grid_indices;
    Array<Array<int> > valid_grid_indices_thread;
    Array<Interval<int> > x_intervals;
    Matrix_MxN<Array<int> > particle_bins;
    TV gravity;

    MPM_Hierarchy *mpm_hierarchy;
    Diff_Hierarchy *diff_hierarchy;

    unsigned MPM_struct_type::* mpm_flags_channel;
    T MPM_struct_type::* mass_channel;
    MPM_Channel_Vector velocity_channels;
    MPM_Channel_Vector velocity_star_channels;
    MPM_Channel_Vector f_channels;


    // Hydrogel variables
    T E;
    T nu;
    T eta;
    T diff_coeff;
    T tau;
    T Fc;
    bool FICKS;
    bool explicit_diffusion;

    // Hydrogel channels
    unsigned Diff_struct_type::* diff_flags_channel;
    T Diff_struct_type::* saturation_channel;
    T Diff_struct_type::* lap_saturation_channel;
    T Diff_struct_type::* void_mass_fluid_channel;
    T Diff_struct_type::* volume_channel;
    T Diff_struct_type::* div_Qc_channel;

    T Diff_struct_type::* diff_rhs_channel;
    T Diff_struct_type::* diff_q_channel;
    T Diff_struct_type::* diff_s_channel;
    T Diff_struct_type::* diff_r_channel;
    T Diff_struct_type::* diff_z_channel;

    // Krylov solver channels
    MPM_Channel_Vector mpm_rhs_channels;
    MPM_Channel_Vector mpm_q_channels;
    MPM_Channel_Vector mpm_s_channels;
    MPM_Channel_Vector mpm_r_channels;
    MPM_Channel_Vector mpm_z_channels;

    MPM_Example();

    ~MPM_Example();


//######################################################################
    virtual void Initialize_Particles(int test_case)=0;
//######################################################################
    void Initialize_SPGrid();
    void Initialize();
    void Reset_Grid_Based_Variables();
    void Reset_Solver_Channels();
    void Populate_Simulated_Particles();
    void Update_Particle_Weights();
    void Group_Particles();
    void Rasterize();
    void Grid_Based_Collision(const bool detect_collision);
    void Update_Constitutive_Model_State();
    void Update_Particle_Velocities_And_Positions(const T dt);
    void Estimate_Particle_Volumes();
    void Ficks_Diffusion(T dt);
    void Non_Ficks_Diffusion(T dt);
    void Apply_Force(const T dt);
    void Apply_Explicit_Force(const T dt);
    void Grid_Based_Collison(const bool detect_collision);
    T    Max_Particle_Velocity() const;
    void Limit_Dt(T& dt,const T time) override;
    void Test();
    void Register_Options() override;
    void Parse_Options() override;
    void Read_Output_Files(const int frame);
    void Write_Output_Files(const int frame) const override;
    void Rasterize_Voxels();
  protected:
    void Compute_Bounding_Box(Range<T,3>& bbox);
//######################################################################
};


}
#endif

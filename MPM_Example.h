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
#include "MPM_Particle.h"
#include "MPM_Plane_Barrier.h"
#include "./Tools/Matrix_MXN.h"
#include "./Tools/Interval.h"
#include "./Tools/Cropped_Range_Interator.h"

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
Vector_To_Flag(Vector<int,2> current_node)
{
    using Struct_type               = MPM_Data<T>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,2>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    std::cout<<Flag_array_mask::Linear_Offset(current_node(0),current_node(1))<<std::endl;
}

template<class T> void
Vector_To_Flag(Vector<int,3> current_node)
{
    using Struct_type               = MPM_Data<T>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,3>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    std::cout<<Flag_array_mask::Linear_Offset(current_node(0),current_node(1),current_node(2))<<std::endl;
}

template<class T,int d>
class MPM_Example: public Example<T,d>
{
    using TV                        = Vector<T,d>;
    using Base                      = Example<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using T_Particle                = MPM_Particle<T,d>;
    using T_Barrier                 = MPM_Plane_Barrier<T,d>;
    using Struct_type               = MPM_Data<T>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using T_Range_Iterator          = Range_Iterator<d,T_INDEX>;
    using T_Cropped_Range_Iterator  = Cropped_Range_Iterator<d,T_INDEX>;

  public:
    using Base::frame_title;using Base::output_directory;using Base::parse_args;using Base::first_frame;

    T flip;
    T cfl;
    T solver_tolerance;
    int solver_iterations;
    int levels,threads;
    T_INDEX counts;
    Range<T,d> domain;
    Range<T,d> bbox;
    Array<T_Particle> particles;
    Array<T_Barrier> barriers;
    Array<int> simulated_particles;
    Array<int> invalid_particles;
    Array<int> valid_grid_indices;
    Array<Array<int> > valid_grid_indices_thread;
    Array<Interval<int> > x_intervals;
    Matrix_MxN<Array<int> > particle_bins;
    TV gravity;
    Hierarchy *hierarchy;

    

    unsigned Struct_type::* flags_channel;
    T Struct_type::* mass_channel;
    Channel_Vector velocity_channels;
    Channel_Vector velocity_star_channels;
    Channel_Vector f_channels;


    // Krylov solver channels
    Channel_Vector rhs_channels;
    Channel_Vector q_channels;
    Channel_Vector s_channels;
    Channel_Vector r_channels;
    Channel_Vector z_channels;

    T Struct_type::* collide_nodes_channel;

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
    void Update_Constitutive_Model_State();
    void Update_Particle_Velocities_And_Positions(const T dt);
    void Estimate_Particle_Volumes();
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
  protected:
    void Compute_Bounding_Box(Range<T,d>& bbox);
    void Rasterize_Voxels(const Range<T,d>& bbox);
//######################################################################
};
}
#endif

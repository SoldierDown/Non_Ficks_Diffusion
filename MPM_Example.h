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

namespace Nova{
template<class T,int d>
class MPM_Example: public Example<T,d>
{
    using TV                        = Vector<T,d>;
    using Base                      = Example<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using T_Particle                = MPM_Particle<T,d>;
    using Struct_type               = MPM_Data<T>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using T_Range_Iterator          = Range_Iterator<d,T_INDEX>;

  public:
    using Base::frame_title;using Base::output_directory;using Base::parse_args;using Base::first_frame;

    T flip;
    T cfl;
    int levels,threads;
    T_INDEX counts;
    Range<T,d> domain;
    Range<T,d> bbox;
    Array<T_Particle> particles;
    Array<int> simulated_particles;
    Array<int> invalid_particles;
    Array<int> valid_grid_indices;
    Array<Array<int> > valid_grid_indices_thread;
    TV gravity;
    Hierarchy *hierarchy;

    T Struct_type::* mass_channel;
    Channel_Vector velocity_channels;
    Channel_Vector velocity_star_channels;
    Channel_Vector f_channels;

    MPM_Example();

    ~MPM_Example();

    static T N(const T x)
    {
        if(fabs(x)<(T)1.) return (T).5*Nova_Utilities::Cube(fabs(x))-Nova_Utilities::Sqr(x)+(T)two_thirds;
        else if(fabs(x)<(T)2.) return (T)-one_sixth*Nova_Utilities::Cube(fabs(x))+Nova_Utilities::Sqr(x)-(T)2.*fabs(x)+(T)four_thirds;
        else return (T)0.;
    }

    inline T N(const TV& X)
    {
        const Grid<T,d>& grid=hierarchy->Lattice(0);T value=(T)1.;
        for(int axis=0;axis<d;++axis) value*=N(X(axis)*grid.one_over_dX(axis));
        return value;
    }

    static T dN(const T x)
    {
        int sign=x>=0?1:-1;
        if(fabs(x)<(T)1.) return (T)1.5*Nova_Utilities::Sqr(x)*sign-(T)2.*x;
        else if(fabs(x)<(T)2.) return (T)-.5*Nova_Utilities::Sqr(x)*sign+(T)2.*x-(T)2.*sign;
        else return (T)0.;  
    }

    inline T dN(const TV X, const int axis)
    {
        const Grid<T,d>& grid=hierarchy->Lattice(0);T value=(T)1.;
        for(int v=0;v<d;++v) 
            if(v==axis) value*=N(X(v)*grid.one_over_dX(v));
            else value*=N(X(v)*grid.one_over_dX(v));
        return value;
    }

    inline TV dN(const TV& X)
    {
        TV value=TV();
        for(int axis=0;axis<d;++axis) value(axis)=dN(X,axis);
        return value;
    }

//######################################################################
    virtual void Initialize_Particles()=0;
//######################################################################
    void Initialize_SPGrid();
    void Initialize();
    void Reset_Grid_Based_Variables();
    void Populate_Simulated_Particles();
    void Rasterize();
    void Update_Constitutive_Model_State();
    void Update_Particle_Velocities_And_Positions(const T dt);
    void Estimate_Particle_Volumes();
    void Apply_Force(const T dt);
    void Apply_Explicit_Force(const T dt);
    void Grid_Based_Collison();
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

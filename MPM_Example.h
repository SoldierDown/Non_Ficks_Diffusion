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
#include "Diff_Data.h"
#include "MPM_Particle.h"
#include "MPM_Plane_Barrier.h"
#include "./Tools/Matrix_MXN.h"
#include "./Tools/Interval.h"
#include "./Tools/Influence_Iterator.h"
#include "./Tools/Cropped_Influence_Iterator.h"

namespace Nova{
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
    using Diff_struct_type              = Diff_Data<T>;
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

    T nbw;
    T flip;
    T cfl;
    int levels,threads;
    int n_fluid,n_solid;
    T_INDEX counts;
    Range<T,2> domain;
    Range<T,2> bbox;
    Array<T_Particle> particles;
    Array<T_Barrier> barriers;
    Array<int> waiting_particles;
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

    Range<T,2> fluid_source;
    Random_Numbers<T> random;

    int mg_levels;
    int cg_max_iterations;
    int cg_restart_iterations;
    T min_dt;
    T max_dt;
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
    void Add_Fluid_Source(const T dt);
    int  Allocate_Particle(bool add_to_simulation=true);
    void Reset_Grid_Based_Variables();
    void Reset_Solver_Channels();
    void Populate_Simulated_Particles();
    void Update_Particle_Weights();
    void Group_Particles();
    void Rasterize();
    void Grid_Based_Collision(const bool detect_collision);
    void Update_Constitutive_Model_State();
    void Update_Particle_Velocities_And_Positions(const T dt);
    void Process_Waiting_Particles();
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
    using Diff_struct_type              = Diff_Data<T>;
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

    T nbw;
    T flip;
    T cfl;
    int iteration_counter;
    int levels,threads;
    int n_fluid,n_solid;
    T_INDEX counts;
    Range<T,3> domain;
    Range<T,3> bbox;
    Array<T_Particle> particles;
    Array<T_Barrier> barriers;
    Array<int> waiting_particles;
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

    Range<T,3> fluid_source;
    Random_Numbers<T> random;

    int mg_levels;
    int cg_max_iterations;
    int cg_restart_iterations;
    // Hydrogel variables
    T E;
    T nu;
    T eta;
    T diff_coeff;
    T tau;
    T Fc;
    bool FICKS;
    bool explicit_diffusion;
    T min_dt;
    T max_dt;
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
    void Add_Fluid_Source(const T dt);
    int  Allocate_Particle(bool add_to_simulation=true);
    void Reset_Grid_Based_Variables();
    void Reset_Solver_Channels();
    void Populate_Simulated_Particles();
    void Update_Particle_Weights();
    void Group_Particles();
    void Rasterize();
    void Grid_Based_Collision(const bool detect_collision);
    void Update_Constitutive_Model_State();
    void Update_Particle_Velocities_And_Positions(const T dt);
    void Process_Waiting_Particles();
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

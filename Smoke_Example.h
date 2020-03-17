//!#####################################################################
//! \file Smoke_Example.h
//!#####################################################################
// Class Smoke_Example
//######################################################################
#ifndef __Smoke_Example__
#define __Smoke_Example__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Tools/Utilities/Example.h>
#include "Poisson_Data.h"

namespace Nova{
template<class T,int d>
class Smoke_Example: public Example<T,d>
{
    using TV                        = Vector<T,d>;
    using Base                      = Example<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Struct_type               = Poisson_Data<T>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Rasterizer      = Hierarchical_Rasterizer<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    using Base::frame_title;using Base::output_directory;using Base::parse_args;using Base::first_frame;

    bool FICKS;
    bool explicit_diffusion;
    T diff_coeff;
    T Fc;
    T tau;
    T_INDEX counts;
    int levels,mg_levels,cg_iterations,cg_restart_iterations;
    T cfl,cg_tolerance;
    Hierarchy *hierarchy;
    Hierarchy_Rasterizer *rasterizer;

    T Struct_type::* density_channel;
    T Struct_type::* pressure_channel;
    Vector<T Struct_type::*,d> face_velocity_channels;
    Vector<T Struct_type::*,d> face_qc_channels;
    Vector<Vector<bool,2>,d> domain_walls;

    Array<Implicit_Object<T,d>*> sources;

    Smoke_Example();

    ~Smoke_Example()
    {if(hierarchy!=nullptr) delete hierarchy;}

//######################################################################
    virtual void Initialize_Rasterizer()=0;
    virtual void Initialize_Fluid_State()=0;
    virtual void Initialize_Sources()=0;
//######################################################################
    void Initialize();
    void Initialize_SPGrid();
    void Limit_Dt(T& dt,const T time) override;
    void Advect_Density(const T dt);
    void Diffuse_Density(const T dt);
    void Ficks_Diffusion(const T dt);
    void Non_Ficks_Diffusion(const T dt);
    void Reset_Solver_Channels();
    void Modify_Density_With_Sources();
    void Advect_Face_Velocities(const T dt);
    void Advect_Face_Qc(const T dt);
    void Set_Neumann_Faces_Inside_Sources();
    void Initialize_Values_At_Boundary_Conditions();
    void Set_Boundary_Conditions();
    void Project(const T dt);
    void Register_Options() override;
    void Parse_Options() override;
    void Read_Output_Files(const int frame);
    void Write_Output_Files(const int frame) const override;
//######################################################################
};
}
#endif

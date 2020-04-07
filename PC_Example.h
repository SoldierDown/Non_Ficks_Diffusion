//!#####################################################################
//! \file PC_Example.h
//!#####################################################################
// Class PC_Example
//######################################################################
#ifndef __PC_Example__
#define __PC_Example__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Tools/Utilities/Example.h>
#include "PC_Data.h"

namespace Nova{
template<class T,int d>
class PC_Example: public Example<T,d>
{
    using TV                        = Vector<T,d>;
    using Base                      = Example<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Struct_type               = PC_Data<T>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Rasterizer      = Hierarchical_Rasterizer<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    using Base::frame_title;using Base::output_directory;using Base::parse_args;using Base::first_frame;

    bool explicit_diffusion;
    // phase change paras
    // for updating density
    T tau_p;
    // for updating face_qs_channels
    T tau_1,SR,K;
    // for updating face_qt_channels
    T tau_2,k;
    // for updating m_channel
    T m_alpha,gamma,Teq;
    // for updating epsilon_channel
    int omega;
    T delta;
    T bv;
    T_INDEX counts;
    int levels,mg_levels,cg_iterations,cg_restart_iterations;
    T cfl,cg_tolerance;
    Hierarchy *hierarchy;
    Hierarchy_Rasterizer *rasterizer;

    T Struct_type::* density_channel;
    T Struct_type::* density_backup_channel;
    T Struct_type::* T_channel;
    T Struct_type::* T_backup_channel;
    T Struct_type::* epsilon_channel;
    
    Channel_Vector face_qs_channels;
    Channel_Vector face_qs_backup_channels;
    Channel_Vector face_qt_channels;
    Channel_Vector face_qt_backup_channels;

    Channel_Vector face_velocity_channels;

    Vector<Vector<bool,2>,d> domain_walls;

    Array<Implicit_Object<T,d>*> sources;

    PC_Example();
    ~PC_Example()
    {if(hierarchy!=nullptr) delete hierarchy;}

//######################################################################
    virtual void Initialize_Rasterizer()=0;
    virtual void Initialize_State()=0;
    virtual void Initialize_Sources()=0;
//######################################################################
    void Initialize();
    void Initialize_SPGrid();
    void Limit_Dt(T& dt,const T time) override;
    
    void Advect_Scalar_Field(const T dt);
    void Advect_Density(const T dt);
    void Advect_Temperature(const T dt);

    void Advect_Face_Vector_Field(const T dt);
    void Advect_Face_Qs(const T dt);
    void Advect_Face_Qt(const T dt);

    void Update_Density(const T dt);
    void Explicitly_Update_Density(const T dt);
    void Add_Divergence_Term_To_Density(const T dt);
    void Add_Differential_Term_To_Density(const T dt);
    void Add_Poly_Term_To_Density(const T dt);
    void Add_Random_Term_To_Density(const T dt);
    void Implicitly_Update_Density(const T dt);

    void Update_Face_Qs(const T dt);
    void Explicitly_Update_Face_Qs(const T dt);
    void Add_Gradient_Term_To_Face_Qs(const T dt);
    void Add_Linear_Term_To_Face_Qs(const T dt);
    void Implicitly_Update_Face_Qs(const T dt);


    void Update_Temperature(const T dt);
    void Explicitly_Update_Temperature(const T dt);
    void Add_Divergence_Term_To_Temperature(const T dt);
    void Add_Constant_Term_To_Temperature(const T dt);
    void Add_dSdt_Term_To_Temperature(const T dt);
    void Implicitly_Update_Temperature(const T dt);

    void Update_Face_Qt(const T dt);
    void Explicitly_Update_Face_Qt(const T dt);
    void Add_Gradient_Term_To_Face_Qt(const T dt);
    void Add_Linear_Term_To_Face_Qt(const T dt);
    void Implicitly_Update_Face_Qt(const T dt);

    void Backup();
    void Backup_Density();
    void Backup_Temperature();
    void Backup_Qs();
    void Backup_Qt();

    void Modify_Density_With_Sources();
    void Add_Source(const T dt);
    // void Advect_Face_Velocities(const T dt);
    void Set_Neumann_Faces_Inside_Sources();
    void Initialize_Velocity_Field();
    // void Project(const T dt);
    void Register_Options() override;
    void Parse_Options() override;
    void Read_Output_Files(const int frame);
    void Write_Output_Files(const int frame) const override;
    void Wrtie_To_File(const int frame);
//######################################################################
};
}
#endif
//!#####################################################################
//! \file MPM_CG_System.h
//!#####################################################################
// Class MPM_CG_System
//######################################################################
#ifndef __MPM_CG_System__
#define __MPM_CG_System__

#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include "MPM_CG_Vector.h"
#include "Clear_Non_Active_Helper.h"
#include "Convergence_Norm_Helper.h"
#include "Inner_Product_Helper.h"
#include "Multiply_Helper.h"
#include "MPM_Flags.h"

namespace Nova{
template<class Struct_type,class Multigrid_struct_type,class T,int d>
class MPM_CG_System: public Krylov_System_Base<T>
{
    using Base                      = Krylov_System_Base<T>;
    using Vector_Base               = Krylov_Vector_Base<T>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;

  public:
    Hierarchy& hierarchy;
    Channel_Vector channel;
    const int boundary_smoothing_iterations,interior_smoothing_iterations,bottom_smoothing_iterations;

    MPM_CG_System(Hierarchy& hierarchy_input,Channel_Vector& channel_input,
              const int boundary_smoothing_iterations_input,const int interior_smoothing_iterations_input,
              const int bottom_smoothing_iterations_input)
        :Base(true,false),hierarchy(hierarchy_input),channel(channel_input),
        boundary_smoothing_iterations(boundary_smoothing_iterations_input),
        interior_smoothing_iterations(interior_smoothing_iterations_input),bottom_smoothing_iterations(bottom_smoothing_iterations_input)
    {}

    ~MPM_CG_System() {}

    void Set_Boundary_Conditions(Vector_Base& v) const {}
    void Project_Nullspace(Vector_Base& v) const {}

    void Multiply(const Vector_Base& x,Vector_Base& result) const
    {
        Channel_Vector x_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(x).channel;
        Channel_Vector result_channel    = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(result).channel;
        
        for(int level=0;level<hierarchy.Levels();++level)
            Multiply_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),x_channel,result_channel,(unsigned)Node_Saturated);
        
    }

    void Project(Vector_Base& v) const
    {
        Channel_Vector v_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            Clear_Non_Active_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channel(v));
    }

    double Inner_Product(const Vector_Base& v1,const Vector_Base& v2) const
    {
        const Hierarchy& v1_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v1);
        const Hierarchy& v2_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v2);
        Channel_Vector const v1_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v1).channel;
        Channel_Vector const v2_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v2).channel;
        assert(&hierarchy == &v1_hierarchy);
        assert(&hierarchy == &v2_hierarchy);

        double result=0;

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            Inner_Product_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channel(v),
                                                  v2_channel(v),result,(unsigned)Node_Saturated);

        return result;
    }

    T Convergence_Norm(const Vector_Base& v) const
    {
        Channel_Vector v_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;
        T max_value=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level) 
            Convergence_Norm_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     v_channel,max_value,(unsigned)Node_Saturated);

        return max_value;
    }

    void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const
    {
        // T Struct_type::* r_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(r).channel;
        // T Struct_type::* z_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(z).channel;

        // multigrid_solver.Initialize_Right_Hand_Side(r_channel);
        // multigrid_solver.Initialize_Guess();
        // multigrid_solver.V_Cycle(boundary_smoothing_iterations,interior_smoothing_iterations,bottom_smoothing_iterations);

        // // clear z
        // for(int level=0;level<hierarchy.Levels();++level)
        //     SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),z_channel);

        // // copy u from multigrid hierarchy
        // multigrid_solver.Copy_Channel_Values(z_channel,multigrid_solver.u_channel,false);
    }
};
}
#endif

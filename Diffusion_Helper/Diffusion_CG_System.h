//!#####################################################################
//! \file Diffusion_CG_System.h
//!#####################################################################
// Class Diffusion_CG_System
//######################################################################
#ifndef __Diffusion_CG_System__
#define __Diffusion_CG_System__

#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>
#include "Diffusion_CG_Vector.h"
#include "Diffusion_Clear_Non_Active_Helper.h"
#include "Diffusion_Convergence_Norm_Helper.h"
#include "Diffusion_Inner_Product_Helper.h"
#include "Diffusion_Multiply_Helper.h"
#include "../MPM_Flags.h"

namespace Nova{
template<class Base_struct_type,class T,int d>
class Diffusion_CG_System: public Krylov_System_Base<T>
{
    using Base                      = Krylov_System_Base<T>;
    using Vector_Base               = Krylov_Vector_Base<T>;
    using Channel_Vector            = Vector<T Base_struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Base_struct_type,T,d>;

  public:
    const bool FICKS;
    T a;
    T four_a_plus_one;
    T coeff1;

    Hierarchy& hierarchy;

    Diffusion_CG_System(Hierarchy& hierarchy_input,bool FICKS_input)
        :Base(true,false),hierarchy(hierarchy_input),FICKS(FICKS_input),a((T)0.),four_a_plus_one((T)0.),coeff1((T)0.)
    {}

    ~Diffusion_CG_System() {}

    void Set_Boundary_Conditions(Vector_Base& v) const {}
    void Project_Nullspace(Vector_Base& v) const {}

    void Multiply(const Vector_Base& v,Vector_Base& result) const
    {
        T Base_struct_type::* v_channel         = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;
        T Base_struct_type::* result_channel    = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(result).channel;
        
        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Multiply_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channel,result_channel,
                                    FICKS,a,four_a_plus_one,coeff1);       
        
    }

    void Project(Vector_Base& v) const
    {
        T Base_struct_type::* v_channel         = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;

        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Clear_Non_Active_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channel);
    }

    double Inner_Product(const Vector_Base& v1,const Vector_Base& v2) const
    {
        const Hierarchy& v1_hierarchy           = Diffusion_CG_Vector<Base_struct_type,T,d>::Hierarchy(v1);
        const Hierarchy& v2_hierarchy           = Diffusion_CG_Vector<Base_struct_type,T,d>::Hierarchy(v2);
        T Base_struct_type::* const v1_channel  = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v1).channel;
        T Base_struct_type::* const v2_channel  = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v2).channel;
        assert(&hierarchy == &v1_hierarchy);
        assert(&hierarchy == &v2_hierarchy);

        double result=0;

        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Inner_Product_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channel,
                                                  v2_channel,result,(unsigned)Node_Saturated);

        return result;
    }

    T Convergence_Norm(const Vector_Base& v) const
    {
        T Base_struct_type::* v_channel         = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;
        T max_value=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Convergence_Norm_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     v_channel,max_value,(unsigned)Node_Saturated);

        return max_value;
    }

    void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const
    {
        // T Base_struct_type::* r_channel         = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(r).channel;
        // T Base_struct_type::* z_channel         = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(z).channel;

        // multigrid_solver.Initialize_Right_Hand_Side(r_channel);
        // multigrid_solver.Initialize_Guess();
        // multigrid_solver.V_Cycle(boundary_smoothing_iterations,interior_smoothing_iterations,bottom_smoothing_iterations);

        // // clear z
        // for(int level=0;level<hierarchy.Levels();++level)
        //     SPGrid::Clear<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),z_channel);

        // // copy u from multigrid hierarchy
        // multigrid_solver.Copy_Channel_Values(z_channel,multigrid_solver.u_channel,false);
    }
};
}
#endif

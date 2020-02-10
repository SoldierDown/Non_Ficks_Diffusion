//!#####################################################################
//! \file Diffusion_CG_System.h
//!#####################################################################
// Class Diffusion_CG_System
//######################################################################
#ifndef __Diffusion_CG_System__
#define __Diffusion_CG_System__

#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>
#include "Diffusion_CG_Vector.h"
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
    T twod_a_plus_one;
    T coeff1;

    Hierarchy& hierarchy;

    Diffusion_CG_System(Hierarchy& hierarchy_input,bool FICKS_input)
        :Base(true,false),hierarchy(hierarchy_input),FICKS(FICKS_input),a((T)0.),twod_a_plus_one((T)0.),coeff1((T)0.)
    {}

    ~Diffusion_CG_System(){}

    void Set_Boundary_Conditions(Vector_Base& v) const {}
    void Project_Nullspace(Vector_Base& v) const {}

    void Multiply(const Vector_Base& v,Vector_Base& result) const
    {
        T Base_struct_type::* v_channel         = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;
        T Base_struct_type::* result_channel    = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(result).channel;
        
        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Multiply_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channel,result_channel,
                                    FICKS,a,twod_a_plus_one,coeff1);       
        
    }

    void Project(Vector_Base& v) const
    {
    }

    double Inner_Product(const Vector_Base& v1,const Vector_Base& v2) const
    {
        const Hierarchy& v1_hierarchy           = Diffusion_CG_Vector<Base_struct_type,T,d>::Hierarchy(v1);
        const Hierarchy& v2_hierarchy           = Diffusion_CG_Vector<Base_struct_type,T,d>::Hierarchy(v2);
        T Base_struct_type::* const v1_channel  = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v1).channel;
        T Base_struct_type::* const v2_channel  = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v2).channel;
        assert(&hierarchy == &v1_hierarchy);    assert(&hierarchy == &v2_hierarchy);

        T result=(T)0;

        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Inner_Product_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channel,
                                                  v2_channel,result,(unsigned)Cell_Saturated);
        return result;
    }

    T Convergence_Norm(const Vector_Base& v) const
    {
        const Hierarchy& v_hierarchy           = Diffusion_CG_Vector<Base_struct_type,T,d>::Hierarchy(v);
        T Base_struct_type::* const v_channel  = Diffusion_CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;
        assert(&hierarchy == &v_hierarchy);
        T result=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Convergence_Norm_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                  v_channel,result,(unsigned)Cell_Saturated);
        return result;
    }

    void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const
    {
    }
};
}
#endif

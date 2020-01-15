//!#####################################################################
//! \file MPM_CG_Vector.h
//!#####################################################################
// Class MPM_CG_Vector
//######################################################################
#ifndef __MPM_CG_Vector__
#define __MPM_CG_Vector__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Arithmetic.h>
#include <nova/Tools/Krylov_Solvers/Krylov_Vector_Base.h>
#include <nova/Tools/Log/Debug_Utilities.h>
#include <cassert>
#include "MPM_Flags.h"

namespace Nova{
template<class Struct_type,class T,int d>
class MPM_CG_Vector: public Krylov_Vector_Base<T>
{
    using Base              = Krylov_Vector_Base<T>;
    using Hierarchy_type    = Grid_Hierarchy<Struct_type,T,d>;
    using Channel_Vector    = Vector<T Struct_type::*,d>;
    
  public:
    Hierarchy_type& hierarchy;
    Channel_Vector&  channel;

    MPM_CG_Vector(Hierarchy_type& hierarchy_input,Channel_Vector& channel_input)
        :hierarchy(hierarchy_input),channel(channel_input)
    {}

    static const MPM_CG_Vector& Cg_Vector(const Base& base)
    {return dynamic_cast<const MPM_CG_Vector&>(base);}

    static const Hierarchy_type& Hierarchy(const Base& base)
    {return dynamic_cast<const MPM_CG_Vector&>(base).hierarchy;}

    static Hierarchy_type& Hierarchy(Base& base)
    {return dynamic_cast<MPM_CG_Vector&>(base).hierarchy;}

    size_t Raw_Size() const
    {FUNCTION_IS_NOT_DEFINED();}

    T& Raw_Get(int i)
    {FUNCTION_IS_NOT_DEFINED();}


    Base& operator+=(const Base& bv)
    {
        const Hierarchy_type& bv_hierarchy      = MPM_CG_Vector::Hierarchy(bv);
        Channel_Vector const bv_channel         = Cg_Vector(bv).channel;
        assert(&hierarchy == &bv_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            SPGrid::Masked_Add<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),channel(v),
                                                bv_channel(v),channel(v),(unsigned)Node_Saturated);
    }

    Base& operator-=(const Base& bv)
    {
        const Hierarchy_type& bv_hierarchy      = MPM_CG_Vector::Hierarchy(bv);
        Channel_Vector const bv_channel         = Cg_Vector(bv).channel;
        assert(&hierarchy == &bv_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),channel(v),
                                                     bv_channel(v),channel(v),(unsigned)Node_Saturated);
    }

    Base& operator*=(const T a)
    {
        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            SPGrid::Masked_Scalar_Multiply<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                            channel(v),a,channel(v),(unsigned)Node_Saturated);
    }

    void Copy(const T c,const Base& bv)
    {
        const Hierarchy_type& bv_hierarchy      = MPM_CG_Vector::Hierarchy(bv);
        Channel_Vector const bv_channel         = Cg_Vector(bv).channel;
        assert(&hierarchy == &bv_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            SPGrid::Masked_Scalar_Multiply<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                            bv_channel(v),c,channel(v),(unsigned)Node_Saturated);
    }

    void Copy(const T c1,const Base& bv1,const Base& bv2)
    {
        const Hierarchy_type& bv1_hierarchy     = MPM_CG_Vector::Hierarchy(bv1);
        const Hierarchy_type& bv2_hierarchy     = MPM_CG_Vector::Hierarchy(bv2);
        Channel_Vector bv1_channel              = Cg_Vector(bv1).channel;
        Channel_Vector bv2_channel              = Cg_Vector(bv2).channel;
        assert(&hierarchy == &bv1_hierarchy);
        assert(&hierarchy == &bv2_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;v++)
            SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),c1,bv1_channel(v),
                                                  bv2_channel(v),channel(v),(unsigned)Node_Saturated);
    }
};
}
#endif

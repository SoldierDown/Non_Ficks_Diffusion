//!#####################################################################
//! \file CG_Vector.h
//!#####################################################################
// Class CG_Vector
//######################################################################
#ifndef __CG_Vector__
#define __CG_Vector__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Arithmetic.h>
#include <nova/Tools/Krylov_Solvers/Krylov_Vector_Base.h>
#include <nova/Tools/Log/Debug_Utilities.h>
#include <cassert>

namespace Nova{
template<class Struct_type,class T,int d>
class CG_Vector: public Krylov_Vector_Base<T>
{
    using Base              = Krylov_Vector_Base<T>;
    using Hierarchy_type    = Grid_Hierarchy<Struct_type,T,d>;

  public:
    Hierarchy_type& hierarchy;
    T Struct_type::* channel;

    CG_Vector(Hierarchy_type& hierarchy_input,T Struct_type::* channel_input)
        :hierarchy(hierarchy_input),channel(channel_input)
    {}

    static const CG_Vector& Cg_Vector(const Base& base)
    {return dynamic_cast<const CG_Vector&>(base);}

    static const Hierarchy_type& Hierarchy(const Base& base)
    {return dynamic_cast<const CG_Vector&>(base).hierarchy;}

    static Hierarchy_type& Hierarchy(Base& base)
    {return dynamic_cast<CG_Vector&>(base).hierarchy;}

    size_t Raw_Size() const
    {FUNCTION_IS_NOT_DEFINED();}

    T& Raw_Get(int i)
    {FUNCTION_IS_NOT_DEFINED();}

    Base& operator+=(const Base& bv)
    {
        const Hierarchy_type& bv_hierarchy      = CG_Vector::Hierarchy(bv);
        T Struct_type::* const bv_channel       = Cg_Vector(bv).channel;
        assert(&hierarchy == &bv_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level)
            SPGrid::Masked_Add<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),channel,
                                                bv_channel,channel,(unsigned)Cell_Type_Interior);
    }

    Base& operator-=(const Base& bv)
    {
        const Hierarchy_type& bv_hierarchy      = CG_Vector::Hierarchy(bv);
        T Struct_type::* const bv_channel       = Cg_Vector(bv).channel;
        assert(&hierarchy == &bv_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level)
            SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),channel,
                                                     bv_channel,channel,(unsigned)Cell_Type_Interior);
    }

    Base& operator*=(const T a)
    {
        for(int level=0;level<hierarchy.Levels();++level)
            SPGrid::Masked_Scalar_Multiply<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                            channel,a,channel,(unsigned)Cell_Type_Interior);
    }

    void Copy(const T c,const Base& bv)
    {
        const Hierarchy_type& bv_hierarchy      = CG_Vector::Hierarchy(bv);
        T Struct_type::* const bv_channel       = Cg_Vector(bv).channel;
        assert(&hierarchy == &bv_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level)
            SPGrid::Masked_Scalar_Multiply<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                            bv_channel,c,channel,(unsigned)Cell_Type_Interior);
    }

    void Copy(const T c1,const Base& bv1,const Base& bv2)
    {
        const Hierarchy_type& bv1_hierarchy     = CG_Vector::Hierarchy(bv1);
        const Hierarchy_type& bv2_hierarchy     = CG_Vector::Hierarchy(bv2);
        T Struct_type::* bv1_channel            = Cg_Vector(bv1).channel;
        T Struct_type::* bv2_channel            = Cg_Vector(bv2).channel;
        assert(&hierarchy == &bv1_hierarchy);
        assert(&hierarchy == &bv2_hierarchy);

        for(int level=0;level<hierarchy.Levels();++level)
            SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),c1,bv1_channel,
                                                  bv2_channel,channel,(unsigned)Cell_Type_Interior);
    }
};
}
#endif

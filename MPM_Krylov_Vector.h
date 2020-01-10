//!#####################################################################
//! \file MPM_Krylov_Vector.h
//!#####################################################################
// Class MPM_Krylov_Vector
//######################################################################
#ifndef __MPM_Krylov_Vector__
#define __MPM_Krylov_Vector__

#include <nova/Tools/Arrays/Array.h>
#include <nova/Tools/Krylov_Solvers/Krylov_Vector_Base.h>
#include <nova/Tools/Log/Debug_Utilities.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class T,int d>
class MPM_Krylov_Vector: public Krylov_Vector_Base<T>
{
    using TV                                = Vector<T,d>;
    using Base                              = Krylov_Vector_Base<T>;

  public:
    Array<TV>& array;

    MPM_Krylov_Vector(Array<TV>& array_input)
        :array(array_input)
    {}

    static Array<TV>& CG_Array(Base& base_array)
    {return dynamic_cast<MPM_Krylov_Vector<T,d>&>(base_array).array;}

    static const Array<TV>& CG_Array(const Base& base_array)
    {return dynamic_cast<const MPM_Krylov_Vector<T,d>&>(base_array).array;}

    Base& operator+=(const Base& bv)
    {array+=CG_Array(bv);return *this;}

    Base& operator-=(const Base& bv)
    {array-=CG_Array(bv);return *this;}

    Base& operator*=(const T a)
    {array*=a;return *this;}

    void Copy(const T c,const Base& bv)
    {Array<TV>::template Copy<T>(c,CG_Array(bv),array);}

    void Copy(const T c1,const Base& bv1,const Base& bv2)
    {Array<TV>::template Copy<T>(c1,CG_Array(bv1),CG_Array(bv2),array);}

    size_t Raw_Size() const
    {FUNCTION_IS_NOT_DEFINED();}

    T& Raw_Get(int i)
    {FUNCTION_IS_NOT_DEFINED();}
};
}
#endif
//!#####################################################################
//! \file Sphere_Implicit_Object.h
//!#####################################################################
// Class Sphere_Implicit_Object
//######################################################################
#ifndef __Sphere_Implicit_Object__
#define __Sphere_Implicit_Object__

#include <nova/Geometry/Basic_Geometry/Sphere.h>
#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Tools/Log/Debug_Utilities.h>

namespace Nova{
template<class T,int d>
class Sphere_Implicit_Object: public Implicit_Object<T,d>, public Sphere<T,d>
{
    using TV                        = Vector<T,d>;
    using Geometry_Base             = Sphere<T,d>;

  public:
    Sphere_Implicit_Object(const TV& center,const T& radius)
        :Geometry_Base(center,radius)
    {}

    ~Sphere_Implicit_Object() {}

//######################################################################
    T Signed_Distance(const TV& X) const override
    {return Geometry_Base::Signed_Distance(X);}
//######################################################################
};
}
#endif

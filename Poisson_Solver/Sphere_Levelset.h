//!#####################################################################
//! \file Sphere_Levelset.h
//!#####################################################################
// Class Sphere_Levelset
//######################################################################
#ifndef __Sphere_Levelset__
#define __Sphere_Levelset__

#include <nova/Geometry/Basic_Geometry/Sphere.h>

namespace Nova{
template<class T,int d>
class Sphere_Levelset: public Analytic_Levelset<T,d>
{
    using TV        = Vector<T,d>;

  public:
    Sphere<T,d> sphere;

    Sphere_Levelset(const TV& center,const T radius)
        :sphere(center,radius)
    {}

    ~Sphere_Levelset() {}

    T Signed_Distance(const TV& X) override
    {return sphere.Signed_Distance(X);}
};
}
#endif

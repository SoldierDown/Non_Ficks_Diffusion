//!#####################################################################
//! \file Sphere_Levelset.h
//!#####################################################################
// Class Sphere_Levelset
//######################################################################
#ifndef __Sphere_Levelset__
#define __Sphere_Levelset__

namespace Nova{
template<class T,int d>
class Sphere_Levelset
{
    using TV                            = Vector<T,d>;
    using T_INDEX                       = Vector<int,d>;
public:
    TV center;
    T radius;
    Sphere_Levelset(TV center_input,T radius_input):center(center_input),radius(radius_input)
    {}
    ~Sphere_Levelset(){}
    T Distance(const TV& location){
        return (location-center).Norm()-radius;
    }
    T Sign(const TV& location){
        return (location-center).Norm()>radius?(T)1.:(T)-1.;
    }
};



}
#endif

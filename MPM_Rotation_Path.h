//!#####################################################################
//! \file MPM_Rotation_Path.h
//!#####################################################################
// Class MPM_Rotation_Path
//######################################################################
#ifndef __MPM_Rotation_Path__
#define __MPM_Rotation_Path__

#include <nova/Tools/Vectors/Vector.h>

namespace Nova{

template<class T,int d>
struct MPM_Rotation_Path
{
    using TV = Vector<T,d>;
    typedef typename TV::SPIN T_SPIN;

    virtual ~MPM_Rotation_Path(){}
    virtual ROTATION<TV> R(const T t) const=0;
    virtual T_SPIN W(const T t) const=0;
};

template<class T,int d>
struct MPM_Rotation_Path_Uniform: public MPM_Rotation_Path<TV>
{
    using TV = Vector<T,d>;
    typedef typename TV::SPIN T_SPIN;

    const ROTATION<TV> r;
    const T_SPIN w;

    MPM_Rotation_Path_Uniform(const ROTATION<TV>& r_input,const T_SPIN& w_input):r(r_input),w(w_input){}
    virtual ~MPM_Rotation_Path_Uniform(){}

    virtual ROTATION<TV> R(const T t) const {return ROTATION<TV>::From_Rotation_Vector(w*t)*r;}
    virtual T_SPIN W(const T t) const {return w;}
};

template<class T,int d>
struct MPM_Rotation_Path_Static: public MPM_Rotation_Path<TV>
{
    using TV = Vector<T,d>;
    typedef typename TV::SPIN T_SPIN;

    const ROTATION<TV> r;

    MPM_Rotation_Path_Static(const ROTATION<TV>& r_input):r(r_input){}
    virtual ~MPM_Rotation_Path_STATIC(){}

    virtual ROTATION<TV> R(const T t) const {return r;}
    virtual T_SPIN W(const T t) const {return T_SPIN();}
};
}
#endif

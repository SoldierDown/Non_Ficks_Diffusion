//!#####################################################################
//! \file MPM_Translation_Path.h
//!#####################################################################
// Class MPM_Translation_Path
//######################################################################
#ifndef __MPM_Translation_Path__
#define __MPM_Translation_Path__

#include <nova/Tools/Vectors/Vector.h>

namespace Nova{

template<class T, int d>
struct MPM_Translation_Path
{
    using TV = Vector<T,d>;

    virtual ~MPM_Translation_Path(){}
    virtual TV X(const T t) const=0;
    virtual TV V(const T t) const=0;
};

template<class T,int d>
struct MPM_Translation_Path_Uniform: public MPM_Translation_Path<TV>
{
    using TV = Vector<T,d>;

    TV x;
    TV v;

    MPM_Translation_Path_Uniform(const TV& x_input,const TV& v_input):x(x_input),v(v_input){}
    virtual ~MPM_Translation_Path_Uniform(){}

    virtual TV X(const T t) const {return x+v*t;}
    virtual TV V(const T t) const {return v;}
};

template<class T,int d>
struct MPM_Translation_Path_Static: public MPM_Translation_Path<TV>
{
    using TV = Vector<T,d>;
    
    TV x;

    MPM_Translation_Path_Static(const TV& x_input):x(x_input){}
    virtual ~MPM_Translation_Path_Static(){}

    virtual TV X(const T t) const {return x;}
    virtual TV V(const T t) const {return TV();}
};
}
#endif

//!#####################################################################
//! \file MPM_Collidable_Object.h
//!#####################################################################
// Class MPM_Collidable_Object
//###################################################################### 
#ifndef __MPM_Collidable_Object__
#define __MPM_Collidable_Object__

#include <nova/Tools/Vectors/Vector.h>
#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include "MPM_Rotation_Path.h"
#include "MPM_Translation_Path.h"
namespace Nova{

template<class T,int d>
class MPM_COLLIDABLE_OBJECT
{
public:

    const T mu; // coefficient of friction

    Implicit_Object<T,d>* ls;
    MPM_Translation_Path<TV>* tp;
    MPM_Rotation_Path<TV>* rp;

    const bool sticky;

    MPM_Collidable_Object(const T mu_input,
        Implicit_Object<TV>* ls_input,
        MPM_Translation_Path<TV>* tp_input,
        MPM_Rotation_Path<TV>* rp_input,
        const bool sticky_input=false)
        :mu(mu_input),ls(ls_input),
        tp(tp_input),rp(rp_input),
        sticky(sticky_input)
    {}

    virtual ~MPM_Collidable_Object()
    {
        delete ls;
        delete tp;
        delete rp;
    }

    bool Collide(const TV& x,const T t,TV& v,TV& n) const
    {
        TV normal;
        if(Detect_Collision(x,t,&normal)){
            TV V=Velocity(x,t);
            v-=V;
            Collide_Static(x,t,normal,v);
            v+=V;
            n=normal;
            return true;}
        return false;
    }

    bool Detect_Collision(const TV& x,const T t,TV* normal=0) const
    {
        TV X=rp->R(t).Inverse_Rotate(x-tp->X(t));
        if(ls->Extended_Phi(X)<0){
            if(normal) *normal=rp->R(t).Rotate(ls->Normal(X));
            return true;}
        return false;
    }

    virtual TV Velocity(const TV& x,const T t) const
    {return TV::Cross_Product(rp->W(t),x-tp->X(t))+tp->V(t);}

protected:

    void Collide_Static(const TV& x,const T t,const TV& normal,TV& v) const
    {
        if(sticky){v=TV();return;}
        const T projection=TV::Dot_Product(normal,v);
        if(projection<0){
        // if(1){
            v-=normal*projection;
            if(-projection*mu<v.Magnitude())
                v+=v.Normalized()*projection*mu;
            else v=TV();}
    }
};

template<class TV>
class MPM_Collidable_Object_Static: public MPM_Collidable_Object<TV>
{
    typedef MPM_Collidable_Object<T, d> BASE;
    using TV = Vector<T,d>;    
public:

    MPM_Collidable_Object_Static(const T mu_input,IMPLICIT_OBJECT<TV>* const ls_input,const TV& tp_input,const ROTATION<TV>& rp_input,const bool sticky_input=false)
        :BASE(mu_input,ls_input,new MPM_Translation_Path_STATIC<TV>(tp_input),new MPM_Rotation_Path_STATIC<TV>(rp_input),sticky_input)
    {}

    virtual ~MPM_Collidable_Object_Static(){}  
};

template<class TV>
class MPM_Deformable_Collidable_Object: public MPM_Collidable_Object<TV>
{
    typedef MPM_Collidable_Object<TV> BASE;
    using TV = Vector<T,d>;    
    using MPM_Collidable_Object<TV>::ls;

public:

    MPM_Deformable_Collidable_Object(const T mu_input,MPM_MOVING_LEVELSET<TV>* const ls_input,const bool sticky_input=false)
        :BASE(mu_input,ls_input,new MPM_Translation_Path_Static<TV>(TV()),new MPM_Rotation_Path_Static<TV>(ROTATION<TV>()),sticky_input)
    {}

    virtual ~MPM_Deformable_Collidable_Object(){}

    virtual TV Velocity(const TV& x,const T t) const
    {
        MPM_MOVING_LEVELSET<TV>* mls=dynamic_cast<MPM_MOVING_LEVELSET<TV>*>(ls);
        return mls->Velocity(x);
    }
};
}
#endif
//!#####################################################################
//! \file Cropped_Influence_Iterator.h
//!#####################################################################
// Class Cropped_Influence_Iterator
//######################################################################
#ifndef __Cropped_Influence_Iterator__
#define __Cropped_Influence_Iterator__

#include <nova/Tools/Vectors/Vector.h>
#include <type_traits>
#include <ostream>
#include <iterator>
#include "Influence_Iterator.h"
#include "Interval.h"

namespace Nova{
template<class T,int d,class TV_INT=Vector<int,d> > 
class Cropped_Influence_Iterator: public Influence_Iterator<T,d,TV_INT>
{
    using Base              = Influence_Iterator<T,d,TV_INT>;
    using T_Particle        = MPM_Particle<T,d>;
    using Base::max_corner; using Base::min_corner; using Base::index; using Base::direction;
    
    int x_max;

  public:
    Cropped_Influence_Iterator(T_Particle& particle_input)
        :Base(particle_input)
    {index=min_corner;}
    
    Cropped_Influence_Iterator(const TV_INT& min_corner_input,const TV_INT& max_corner_input,const Interval<int>& x_interval_input,T_Particle& particle_input)
        :Base(min_corner_input,max_corner_input,particle_input)
    {   
        min_corner(0)=std::max(min_corner(0),x_interval_input.min_corner);
        x_max=std::min(max_corner(0),x_interval_input.max_corner);
        index=min_corner;
    }
    
    ~Cropped_Influence_Iterator() {}

    bool Valid() const
    {return index(0)<=x_max;}

};

}
#endif

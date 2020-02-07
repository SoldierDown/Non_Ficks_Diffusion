//!#####################################################################
//! \file Influence_Iterator.h
//!#####################################################################
// Class Influence_Iterator
//######################################################################
#ifndef __Influence_Iterator__
#define __Influence_Iterator__

#include <nova/Tools/Vectors/Vector.h>
#include <type_traits>
#include <ostream>
#include <iterator>
#include "../MPM_Particle.h"
namespace Nova{
template<class T,int d,class TV_INT=Vector<int,d> > 
class Influence_Iterator
{
    using T_Particle                = MPM_Particle<T,d>;
    using TV                        = Vector<T,d>;

public:
    TV_INT min_corner;
    TV_INT max_corner;
    TV_INT direction;
    TV_INT index;
    const T_Particle& particle;

  public:
    Influence_Iterator(const T_Particle& particle_input)
        :min_corner(TV_INT()),max_corner(TV_INT()),direction(TV_INT(1)),particle(particle_input)
    {Reset();}

    template< class RANGE_TYPE>
    Influence_Iterator(const RANGE_TYPE& range,const T_Particle& particle_input)
        :min_corner(range.min_corner),max_corner(range.max_corner),direction(initializeDirection(range.min_corner,range.max_corner)),particle(particle_input)
    {Reset();}
    
    Influence_Iterator(const TV_INT& min_corner_input,const TV_INT& max_corner_input,const T_Particle& particle_input)
        :min_corner(min_corner_input),max_corner(max_corner_input),direction(initializeDirection(min_corner,max_corner)),particle(particle_input)
    {Reset();}
    
    ~Influence_Iterator() {}
    
    static TV_INT initializeDirection(const TV_INT& min_corner,const TV_INT& max_corner)
    {
        TV_INT direction;
        for(int i=0;i<d;i++) direction[i]=min_corner[i]<max_corner[i]?1:-1;
        return direction;
    }

    Influence_Iterator operator++()
    {Next();if(!Valid()) SetEnd();return *this;}

    const TV_INT& operator*()
    {return Index();}

    bool operator!=(const Influence_Iterator& other)
    {return !(min_corner==other.min_corner && max_corner==other.max_corner && index==other.index);}

    void SetEnd()
    {index=max_corner;index[0]+=direction[0];}

    void Reset()
    {index=min_corner;}

    bool Valid() const
    {return (direction[0]*index[0])<=(direction[0]*max_corner[0]);}

    void Next()
    {index=Upcoming();}

    TV_INT Upcoming() const
    {
        TV_INT index_tmp = index;
        for(int i=d-1;i>=0;i--) if((direction[i]*index_tmp[i])<(direction[i]*max_corner[i]) || i==0){index_tmp[i]+=direction[i];return index_tmp;} else index_tmp[i]=min_corner[i];
        return index_tmp;
    }

    const TV_INT& Index() const
    {return index;}

    const TV_INT Current_Cell() const 
    {return particle.closest_cell+Index();}

    T Weight() const
    {
        T weight=(T)1.;
        for(int i=0;i<d;++i) weight*=particle.weights(index(i)+1,i);
        return weight;
    }

    TV Weight_Gradient() const
    {
        TV weight_gradient(TV()+1);
        for(int i=0;i<d;++i)
            for(int j=0;j<d;++j)
                weight_gradient(i)*=(i==j)?particle.dweights(index(j)+1,j):particle.weights(index(j)+1,j);
        return weight_gradient;
    }
};

template<class T,int d, class TV_INT>
std::ostream& operator<<(std::ostream& out,const Influence_Iterator<T,d,TV_INT>& iter)
{
    out<<"Index: "<<iter.Index()<<", Next: "<<iter.Upcoming()<<", Valid: "<<iter.Valid();
    return out;
}
}
#endif

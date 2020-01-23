//!#####################################################################
//! \file Interval.h
//!#####################################################################
// Class Interval
//######################################################################
#ifndef __Interval__
#define __Interval__

#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class T>
class Interval
{
public:
    T min_corner,max_corner;
    Interval(){}
    Interval(T min_corner_input,T max_corner_input): min_corner(min_corner_input),max_corner(max_corner_input){}

    ~Interval() {}

    bool Intersection(const Interval interval) const
    {return min_corner<=interval.max_corner && interval.min_corner<=max_corner;}
};
}
#endif

//!#####################################################################
//! \file Cropped_Range_Iterator.h
//!#####################################################################
// Class Cropped_Range_Iterator
//######################################################################
#ifndef __Cropped_Range_Iterator__
#define __Cropped_Range_Iterator__

#include <nova/Tools/Utilities/Range_Iterator.h>

#include <nova/Tools/Vectors/Vector.h>
#include <type_traits>
#include <ostream>
#include <iterator>

#include "Interval.h"

namespace Nova{
template<int d,class TV_INT=Vector<int,d> > 
class Cropped_Range_Iterator: public Range_Iterator<d,TV_INT>
{
    using Base              = Range_Iterator<d,TV_INT>;
    using Base::max_corner; using Base::min_corner; using Base::index; using Base::direction;

  public:
    Cropped_Range_Iterator()
        :Base()
    {Base::Reset();}
    
    Cropped_Range_Iterator(const TV_INT& min_corner_input,const TV_INT& max_corner_input,const Interval<int>& x_interval_input)
        :Base(min_corner_input,max_corner_input)
    {   
        min_corner(0)=std::max(min_corner(0),x_interval_input.min_corner);
        max_corner(0)=std::min(max_corner(0),x_interval_input.max_corner);
        direction=Base::initializeDirection(min_corner,max_corner);
        Base::Reset();}
    
    ~Cropped_Range_Iterator() {}

    bool Valid() const
    {return (direction[0]*index[0])<=(direction[0]*max_corner[0]);}

};

}
#endif

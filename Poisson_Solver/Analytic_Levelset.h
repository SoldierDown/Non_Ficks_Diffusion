//!#####################################################################
//! \file Analytic_Levelset.h
//!#####################################################################
// Class Analytic_Levelset
//######################################################################
#ifndef __Analytic_Levelset__
#define __Analytic_Levelset__

#include <nova/Tools/Vectors/Vector.h>
#include <limits>

namespace Nova{
template<class T,int d>
class Analytic_Levelset
{
    using TV            = Vector<T,d>;

  public:

    Analytic_Levelset() {}

    virtual ~Analytic_Levelset() {}

    virtual T Signed_Distance(const TV& X)
    {return std::numeric_limits<T>::max();}
};
}
#endif

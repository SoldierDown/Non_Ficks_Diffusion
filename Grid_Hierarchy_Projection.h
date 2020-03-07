//!#####################################################################
//! \file Grid_Hierarchy_Projection.h
//!#####################################################################
// Class Grid_Hierarchy_Projection
//######################################################################
#ifndef __Grid_Hierarchy_Projection__
#define __Grid_Hierarchy_Projection__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Projection/Ghost_Value_Propagate.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Grid_Hierarchy_Projection
{
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;

  public:
    Grid_Hierarchy_Projection() {}
    ~Grid_Hierarchy_Projection() {}

    static void Propagate_Ghost_Values(Hierarchy& hierarchy,T Struct_type::* v_channel)
    {
        for(int level=hierarchy.Levels()-2;level>=0;--level)
            Ghost_Value_Propagate<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),v_channel,v_channel,level);
    }
};
}
#endif

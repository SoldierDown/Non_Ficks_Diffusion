//!#####################################################################
//! \file Central_Rasterizer.h
//!#####################################################################
// Class Central_Rasterizer
//######################################################################
#ifndef __Central_Rasterizer__
#define __Central_Rasterizer__

#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Central_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using TV            = Vector<T,d>;
    using Base          = Hierarchical_Rasterizer<Struct_type,T,d>;
    using T_INDEX       = Vector<int,d>;
    using T_CELL        = std::pair<unsigned,T_INDEX>;
    using Hierarchy     = Grid_Hierarchy<Struct_type,T,d>;

    Range<T,d> domain;

  public:
    using Base::hierarchy;

    Central_Rasterizer(Hierarchy& hierarchy_input)
        :Base(hierarchy_input),domain(TV(-.125),TV(.125))
    {}

    bool Consume(const T_CELL& cell) override
    {
        const unsigned level=cell.first;
        const T_INDEX& index=cell.second;
        const TV dX=hierarchy.Lattice(level).dX;
        const TV X=hierarchy.Lattice(level).Node(index);

        if(level==0){hierarchy.Activate_Cell(level,index,Cell_Type_Interior);return false;}
        else if(domain.Inside(X) && domain.Inside(X+dX)){hierarchy.Activate_Cell(level,index,Cell_Type_Interior);return false;}

        return true;
    }
};
}
#endif

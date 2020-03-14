//!#####################################################################
//! \file Randomized_Rasterizer.h
//!#####################################################################
// Class Randomized_Rasterizer
//######################################################################
#ifndef __Randomized_Rasterizer__
#define __Randomized_Rasterizer__

#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Tools/Random_Numbers/Random_Numbers.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Randomized_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using Base          = Hierarchical_Rasterizer<Struct_type,T,d>;
    using T_INDEX       = Vector<int,d>;
    using T_CELL        = std::pair<unsigned,T_INDEX>;
    using Hierarchy     = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::hierarchy;

    Random_Numbers<T> random;

    Randomized_Rasterizer(Hierarchy& hierarchy_input)
        :Base(hierarchy_input)
    {random.Set_Seed(7);}

    bool Consume(const T_CELL& cell) override
    {
        const unsigned level=cell.first;
        const T_INDEX& index=cell.second;

        T value=random.Get_Uniform_Number((T)0.,(T)1.);

        if(level==0){hierarchy.Activate_Cell(level,index,Cell_Type_Interior);return false;}
        else if(value>(T).7){hierarchy.Activate_Cell(level,index,Cell_Type_Interior);return false;}

        return true;
    }
};
}
#endif

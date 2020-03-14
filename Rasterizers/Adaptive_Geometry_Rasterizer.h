//!#####################################################################
//! \file Adaptive_Geometry_Rasterizer.h
//!#####################################################################
// Class Adaptive_Geometry_Rasterizer
//######################################################################
#ifndef __Adaptive_Geometry_Rasterizer__
#define __Adaptive_Geometry_Rasterizer__

#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Rigid_Bodies/Rigid_Body.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Adaptive_Geometry_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using TV            = Vector<T,d>;
    using Base          = Hierarchical_Rasterizer<Struct_type,T,d>;
    using T_INDEX       = Vector<int,d>;
    using T_CELL        = std::pair<unsigned,T_INDEX>;
    using Hierarchy     = Grid_Hierarchy<Struct_type,T,d>;

    const Rigid_Body<T,d>& rigid_body;
    int shell_width;

  public:
    using Base::hierarchy;

    Adaptive_Geometry_Rasterizer(Hierarchy& hierarchy_input,const Rigid_Body<T,d>& rigid_body_input,const int shell_width_input=3)
        :Base(hierarchy_input),rigid_body(rigid_body_input),shell_width(shell_width_input)
    {}

    bool Consume(const T_CELL& cell) override
    {
        const unsigned level=cell.first;
        const T_INDEX& index=cell.second;
        const TV X=hierarchy.Lattice(level).Center(index);

        if(level==0){hierarchy.Activate_Cell(level,index,Cell_Type_Interior);
            return false;}
        else if(rigid_body.implicit_object->Signed_Distance(rigid_body.Object_Space_Point(X)) > shell_width*hierarchy.Lattice(level).dX.Max()){
            hierarchy.Activate_Cell(level,index,Cell_Type_Interior);return false;}

        return true;
    }
};
}
#endif

//!#####################################################################
//! \file Adaptive_Sphere_Rasterizer.h
//!#####################################################################
// Class Adaptive_Sphere_Rasterizer
//######################################################################
#ifndef __Adaptive_Sphere_Rasterizer__
#define __Adaptive_Sphere_Rasterizer__

#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Geometry/Basic_Geometry/Sphere.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Adaptive_Sphere_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using TV            = Vector<T,d>;
    using Base          = Hierarchical_Rasterizer<Struct_type,T,d>;
    using T_INDEX       = Vector<int,d>;
    using T_CELL        = std::pair<unsigned,T_INDEX>;
    using Hierarchy     = Grid_Hierarchy<Struct_type,T,d>;

    Sphere<T,d> sphere;
    int shell_width;

  public:
    using Base::hierarchy;

    Adaptive_Sphere_Rasterizer(Hierarchy& hierarchy_input,const TV& center,const T radius,const int shell_width_input=3)
        :Base(hierarchy_input),sphere(center,radius),shell_width(shell_width_input)
    {}

    bool Consume(const T_CELL& cell) override
    {
        const unsigned level=cell.first;
        const T_INDEX& index=cell.second;
        const TV X=hierarchy.Lattice(level).Center(index);

        if(level==0){
            if(sphere.Signed_Distance(X) >= sphere.radius) hierarchy.Activate_Cell(level,index,Cell_Type_Interior);
            return false;}
        else if(sphere.Signed_Distance(X) > sphere.radius+shell_width*hierarchy.Lattice(level).dX.Max()){
            hierarchy.Activate_Cell(level,index,Cell_Type_Interior);return false;}

        return true;
    }
};
}
#endif

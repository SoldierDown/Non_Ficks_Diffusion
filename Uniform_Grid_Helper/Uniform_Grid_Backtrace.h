//!#####################################################################
//! \file Uniform_Grid_Backtrace_Helper.h
//!#####################################################################
// Class Uniform_Grid_Backtrace_Helper
//######################################################################
#ifndef __Uniform_Grid_Backtrace_Helper__
#define __Uniform_Grid_Backtrace_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Tools/Log/Log.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Uniform_Grid_Backtrace_Helper
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Hierarchy_Lookup          = Grid_Hierarchy_Lookup<Struct_type,T,d>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    // Checked
    static TV Uniform_Grid_Cell_Backtrace(Hierarchy& hierarchy,int& level,const T_INDEX& cell_index,uint64_t& cell_offset,const TV& dX)
    {
        const int original_level=level;
        const uint64_t original_cell_offset=cell_offset;
        const Grid<T,d>& grid=hierarchy.Lattice(level);

        Vector<int,d> base_index(Flag_array_mask::LinearToCoord(cell_offset));

        // get raw backtraced location index
        TV d_i=dX*hierarchy.Lattice(level).one_over_dX;           // delta in cell coord = intra_cell_dX + dX/h (dX = -vel*dt)
        T_INDEX d_index;
        for(int axis=0;axis<d;++axis) d_index(axis)=std::floor(d_i(axis));
        TV weights=d_i-TV(d_index);
        T_INDEX backtraced_cell_index=cell_index+d_index;
        TV backtraced_cell_location=grid.Center(cell_index)+dX;
        // clamp index, weights to grid
        const TV cell_width=grid.dX;
        const TV min_corner=grid.domain.min_corner+TV(.5)*cell_width;
        const TV max_corner=grid.domain.max_corner-TV(.5)*cell_width;
        const Range<T,d> cell_domain(min_corner,max_corner);
        const Range<int,d> cell_indices=hierarchy.Lattice(level).Cell_Indices();        
        for(int axis=0;axis<d;++axis){
            if(backtraced_cell_location(axis)<=cell_domain.min_corner(axis)){
                backtraced_cell_index(axis)=cell_indices.min_corner(axis);
                weights(axis)=(T)0.;}
            if(backtraced_cell_location(axis)>=cell_domain.max_corner(axis)){
                backtraced_cell_index(axis)=cell_indices.max_corner(axis)-1;
                weights(axis)=(T)1.;}}
        cell_offset=Flag_array_mask::Linear_Offset(backtraced_cell_index._data);
        return weights;
    }

    static TV Uniform_Grid_Face_Backtrace(Hierarchy& hierarchy,int& level,const T_INDEX& face_index,uint64_t& face_offset,const TV& dX,const int& axis)
    {
        const int original_level=level;
        const uint64_t original_cell_offset=face_offset;
        const Grid<T,d>& grid=hierarchy.Lattice(level);

        // get raw backtraced location index
        TV d_i=dX*hierarchy.Lattice(level).one_over_dX;
        T_INDEX d_index;
        for(int axis=0;axis<d;++axis) d_index(axis)=std::floor(d_i(axis));
        TV weights=d_i-TV(d_index);
        T_INDEX backtraced_face_index=face_index+d_index;
        TV backtraced_face_location=grid.Face(axis,face_index)+dX;

        // clamp index, weights to grid
        const TV cell_width=grid.dX;
        const TV min_corner=grid.domain.min_corner+(TV(.5)-TV::Axis_Vector(axis)*(T).5)*cell_width;
        const TV max_corner=grid.domain.max_corner-(TV(.5)-TV::Axis_Vector(axis)*(T).5)*cell_width;
        const Range<T,d> face_domain(min_corner,max_corner);
        const Range<int,d> face_indices=Range<int,d>(hierarchy.Lattice(level).Cell_Indices().min_corner,
                                        hierarchy.Lattice(level).Cell_Indices().max_corner+Vector<int,d>::Axis_Vector(axis));    
        for(int axis=0;axis<d;++axis){
            if(backtraced_face_location(axis)<=face_domain.min_corner(axis)){
                backtraced_face_index(axis)=face_indices.min_corner(axis);
                weights(axis)=(T)0.;}
            if(backtraced_face_location(axis)>=face_domain.max_corner(axis)){
                backtraced_face_index(axis)=face_indices.max_corner(axis)-1;
                weights(axis)=(T)1.;}}
        face_offset=Flag_array_mask::Linear_Offset(backtraced_face_index._data);
        return weights;
    }
};
}
#endif

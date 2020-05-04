//!#####################################################################
//! \file Initialize_Neumann_Cells.h
//!#####################################################################
// Class Initialize_Neumann_Cells
//######################################################################
#ifndef __Initialize_Neumann_Cells__
#define __Initialize_Neumann_Cells__

#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Grids/Grid_Iterator_Cell.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Initialize_Neumann_Cells
{
    using T_INDEX                   = Vector<int,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Cell_Iterator             = Grid_Iterator_Cell<T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Initialize_Neumann_Cells(Hierarchy& hierarchy,const Vector<Vector<bool,2>,d>& domain_walls)
    {Run(hierarchy,domain_walls);}

    void Run(Hierarchy& hierarchy,const Vector<Vector<bool,2>,d>& domain_walls) const
    {
        uint64_t face_neighbor_offsets[number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        for(int level=0;level<hierarchy.Levels();++level){
            auto flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
            for(int axis=0;axis<d;++axis) for(int axis_side=1;axis_side<=2;++axis_side){int side=2*axis+axis_side;
                if(!domain_walls(axis)(axis_side-1)) for(Cell_Iterator iterator(hierarchy.Lattice(level),0,Grid<T,d>::Boundary_Interior_Region,side);iterator.Valid();iterator.Next()){
                    const T_INDEX& index=iterator.Cell_Index();uint64_t offset=Flag_array_mask::Linear_Offset(index._data);
                    if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(offset,Cell_Type_Interior))
                        flags(offset)&=~Cell_Type_Interior;}}}
    }
};
}
#endif

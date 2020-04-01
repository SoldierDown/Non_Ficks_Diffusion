//!#####################################################################
//! \file Write_To_File_Helper.h
//!#####################################################################
// Class Write_To_File_Helper
//######################################################################
#ifndef __Write_To_File_Helper__
#define __Write_To_File_Helper__

#include <iostream>
#include<stdlib.h>
#include <string>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

using namespace std;
namespace Nova{
template<class Struct_type,class T,int d>
class Write_To_File_Helper
{
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;
    using Block_Iterator        = SPGrid::SPGrid_Block_Iterator<Flag_array_mask>;

  public:
    Write_To_File_Helper(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel,const std::string& filename)
    {Run(hierarchy,allocator,blocks,density_channel,filename);}

    void Run(Grid_Hierarchy<Struct_type,T,2>& hierarchy,SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel,const std::string& filename) const
    {
        FILE* fp=fopen(filename.c_str(),"w");
        auto density=allocator.template Get_Const_Array<Struct_type,T>(density_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        const Grid<T,2>& grid=hierarchy.Lattice(0);
        auto write_to_file_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&(Cell_Type_Interior|Cell_Type_Dirichlet)){
                    Vector<int,2> index(Flag_array_mask::LinearToCoord(offset)); Vector<T,2> location=grid.Center(index);
                    fprintf(fp, "%.6f %.6f %.6f\n", location(0),location(1),density(offset));}
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            write_to_file_helper(offset);}
        fclose(fp);
    }

    void Run(Grid_Hierarchy<Struct_type,T,3>& hierarchy,SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* density_channel,const std::string& filename) const
    {
        FILE* fp=fopen(filename.c_str(),"w");
        auto density=allocator.template Get_Const_Array<Struct_type,T>(density_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        const Grid<T,3>& grid=hierarchy.Lattice(0);
        auto write_to_file_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&(Cell_Type_Interior|Cell_Type_Dirichlet)){
                    Vector<int,3> index(Flag_array_mask::LinearToCoord(offset)); Vector<T,3> location=grid.Center(index);
                    fprintf(fp, "%.6f %.6f %.6f %.6f\n", location(0),location(1),location(2),density(offset));}
        };
        for(Block_Iterator iterator(blocks);iterator.Valid();iterator.Next_Block()){
            uint64_t offset=iterator.Offset();
            write_to_file_helper(offset);}
        fclose(fp);
    }

};
}
#endif

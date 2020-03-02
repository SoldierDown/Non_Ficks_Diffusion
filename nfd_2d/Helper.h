#include <chrono>
// for test
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include "../Multigrid_Solver/Multigrid_Data.h"
#include "../Multigrid_Solver/Multigrid_Smoother.h"
#include "../Multigrid_Solver/Initial_Guess_Helper.h"
#include "../Multigrid_Solver/Sphere_Levelset.h"
#include "../Multigrid_Solver/Levelset_Initializer.h"
#include "../Multigrid_Solver/Neumann_BC_Initializer.h"
#include "../Multigrid_Solver/Mark_Boundary.h"
#include "../Multigrid_Solver/Multigrid_Solver.h"
#include <omp.h>
using namespace Nova;
using namespace std::chrono;

template<class Struct_type,class T,int d>
void Initialize_Guess(Grid_Hierarchy<Struct_type,T,d>& hierarchy,T Struct_type::* u_channel,const bool random_guess)
{
    using TV                                    = Vector<T,d>;
    using T_INDEX                               = Vector<int,d>;
    using Flags_type                            = typename Struct_type::Flags_type;
    using Allocator_type                        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
    
    Random_Numbers<T> random;
    random.Set_Seed(0);
    
    for(int level=0;level<hierarchy.Levels();++level)
        SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),u_channel);

    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(u_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        TV max_corner=hierarchy.Lattice(level).domain.max_corner; TV min_corner=hierarchy.Lattice(level).domain.min_corner;
        TV domain_range=max_corner-min_corner;
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            std::array<int,d> base_index_s=Flag_array_mask::LinearToCoord(offset);
            T_INDEX base_index=*reinterpret_cast<T_INDEX*>(&base_index_s);

            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                const T_INDEX index=base_index+range_iterator.Index();
                if(flags(offset)&Cell_Type_Interior){
                    if(random_guess){
                        data(offset)=random.Get_Uniform_Number(-1,1);}
                    else data(offset)=(T)0.;
                    if(false){
                        const TV X=hierarchy.Lattice(level).Center(index);
                        TV r_X=(X-min_corner)/domain_range;
                        data(offset)=(T)sin(two_pi*r_X(0))*sin(two_pi*r_X(1));}}
                range_iterator.Next();}}}
}

template<class Struct_type,class T,int d>
void Compute_Right_Hand_Side(Grid_Hierarchy<Struct_type,T,d>& hierarchy,T Struct_type::* x_channel,T Struct_type::* b_channel,
                                const bool FICKS,const T dt,const T diff_coeff,const T Fc,const T tau)
{
    using Flags_type                            = typename Struct_type::Flags_type;
    using Allocator_type                        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper                       = Grid_Topology_Helper<Flag_array_mask>;
    if(FICKS){
        for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
            auto x=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(x_channel);
            auto rhs=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(b_channel);
            auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
            const T one_over_dx2=hierarchy.Lattice(level).one_over_dX(0)*hierarchy.Lattice(level).one_over_dX(1);
            const T a=diff_coeff*dt*one_over_dx2; 
            uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
            Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
            for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
                for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                    if(flags(offset)&Cell_Type_Interior){
                        rhs(offset)=x(offset);
                        for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                            int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                            if(flags(neighbor_offset)&Cell_Type_Dirichlet) rhs(offset)+=a*x(neighbor_offset);}
                            // rhs(offset)=x(offset);
                            }}}}
    else{for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto x=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(x_channel);
        auto rhs=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(b_channel);
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        const T one_over_dx2=hierarchy.Lattice(level).one_over_dX(0)*hierarchy.Lattice(level).one_over_dX(1);
        const T coeff1=dt*diff_coeff*(Fc*tau+dt)*one_over_dx2/(dt+tau);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){rhs(offset)=x(offset);//-coeff2*div_Qc(offset);
                for(int face=0;face<Topology_Helper::number_of_faces_per_cell;++face){
                        int64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(flags(neighbor_offset)&Cell_Type_Dirichlet) rhs(offset)+=coeff1*x(neighbor_offset);}
                        
                        // rhs(offset)=x(offset);
                        }}}}    
}

template<class Struct_type,class T,int d>
void Initialize_Dirichlet_Cells(Grid_Hierarchy<Struct_type,T,d>& hierarchy,T Struct_type::* levelset_channel)
{
    using Flags_type                            = typename Struct_type::Flags_type;
    using Allocator_type                        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;

    for(int level=0;level<hierarchy.Levels();++level){auto blocks=hierarchy.Blocks(level);
        auto levelset=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(levelset_channel);
        auto flags=hierarchy.Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);

        for(int b=0;b<blocks.second;b++){uint64_t offset=blocks.first[b];
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior){
                    flags(offset)|=Cell_Type_Dirichlet;
                    flags(offset)&=~Cell_Type_Interior;}}}
}
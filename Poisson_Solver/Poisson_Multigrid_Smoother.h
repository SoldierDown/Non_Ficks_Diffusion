//!#####################################################################
//! \file Poisson_Multigrid_Smoother.h
//!#####################################################################
// Class Poisson_Multigrid_Smoother
//######################################################################
#ifndef __Poisson_Multigrid_Smoother__
#define __Poisson_Multigrid_Smoother__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/SPGrid/Tools/SPGrid_Arithmetic.h>
#include <nova/Tools/Utilities/Constants.h>
#include "Poisson_Multiply_Inverse_Diagonal.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Poisson_Multigrid_Smoother
{
    using Channel_Vector                        = Vector<T Struct_type::*,d>;
    using Hierarchy                             = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Projection                  = Grid_Hierarchy_Projection<Struct_type,T,d>;

  public:
    Poisson_Multigrid_Smoother() {}
    ~Poisson_Multigrid_Smoother() {}

    static void Multiply_With_System_Matrix(Hierarchy& hierarchy,T Struct_type::* u_channel,T Struct_type::* Lu_channel)
    {
        Hierarchy_Projection::Compute_Laplacian(hierarchy,u_channel,Lu_channel);
    }

    static void Compute_Residual(Hierarchy& hierarchy,Channel_Vector& gradient_channels,T Struct_type::* u_channel,
                                 T Struct_type::* b_channel,T Struct_type::* temp_channel,const unsigned mask)
    {
        const int levels=hierarchy.Levels();

        // clear temporary channel
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);

        // compute laplace
        Multiply_With_System_Matrix(hierarchy,u_channel,temp_channel);

        // subtract from right hand side
        for(int level=0;level<levels;++level)
            SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     b_channel,temp_channel,temp_channel,mask);
    }

    static void Compute_Residual(Hierarchy& hierarchy,Channel_Vector& gradient_channels,const std::pair<const uint64_t*,unsigned>& blocks,
                                 const int current_level,T Struct_type::* u_channel,T Struct_type::* b_channel,T Struct_type::* temp_channel,
                                 const unsigned mask)
    {
        const int levels=hierarchy.Levels();

        // clear temporary channel
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);

        // compute laplace
        Multiply_With_System_Matrix(hierarchy,u_channel,temp_channel);

        // subtract from right hand side
        SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy.Allocator(current_level),blocks,
                                                 b_channel,temp_channel,temp_channel,mask);
    }

    static void Jacobi_Iteration(Hierarchy& hierarchy,Channel_Vector& gradient_channels,const std::pair<const uint64_t*,unsigned>& blocks,
                                 const int level,T Struct_type::* u_channel,T Struct_type::* b_channel,T Struct_type::* temp_channel,
                                 const int iterations,const unsigned mask,const T omega=(T)two_thirds)
    {
        for(int i=0;i<iterations;++i){
            Compute_Residual(hierarchy,gradient_channels,blocks,level,u_channel,b_channel,temp_channel,mask);
            // residual <-- residual/diagonal
            Multiply_Inverse_Diagonal<Struct_type,T,d>(hierarchy,blocks,temp_channel,temp_channel,mask,level);
            // u <-- u + omega*(residual/diagonal)
            SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy.Allocator(level),blocks,omega,
                                                  temp_channel,u_channel,u_channel,mask);}
    }

    static void Exact_Solve(Hierarchy& hierarchy,Channel_Vector& gradient_channels,T Struct_type::* u_channel,T Struct_type::* b_channel,
                            T Struct_type::* temp_channel,const int iterations,const unsigned mask,const T omega=(T)two_thirds)
    {
        const int levels=hierarchy.Levels();

        for(int i=0;i<iterations;++i){
            Compute_Residual(hierarchy,gradient_channels,u_channel,b_channel,temp_channel,mask);
            // residual <-- residual/diagonal
            for(int level=0;level<levels;++level)
                Multiply_Inverse_Diagonal<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),temp_channel,temp_channel,mask,level);
            // u <-- u + omega*(residual/diagonal)
            for(int level=0;level<levels;++level)
                SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                      omega,temp_channel,u_channel,u_channel,mask);}
    }
};
}
#endif

//!#####################################################################
//! \file Multigrid_Smoother.h
//!#####################################################################
// Class Multigrid_Smoother
//######################################################################
#ifndef __Multigrid_Smoother__
#define __Multigrid_Smoother__
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/SPGrid/Tools/SPGrid_Arithmetic.h>
#include <nova/Tools/Utilities/Constants.h>
#include "Multiply_Inverse_Diagonal.h"
#include "../Diffusion_Helper/Diffusion_Multiply_Helper.h"

namespace Nova{
template<class Struct_type,class T,int d>
class Multigrid_Smoother
{
    using Channel_Vector                        = Vector<T Struct_type::*,d>;
    using Hierarchy                             = Grid_Hierarchy<Struct_type,T,d>;
    // using Hierarchy_Projection                  = Grid_Hierarchy_Projection<Struct_type,T,d>;

  public:
    Multigrid_Smoother() {}
    ~Multigrid_Smoother() {}

    static void Multiply_With_System_Matrix(Hierarchy& hierarchy,Channel_Vector& gradient_channels,T Struct_type::* saturation_channel,T Struct_type::* result_channel,
                                                bool FICKS,const T a,const T twod_a_plus_one,const T coeff1)
    {
        for(int level=0;level<hierarchy.Levels();++level)
            Diffusion_Multiply_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),saturation_channel,result_channel,
                                    FICKS,a,twod_a_plus_one,coeff1);  
        // Hierarchy_Projection::Compute_Laplacian(hierarchy,gradient_channels,u_channel,Lu_channel);
    }

    static void Compute_Residual(Hierarchy& hierarchy,Channel_Vector& gradient_channels,T Struct_type::* saturation_channel,
                                 T Struct_type::* rhs_channel,T Struct_type::* result_channel,const unsigned mask,
                                 bool FICKS,const T a,const T twod_a_plus_one,const T coeff1)
    {
        const int levels=hierarchy.Levels();

        // clear temporary channel
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),result_channel);

        // compute laplace
        Multiply_With_System_Matrix(hierarchy,gradient_channels,saturation_channel,result_channel,FICKS,a,twod_a_plus_one,coeff1);

        // subtract from right hand side
        for(int level=0;level<levels;++level)
            SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     rhs_channel,result_channel,result_channel,mask);
    }

    static void Compute_Residual(Hierarchy& hierarchy,Channel_Vector& gradient_channels,const std::pair<const uint64_t*,unsigned>& blocks,
                                 const int current_level,T Struct_type::* saturation_channel,T Struct_type::* rhs_channel,T Struct_type::* result_channel,const unsigned mask,
                                 bool FICKS,const T a,const T twod_a_plus_one,const T coeff1)
    {
        const int levels=hierarchy.Levels();
        // clear temporary channel
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),result_channel);
        // compute laplace
        Multiply_With_System_Matrix(hierarchy,gradient_channels,saturation_channel,result_channel,FICKS,a,twod_a_plus_one,coeff1);

        // subtract from right hand side
        SPGrid::Masked_Subtract<Struct_type,T,d>(hierarchy.Allocator(current_level),blocks,
                                                 rhs_channel,result_channel,result_channel,mask);
    }

    static void Jacobi_Iteration(Hierarchy& hierarchy,Channel_Vector& gradient_channels,const std::pair<const uint64_t*,unsigned>& blocks,
                                 const int level,T Struct_type::* saturation_channel,T Struct_type::* rhs_channel,T Struct_type::* result_channel,
                                 const int iterations,const unsigned mask,
                                 bool FICKS,const T a,const T twod_a_plus_one,const T coeff1,
                                 const T omega=(T)two_thirds)
    {
        for(int i=0;i<iterations;++i){
            Compute_Residual(hierarchy,gradient_channels,blocks,level,saturation_channel,rhs_channel,result_channel,mask,FICKS,a,twod_a_plus_one,coeff1);
            // residual <-- residual/diagonal
            Multiply_Inverse_Diagonal<Struct_type,T,d>(hierarchy,blocks,result_channel,result_channel,mask,level,FICKS,twod_a_plus_one,coeff1);
            // u <-- u + omega*(residual/diagonal)
            SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy.Allocator(level),blocks,omega,
                                                  result_channel,saturation_channel,saturation_channel,mask);}
    }

    static void Exact_Solve(Hierarchy& hierarchy,Channel_Vector& gradient_channels,T Struct_type::* saturation_channel,T Struct_type::* rhs_channel,
                            T Struct_type::* result_channel,const int iterations,const unsigned mask,
                            bool FICKS,const T a,const T twod_a_plus_one,const T coeff1,
                            const T omega=(T)two_thirds)
    {
        const int levels=hierarchy.Levels();

        for(int i=0;i<iterations;++i){
            Compute_Residual(hierarchy,gradient_channels,saturation_channel,rhs_channel,result_channel,mask,FICKS,a,twod_a_plus_one,coeff1);
            // residual <-- residual/diagonal
            for(int level=0;level<levels;++level)
                Multiply_Inverse_Diagonal<Struct_type,T,d>(hierarchy,hierarchy.Blocks(level),result_channel,result_channel,mask,level,FICKS,twod_a_plus_one,coeff1);
            // u <-- u + omega*(residual/diagonal)
            for(int level=0;level<levels;++level)
                SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                      omega,result_channel,saturation_channel,saturation_channel,mask);}
    }
};
}
#endif

//!#####################################################################
//! \file MPM_CG_System.h
//!#####################################################################
// Class MPM_CG_System
//######################################################################
#ifndef __MPM_CG_System__
#define __MPM_CG_System__

#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Matrices/Matrix_2x2.h>
#include <nova/Tools/Matrices/Matrix_3x3.h>
#include <nova/Tools/Matrices/Diagonal_Matrix.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_2x2.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_3x3.h>

#include "MPM_CG_Vector.h"
#include "Clear_Non_Active_Helper.h"
#include "Convergence_Norm_Helper.h"
#include "Inner_Product_Helper.h"
#include "Multiply_Helper.h"
#include "../MPM_Example.h"
#include "../MPM_Particle.h"
#include "../MPM_Flags.h"

namespace Nova{
template<class Struct_type,class T,int d>
class MPM_CG_System: public Krylov_System_Base<T>
{
    using Base                      = Krylov_System_Base<T>;
    using Vector_Base               = Krylov_Vector_Base<T>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using T_Particle                = MPM_Particle<T,d>;
    using T_Range_Iterator          = Range_Iterator<d,T_INDEX>;

  public:
    Hierarchy& hierarchy;
    Array<T_Particle>& particles;
    Array<int> simulated_particles;
    T Struct_type::* mass_channel;
    const T trapezoidal;
    const T dt;

    MPM_CG_System(Hierarchy& hierarchy_input,Array<int>& simulated_particles_input,Array<T_Particle>& particles_input,T Struct_type::* mass_channel_input,const T trapezoidal_input,const T dt_input)
        :Base(true,false),hierarchy(hierarchy_input),simulated_particles(simulated_particles_input),particles(particles_input),mass_channel(mass_channel_input),trapezoidal(trapezoidal_input),dt(dt_input)
    {}
    ~MPM_CG_System() {}

    void Set_Boundary_Conditions(Vector_Base& v) const {}
    void Project_Nullspace(Vector_Base& v) const {}

    void Multiply(const Vector_Base& x,Vector_Base& result) const
    {
        Channel_Vector& x_channels           = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(x).channel;
        Channel_Vector& result_channels      = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(result).channel;
        Force(result_channels,x_channels);
        const T scaled_dt_squared=dt*dt/((T)2.+trapezoidal);
        for(int level=0;level<hierarchy.Levels();++level)
            Multiply_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),x_channels,result_channels,mass_channel,scaled_dt_squared,(unsigned)Node_Saturated);
    }
    
    void Force(Channel_Vector f,const Channel_Vector x) const
    {
        const Grid<T,d>& grid=hierarchy.Lattice(0);
        const TV one_over_dX=grid.one_over_dX;
    #pragma omp parallel for
        for(int i=0;i<simulated_particles.size();++i){
            int id=simulated_particles(i); T_Particle& p=particles(id);
            Matrix<T,d> tmp_mat=Matrix<T,d>();
            T_INDEX closest_node=grid.Closest_Node(p.X); 
            for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){
                T_INDEX current_node=closest_node+iterator.Index();
                if(grid.Node_Indices().Inside(current_node)){
                    const TV current_node_location=grid.Node(current_node); TV weight_grad=dN2<T,d>(p.X-current_node_location,one_over_dX);
                    tmp_mat+=Matrix<T,d>::Outer_Product(current_node_location,weight_grad);}}
            Matrix<T,d> F=p.constitutive_model.Fe; const T kp=(T)1e4; const T saturation=p.saturation; const T eta=p.constitutive_model.eta;
            tmp_mat=F.Times_Transpose(p.constitutive_model.Times_dP_dF(tmp_mat.Transpose_Times(F))
                                        -eta*kp*saturation*Times_Cofactor_Matrix_Derivative(F,tmp_mat.Transpose_Times(F)));
            p.scp=p.volume*tmp_mat;}
        
        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),f(v));  
        
        for(int i=0;i<simulated_particles.size();++i){
            int id=simulated_particles(i); T_Particle& p=particles(id);
            T_INDEX closest_node=grid.Closest_Node(p.X); 
            for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){
                T_INDEX current_node=closest_node+iterator.Index();
                if(grid.Node_Indices().Inside(current_node)){
                    const TV current_node_location=grid.Node(current_node); TV weight_grad=dN2<T,d>(p.X-current_node_location,one_over_dX);
                    TV tmp_vec=p.scp.Transpose_Times(weight_grad); 
                    for(int v=0;v<d;++v) hierarchy.Channel(0,f(v))(current_node._data)+=tmp_vec(v);}}}
    }

    void Project(Vector_Base& v) const
    {
        Channel_Vector v_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            Clear_Non_Active_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channel(v));
    }

    double Inner_Product(const Vector_Base& v1,const Vector_Base& v2) const
    {
        const Hierarchy& v1_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v1);
        const Hierarchy& v2_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v2);
        Channel_Vector const v1_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v1).channel;
        Channel_Vector const v2_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v2).channel;
        assert(&hierarchy == &v1_hierarchy);
        assert(&hierarchy == &v2_hierarchy);

        double result=0;

        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v)
            Inner_Product_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channel(v),
                                                  v2_channel(v),result,(unsigned)Node_Saturated);

        return result;
    }

    T Convergence_Norm(const Vector_Base& v) const
    {
        Channel_Vector v_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;
        T max_value=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level) 
            Convergence_Norm_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     v_channel,max_value,(unsigned)Node_Saturated);

        return max_value;
    }

    void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const
    {
        // T Struct_type::* r_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(r).channel;
        // T Struct_type::* z_channel         = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(z).channel;

        // multigrid_solver.Initialize_Right_Hand_Side(r_channel);
        // multigrid_solver.Initialize_Guess();
        // multigrid_solver.V_Cycle(boundary_smoothing_iterations,interior_smoothing_iterations,bottom_smoothing_iterations);

        // // clear z
        // for(int level=0;level<hierarchy.Levels();++level)
        //     SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),z_channel);

        // // copy u from multigrid hierarchy
        // multigrid_solver.Copy_Channel_Values(z_channel,multigrid_solver.u_channel,false);
    }
};
}
#endif

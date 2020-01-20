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
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Tools/Matrices/Matrix_2x2.h>
#include <nova/Tools/Matrices/Matrix_3x3.h>
#include <nova/Tools/Matrices/Diagonal_Matrix.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_2x2.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_3x3.h>

#include <chrono>

#include "MPM_CG_Vector.h"
#include "Clear_Non_Active_Helper.h"
#include "Convergence_Norm_Helper.h"
#include "Inner_Product_Helper.h"
#include "Multiply_Helper.h"
#include "../MPM_Plane_Barrier.h"
#include "../MPM_Example.h"
#include "../MPM_Particle.h"
#include "../MPM_Flags.h"

using namespace std::chrono;

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
    using T_Barrier                 = MPM_Plane_Barrier<T,d>;



  public:
    Hierarchy& hierarchy;
    Array<T_Particle>& particles;
    Array<int>& simulated_particles;
    Array<T_Barrier>& barriers;
    const T trapezoidal;
    const T dt;
    bool print_running_time;


    MPM_CG_System(Hierarchy& hierarchy_input,Array<int>& simulated_particles_input,Array<T_Particle>& particles_input,Array<T_Barrier>& barriers_input,const T trapezoidal_input,const T dt_input)
        :Base(true,false),hierarchy(hierarchy_input),simulated_particles(simulated_particles_input),particles(particles_input),barriers(barriers_input),trapezoidal(trapezoidal_input),dt(dt_input),
        print_running_time(false)
    {}
    ~MPM_CG_System() {}

    void Set_Boundary_Conditions(Vector_Base& v) const {}
    void Project_Nullspace(Vector_Base& v) const {}

    void Multiply(const Vector_Base& x,Vector_Base& result) const
    {
        
        Channel_Vector& x_channels           = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(x).channel;
        Channel_Vector& result_channels      = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(result).channel;
        Force(result_channels,x_channels);
        high_resolution_clock::time_point tb=high_resolution_clock::now();
        const T scaled_dt_squared=dt*dt/((T)1.+trapezoidal);
        for(int level=0;level<hierarchy.Levels();++level)
            Multiply_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),x_channels,result_channels,scaled_dt_squared,(unsigned)Node_Saturated);
        if(print_running_time){ high_resolution_clock::time_point te=high_resolution_clock::now();
        duration<double> dur=duration_cast<duration<double> >(te-tb);
        Log::cout<< "Multiply duration: "<<dur.count()<<std::endl;}

    }
    
    void Force(Channel_Vector& f,const Channel_Vector& x) const
    {
        high_resolution_clock::time_point tb1=high_resolution_clock::now();
        unsigned Struct_type::* flags=&Struct_type::flags;
        const Grid<T,d>& grid=hierarchy.Lattice(0);
        const TV one_over_dX=grid.one_over_dX;
#pragma omp parallel for
        for(int i=0;i<simulated_particles.size();++i){
            int id=simulated_particles(i); T_Particle& p=particles(id);
            Matrix<T,d> tmp_mat;
            T_INDEX closest_node=grid.Closest_Node(p.X); 
            for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){
                T_INDEX current_node=closest_node+iterator.Index();
                if(grid.Node_Indices().Inside(current_node)){       
                    const TV current_node_location=grid.Node(current_node);
                    T weight=N2<T,d>(p.X-current_node_location,one_over_dX);                                
                    if(weight>(T)0.){
                        TV weight_grad=dN2<T,d>(p.X-current_node_location,one_over_dX),v_vec;
                        for(int v=0;v<d;++v) v_vec(v)=hierarchy.Channel(0,x(v))(current_node._data);
                        tmp_mat+=Matrix<T,d>::Outer_Product(weight_grad,v_vec);}}}
            Matrix<T,d> F=p.constitutive_model.Fe; const T kp=(T)1e4; const T saturation=p.saturation; const T eta=p.constitutive_model.eta*saturation;
            tmp_mat=F.Times_Transpose(p.constitutive_model.Times_dP_dF(tmp_mat.Transpose_Times(F))-eta*kp*saturation*Times_Cofactor_Matrix_Derivative(F,tmp_mat.Transpose_Times(F)));
            p.scp=p.volume*tmp_mat;}
        if(print_running_time){high_resolution_clock::time_point te1=high_resolution_clock::now();
        duration<double> dur1=duration_cast<duration<double> >(te1-tb1);
        Log::cout<< "Collect to particles duration: "<<dur1.count()<<std::endl;}


        high_resolution_clock::time_point tb2=high_resolution_clock::now();
        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),f(v));  
        if(print_running_time){high_resolution_clock::time_point te2=high_resolution_clock::now();
        duration<double> dur2=duration_cast<duration<double> >(te2-tb2);
        Log::cout<< "Clear duration: "<<dur2.count()<<std::endl;}
      

        high_resolution_clock::time_point tb3=high_resolution_clock::now();
        for(int i=0;i<simulated_particles.size();++i){
            int id=simulated_particles(i); T_Particle& p=particles(id);
            T_INDEX closest_node=grid.Closest_Node(p.X); 
            for(T_Range_Iterator iterator(T_INDEX(-2),T_INDEX(2));iterator.Valid();iterator.Next()){
                T_INDEX current_node=closest_node+iterator.Index();
                if(grid.Node_Indices().Inside(current_node)){
                    const TV current_node_location=grid.Node(current_node);
                    T weight=N2<T,d>(p.X-current_node_location,one_over_dX);
                    if(weight>(T)0.){
                        TV weight_grad=dN2<T,d>(p.X-current_node_location,one_over_dX);
                        TV tmp_vec=p.scp.Transpose_Times(weight_grad); 
                        for(int v=0;v<d;++v) hierarchy.Channel(0,f(v))(current_node._data)+=tmp_vec(v);}}}}
        if(print_running_time){high_resolution_clock::time_point te3=high_resolution_clock::now();
        duration<double> dur3=duration_cast<duration<double> >(te3-tb3);
        Log::cout<< "Collect to grid duration: "<<dur3.count()<<std::endl;}

    }

    void Project(Vector_Base& v) const
    {
        high_resolution_clock::time_point tb=high_resolution_clock::now();
        Channel_Vector& v_channels              = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;

        for(int level=0;level<hierarchy.Levels();++level)
            Clear_Non_Active_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channels,barriers);
        if(print_running_time){ high_resolution_clock::time_point te=high_resolution_clock::now();
        duration<double> dur=duration_cast<duration<double> >(te-tb);
        Log::cout<< "Project duration: "<<dur.count()<<std::endl;}

    }

    double Inner_Product(const Vector_Base& v1,const Vector_Base& v2) const
    {
        high_resolution_clock::time_point tb=high_resolution_clock::now();
        const Hierarchy& v1_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v1);
        const Hierarchy& v2_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v2);
        Channel_Vector const v1_channels        = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v1).channel;
        Channel_Vector const v2_channels        = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v2).channel;
        assert(&hierarchy == &v1_hierarchy);
        assert(&hierarchy == &v2_hierarchy);

        double result=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level)
            Inner_Product_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channels,v2_channels,result,(unsigned)Node_Saturated);
        if(print_running_time){ high_resolution_clock::time_point te=high_resolution_clock::now();
        duration<double> dur=duration_cast<duration<double> >(te-tb);
        Log::cout<< "Inner Product duration: "<<dur.count()<<std::endl;}

        return result;
    }

    T Convergence_Norm(const Vector_Base& v) const
    {
        high_resolution_clock::time_point tb=high_resolution_clock::now();
        Channel_Vector& v_channels              = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;
        T result=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level)
            Convergence_Norm_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     v_channels,result,(unsigned)Node_Saturated);
        if(print_running_time){ high_resolution_clock::time_point te=high_resolution_clock::now();
        duration<double> dur=duration_cast<duration<double> >(te-tb);
        Log::cout<< "Convergence Norm duration: "<<dur.count()<<std::endl;}
        return result;
    }

    void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const
    {
    }
};
}
#endif

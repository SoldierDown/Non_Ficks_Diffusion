//!#####################################################################
//! \file MPM_CG_System.h
//!#####################################################################
// Class MPM_CG_System
//######################################################################
#ifndef __MPM_CG_System__
#define __MPM_CG_System__

#include <chrono>

#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/Tools/Matrices/Matrix_2x2.h>
#include <nova/Tools/Matrices/Matrix_3x3.h>
#include <nova/Tools/Matrices/Diagonal_Matrix.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_2x2.h>
#include <nova/Tools/Matrices/Symmetric_Matrix_3x3.h>

#include "MPM_CG_Vector.h"
#include "Project_Helper.h"
#include "Convergence_Norm_Helper.h"
#include "Inner_Product_Helper.h"
#include "Multiply_Helper.h"
#include "../MPM_Plane_Barrier.h"
#include "../MPM_Example.h"
#include "../MPM_Particle.h"
#include "../MPM_Flags.h"
#include "../Tools/Influence_Iterator.h"
#include "../Tools/Cropped_Influence_Iterator.h"
#include "../Tools/Interval.h"
#include "../Tools/Matrix_MXN.h"
#include "../Channel_Vector_Traverse_Helper.h"

using namespace std::chrono;

namespace Nova{
template<class Struct_type,class T,int d>
class MPM_CG_System: public Krylov_System_Base<T>
{
    using Base                          = Krylov_System_Base<T>;
    using Vector_Base                   = Krylov_Vector_Base<T>;
    using Channel_Vector                = Vector<T Struct_type::*,d>;
    using Hierarchy                     = Grid_Hierarchy<Struct_type,T,d>;
    using TV                            = Vector<T,d>;
    using T_INDEX                       = Vector<int,d>;
    using T_Particle                    = MPM_Particle<T,d>;
    using T_Influence_Iterator          = Influence_Iterator<T,d,T_INDEX>;
    using T_Cropped_Influence_Iterator  = Cropped_Influence_Iterator<T,d,T_INDEX>;
    using T_Barrier                     = MPM_Plane_Barrier<T,d>;

  public:
    Hierarchy& hierarchy;
    Array<T_Particle>& particles;
    Array<int>& simulated_particles;
    Array<T_Barrier>& barriers;
    Array<Interval<int> >& x_intervals;
    Matrix_MxN<Array<int> >& particle_bins;
    
    const T trapezoidal;
    const T dt;
    const int threads;


    MPM_CG_System(Hierarchy& hierarchy_input,Array<int>& simulated_particles_input,Array<T_Particle>& particles_input,Matrix_MxN<Array<int> >& particle_bins_input,Array<Interval<int> >& x_intervals_input,Array<T_Barrier>& barriers_input,const T trapezoidal_input,const T dt_input,const int threads_input)
        :Base(true,false),hierarchy(hierarchy_input),simulated_particles(simulated_particles_input),particles(particles_input),particle_bins(particle_bins_input),x_intervals(x_intervals_input),barriers(barriers_input),trapezoidal(trapezoidal_input),dt(dt_input),threads(threads_input)
    {}
    ~MPM_CG_System() {}

    void Set_Boundary_Conditions(Vector_Base& v) const {}
    void Project_Nullspace(Vector_Base& v) const {}

    void Multiply(const Vector_Base& x,Vector_Base& result) const
    {
        Channel_Vector& x_channels           = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(x).channel;
        Channel_Vector& result_channels      = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(result).channel;
        Force(result_channels,x_channels);
        // Log::cout<<"After Force() x: "<<Convergence_Norm(x)<<", result: "<<Convergence_Norm(result)<<std::endl;
        const T scaled_dt_squared=dt*dt/(1+trapezoidal);
        for(int level=0;level<hierarchy.Levels();++level)
            Multiply_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),x_channels,result_channels,scaled_dt_squared,(unsigned)Node_Saturated);
    }
    
    void Force(Vector<T Struct_type::*,2>& f,const Vector<T Struct_type::*,2>& x) const
    {
        bool print_running_time=false;
        using TV                            = Vector<T,2>;
        using T_INDEX                       = Vector<int,2>;
        using T_Particle                    = MPM_Particle<T,2>;
        using T_Influence_Iterator          = Influence_Iterator<T,2,T_INDEX>;
        using T_Cropped_Influence_Iterator  = Cropped_Influence_Iterator<T,2,T_INDEX>;
        using T_Barrier                     = MPM_Plane_Barrier<T,2>; 
        
        high_resolution_clock::time_point tb1 = high_resolution_clock::now();
        auto v0=hierarchy.Channel(0,x(0)); auto v1=hierarchy.Channel(0,x(1));
        auto f0=hierarchy.Channel(0,f(0)); auto f1=hierarchy.Channel(0,f(1));
        const Grid<T,2>& grid=hierarchy.Lattice(0);
#pragma omp parallel for
        for(int i=0;i<simulated_particles.size();++i){
            high_resolution_clock::time_point tb = high_resolution_clock::now();    
            int id=simulated_particles(i); T_Particle& p=particles(id); Matrix<T,2>& tmp_mat=p.scp; tmp_mat=Matrix<T,2>();
            for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){          
                auto data=iterator.Current_Node()._data; TV v_vec({v0(data),v1(data)});
                tmp_mat+=Matrix<T,d>::Outer_Product(iterator.Weight_Gradient(),v_vec);}
            Matrix<T,d> F=p.constitutive_model.Fe;
            const T k_p=(T)1e4; const T saturation=p.saturation; const T eta=p.constitutive_model.eta*saturation;
            tmp_mat=F.Times_Transpose(p.constitutive_model.Times_dP_dF(tmp_mat.Transpose_Times(F))
                                        -eta*k_p*saturation*Times_Cofactor_Matrix_Derivative(F,tmp_mat.Transpose_Times(F)));

            tmp_mat*=p.volume;
            if(false){ high_resolution_clock::time_point te = high_resolution_clock::now();
    	    duration<double> dur = duration_cast<duration<double>>(te-tb);
    	    Log::cout<<"single particle: "<<dur.count()<<std::endl;}}
        if(print_running_time){ high_resolution_clock::time_point te1 = high_resolution_clock::now();
    	duration<double> dur1 = duration_cast<duration<double>>(te1-tb1);
    	Log::cout<<"Collect to scp: "<<dur1.count()<<std::endl;}

        high_resolution_clock::time_point tb2 = high_resolution_clock::now();
        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),f(v));
        if(print_running_time){high_resolution_clock::time_point te2 = high_resolution_clock::now();
	    duration<double> dur2 = duration_cast<duration<double>>(te2-tb2);
	    Log::cout<<"Clear force: "<<dur2.count()<<std::endl;}

        high_resolution_clock::time_point tb3 = high_resolution_clock::now();
#pragma omp parallel for
        for(int tid_process=0;tid_process<threads;++tid_process){
            const Interval<int> thread_x_interval=x_intervals(tid_process);
            for(int tid_collect=0;tid_collect<threads;++tid_collect){
                const Array<int>& index=particle_bins(tid_process,tid_collect);
                for(int i=0;i<index.size();++i){
                    T_Particle& p=particles(index(i));T_INDEX& closest_node=p.closest_node;
                    const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_node(0),thread_x_interval.max_corner-closest_node(0));
            for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){ 
                auto data=iterator.Current_Node()._data;
                TV tmp_vec=p.scp.Transpose_Times(iterator.Weight_Gradient()); 
                f0(data)+=tmp_vec(0); f1(data)+=tmp_vec(1);}}}}

        if (print_running_time){high_resolution_clock::time_point te3 = high_resolution_clock::now();
	    duration<double> dur3 = duration_cast<duration<double>>(te3-tb3);
	    Log::cout<<"Collect to grid: "<<dur3.count()<<std::endl;}
    }

    void Force(Vector<T Struct_type::*,3>& f,const Vector<T Struct_type::*,3>& x) const
    {
        bool print_running_time=false;
        using TV                            = Vector<T,3>;
        using T_INDEX                       = Vector<int,3>;
        using T_Particle                    = MPM_Particle<T,3>;
        using T_Influence_Iterator          = Influence_Iterator<T,3,T_INDEX>;
        using T_Cropped_Influence_Iterator  = Cropped_Influence_Iterator<T,3,T_INDEX>;
        using T_Barrier                     = MPM_Plane_Barrier<T,3>; 
        
        high_resolution_clock::time_point tb1 = high_resolution_clock::now();
        auto v_0=hierarchy.Channel(0,x(0)); auto v_1=hierarchy.Channel(0,x(1)); auto v_2=hierarchy.Channel(0,x(2));
        auto f_0=hierarchy.Channel(0,f(0)); auto f_1=hierarchy.Channel(0,f(1)); auto f_2=hierarchy.Channel(0,f(2));
        const Grid<T,3>& grid=hierarchy.Lattice(0);
#pragma omp parallel for
        for(int i=0;i<simulated_particles.size();++i){
            int id=simulated_particles(i); T_Particle& p=particles(id); Matrix<T,d> tmp_mat; 
            for(T_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),p);iterator.Valid();iterator.Next()){
                auto data=iterator.Current_Node()._data; TV v_vec;
                v_vec(0)=v_0(data); v_vec(1)=v_1(data); v_vec(2)=v_2(data);
                tmp_mat+=Matrix<T,d>::Outer_Product(iterator.Weight_Gradient(),v_vec);}
            Matrix<T,d> F=p.constitutive_model.Fe;
            const T k_p=(T)1e4; const T saturation=p.saturation; const T eta=p.constitutive_model.eta*saturation;
            tmp_mat=F.Times_Transpose(p.constitutive_model.Times_dP_dF(tmp_mat.Transpose_Times(F))
                    -eta*k_p*saturation*Times_Cofactor_Matrix_Derivative(F,tmp_mat.Transpose_Times(F)));
            p.scp=p.volume*tmp_mat;}
        if (print_running_time){high_resolution_clock::time_point te1 = high_resolution_clock::now();
    	duration<double> dur1 = duration_cast<duration<double>>(te1-tb1);
    	Log::cout<<"Collect to scp: "<<dur1.count()<<std::endl;}

        high_resolution_clock::time_point tb2 = high_resolution_clock::now();
        for(int level=0;level<hierarchy.Levels();++level) for(int v=0;v<d;++v) SPGrid::Clear<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),f(v));  
        if (print_running_time){high_resolution_clock::time_point te2 = high_resolution_clock::now();
	    duration<double> dur2 = duration_cast<duration<double>>(te2-tb2);
	    Log::cout<<"Clear force: "<<dur2.count()<<std::endl;}

        high_resolution_clock::time_point tb3 = high_resolution_clock::now();
#pragma omp parallel for
    for(int tid_process=0;tid_process<threads;++tid_process){
        const Interval<int>& thread_x_interval=x_intervals(tid_process);
        for(int tid_collect=0;tid_collect<threads;++tid_collect){
            const Array<int>& index=particle_bins(tid_process,tid_collect);
            for(int i=0;i<index.size();++i){
                T_Particle& p=particles(index(i));T_INDEX& closest_node=p.closest_node;
                const Interval<int> relative_interval=Interval<int>(thread_x_interval.min_corner-closest_node(0),thread_x_interval.max_corner-closest_node(0));
        for(T_Cropped_Influence_Iterator iterator(T_INDEX(-1),T_INDEX(1),relative_interval,p);iterator.Valid();iterator.Next()){ 
            T_INDEX dindex=iterator.Index();T_INDEX current_node=closest_node+iterator.Index();
            if(grid.Node_Indices().Inside(current_node)){ T weight=p.Weight(dindex);
                if(weight>(T)0.){ TV weight_grad=p.Weight_Gradient(dindex); TV tmp_vec=p.scp.Transpose_Times(weight_grad); 
                    f_0(current_node._data)+=tmp_vec(0);f_1(current_node._data)+=tmp_vec(1);f_2(current_node._data)+=tmp_vec(2);}}}}}}
    if (print_running_time){high_resolution_clock::time_point te3 = high_resolution_clock::now();
	duration<double> dur3 = duration_cast<duration<double>>(te3-tb3);
	Log::cout<<"Collect to grid: "<<dur3.count()<<std::endl;}
    }




    void Project(Vector_Base& v) const
    {
        Channel_Vector& v_channels              = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;

        for(int level=0;level<hierarchy.Levels();++level)
            Project_Helper<Struct_type,T,d>(hierarchy.Lattice(0),hierarchy.Allocator(level),hierarchy.Blocks(level),v_channels,barriers);
        // for(int level=0;level<hierarchy.Levels();++level)
        //     Channel_Vector_Traverse_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channels,(unsigned)Node_Saturated);
    }

    double Inner_Product(const Vector_Base& v1,const Vector_Base& v2) const
    {
        const Hierarchy& v1_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v1);
        const Hierarchy& v2_hierarchy           = MPM_CG_Vector<Struct_type,T,d>::Hierarchy(v2);
        Channel_Vector const v1_channels        = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v1).channel;
        Channel_Vector const v2_channels        = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v2).channel;
        assert(&hierarchy == &v1_hierarchy);
        assert(&hierarchy == &v2_hierarchy);

        // for(int level=0;level<hierarchy.Levels();++level)
        //     Channel_Vector_Traverse_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channels,(unsigned)Node_Saturated);
        // Log::cout<<std::endl;
        // for(int level=0;level<hierarchy.Levels();++level)
        //     Channel_Vector_Traverse_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v2_channels,(unsigned)Node_Saturated);

        T result=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level)
            Inner_Product_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channels,v2_channels,result,(unsigned)Node_Saturated);
        return result;
    }

    T Convergence_Norm(const Vector_Base& v) const
    {
        Channel_Vector& v_channels              = MPM_CG_Vector<Struct_type,T,d>::Cg_Vector(v).channel;
        T result=(T)0.;
        
        // for(int level=0;level<hierarchy.Levels();++level)
        //     Channel_Vector_Traverse_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channels,(unsigned)Node_Saturated);
        
        for(int level=0;level<hierarchy.Levels();++level)
            Convergence_Norm_Helper<Struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     v_channels,result,(unsigned)Node_Saturated);
        return result;
    }

    void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const
    {
    }
};
}
#endif

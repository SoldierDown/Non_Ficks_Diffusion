//!#####################################################################
//! \file Multigrid_Solver.h
//!#####################################################################
// Class Multigrid_Solver
//######################################################################
#ifndef __Multigrid_Solver__
#define __Multigrid_Solver__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include "Clear_Mask.h"
#include "Convergence_Norm_Helper.h"
// #include "../Diffusion_Helper/Diffusion_Convergence_Norm_Helper.h"
#include "Copy_Channel.h"
#include "Initialize_Mask.h"
#include "Mark_Boundary.h"
#include "Multigrid_Refinement.h"
#include "Multigrid_Smoother.h"
#include "Initial_Guess_Helper.h"
#include "../Traverse_Helper.h"

namespace Nova{
template<class Base_struct_type,class Multigrid_struct_type,class T,int d>
class Multigrid_Solver
{
    using T_INDEX                               = Vector<int,d>;
    using Hierarchy_Base                        = Grid_Hierarchy<Base_struct_type,T,d>;
    using Base_allocator_type                   = SPGrid::SPGrid_Allocator<Base_struct_type,d>;
    using Base_flag_array_mask                  = typename Base_allocator_type::template Array_mask<unsigned>;
    using Base_block_iterator                   = SPGrid::SPGrid_Block_Iterator<Base_flag_array_mask>;
    using Hierarchy_Multigrid                   = Grid_Hierarchy<Multigrid_struct_type,T,d>;
    using Multigrid_allocator_type              = SPGrid::SPGrid_Allocator<Multigrid_struct_type,d>;
    using Multigrid_flag_array_mask             = typename Multigrid_allocator_type::template Array_mask<unsigned>;

  public:
    Hierarchy_Base& hierarchy;
    Array<Hierarchy_Multigrid*> multigrid_hierarchy;

    T Multigrid_struct_type::* x_channel;
    T Multigrid_struct_type::* b_channel;
    T Multigrid_struct_type::* temp_channel;

    bool FICKS;
    const T dt;
    const T diff_coeff;
    const T Fc;
    const T tau;


    Multigrid_Solver(Hierarchy_Base& hierarchy_input,const int mg_levels,const T FICKS_input,const T dt_input,const T diff_coeff_input,const T Fc_input,const T tau_input)
        :hierarchy(hierarchy_input),FICKS(FICKS_input),dt(dt_input),diff_coeff(diff_coeff_input),Fc(Fc_input),tau(tau_input)
    {
        x_channel                       = &Multigrid_struct_type::ch0;
        b_channel                       = &Multigrid_struct_type::ch1;
        temp_channel                    = &Multigrid_struct_type::ch2;

        // clean up
        for(size_t i=0;i<multigrid_hierarchy.size();++i)
            if(multigrid_hierarchy(i)!=nullptr) delete multigrid_hierarchy(i);
        multigrid_hierarchy.resize(mg_levels);
        multigrid_hierarchy.Fill(nullptr);

        for(int mg_level=0;mg_level<mg_levels;++mg_level)
            multigrid_hierarchy(mg_level)=new Hierarchy_Multigrid(hierarchy.Lattice(0),std::max(hierarchy.Levels(),mg_levels));

        {
            Log::Scope scope("Multigrid_Solver::Activate_Multigrid_Pages");
            // activate page for all multigrid levels
            for(int spgrid_level=0;spgrid_level<hierarchy.Levels();++spgrid_level)
                for(Base_block_iterator iterator(hierarchy.Blocks(spgrid_level));iterator.Valid();iterator.Next_Block()){uint64_t offset=iterator.Offset();
                    uint64_t translated_offset=Multigrid_flag_array_mask::template Translate_Linear_Offset<Base_flag_array_mask>(offset);
                    for(int mg_level=0;mg_level<mg_levels;++mg_level){
                        if(mg_level<=spgrid_level) multigrid_hierarchy(mg_level)->Page_Map(spgrid_level).Set_Page(translated_offset);
                        else{offset=Base_flag_array_mask::DownsampleOffset(offset);
                            translated_offset=Multigrid_flag_array_mask::template Translate_Linear_Offset<Base_flag_array_mask>(offset);
                            multigrid_hierarchy(mg_level)->Page_Map(mg_level).Set_Page(translated_offset);}}}

            for(int i=0;i<mg_levels;++i) multigrid_hierarchy(i)->Update_Block_Offsets();
        }


        {
            Log::Scope scope("Multigrid_Solver::Copy_Interior");
            // copy interior information
            for(int level=0;level<hierarchy.Levels();++level)
                Initialize_Mask<Base_struct_type,Multigrid_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                                            multigrid_hierarchy,level,(unsigned)Cell_Type_Interior);
        }

        {
            Log::Scope scope("Multigrid_Solver::Copy_Dirichlet");
            // copy dirichlet information
            for(int level=0;level<hierarchy.Levels();++level)
                Initialize_Mask<Base_struct_type,Multigrid_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                                            multigrid_hierarchy,level,(unsigned)Cell_Type_Dirichlet);
        }

        // // generate auxiliary flags
        // for(int mg_level=0;mg_level<mg_levels;++mg_level)
        //     Grid_Hierarchy_Initializer<Multigrid_struct_type,T,d>::Flag_Ghost_Cells(*multigrid_hierarchy(mg_level));

        // clear all channels
        for(int i=0;i<mg_levels;++i) for(int level=0;level<multigrid_hierarchy(i)->Levels();++level){
            SPGrid::Clear<Multigrid_struct_type,T,d>(multigrid_hierarchy(i)->Allocator(level),multigrid_hierarchy(i)->Blocks(level),x_channel);
            SPGrid::Clear<Multigrid_struct_type,T,d>(multigrid_hierarchy(i)->Allocator(level),multigrid_hierarchy(i)->Blocks(level),b_channel);
            SPGrid::Clear<Multigrid_struct_type,T,d>(multigrid_hierarchy(i)->Allocator(level),multigrid_hierarchy(i)->Blocks(level),temp_channel);}
    }

    ~Multigrid_Solver()
    {
        for(int i=0;i<(int)multigrid_hierarchy.size();++i) if(multigrid_hierarchy(i)!=nullptr) delete multigrid_hierarchy(i);
    }

    void Initialize_Boundary_Flags(const int boundary_radius,const unsigned mask)
    {
        const int mg_levels=(int)multigrid_hierarchy.size();
        // clear mask
        for(int mg_level=0;mg_level<mg_levels;++mg_level){Hierarchy_Multigrid& current_hierarchy=*multigrid_hierarchy(mg_level);
            for(int level=0;level<current_hierarchy.Levels();++level)
                Clear_Mask<Multigrid_struct_type,T,d>(current_hierarchy.Allocator(level),current_hierarchy.Blocks(level),mask);}

        // set up acceleration structure
        Array<uint64_t> neighbor_offsets;
        const T_INDEX boundary_radius_vector(boundary_radius),zero_vector=T_INDEX();
        for(Range_Iterator<d> iterator(-boundary_radius_vector,boundary_radius_vector);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Index();
            if(index!=zero_vector) neighbor_offsets.Append(Multigrid_flag_array_mask::Linear_Offset(index._data));}

        // mark boundary
        for(int mg_level=0;mg_level<mg_levels;++mg_level){Hierarchy_Multigrid& current_hierarchy=*multigrid_hierarchy(mg_level);
            for(int level=0;level<current_hierarchy.Levels();++level)
                Mark_Boundary<Multigrid_struct_type,T,d>(current_hierarchy,current_hierarchy.Blocks(level),neighbor_offsets,level,mask);
            current_hierarchy.Initialize_Boundary_Blocks(mask);}
    }

    void Initialize_Guess() const
    {
        for(int level=0;level<multigrid_hierarchy(0)->Levels();++level)
            SPGrid::Clear<Multigrid_struct_type,T,d>(multigrid_hierarchy(0)->Allocator(level),multigrid_hierarchy(0)->Blocks(level),x_channel);        
    }

    void Initialize(const int boundary_radius=3)
    {Initialize_Boundary_Flags(boundary_radius,(unsigned)MG_Boundary);}

    void Copy_Channel_Values(T Base_struct_type::* channel_base,T Multigrid_struct_type::* channel_multigrid,const bool copy_to=true) const
    {
        for(int level=0;level<hierarchy.Levels();++level)
            Copy_Channel<Base_struct_type,Multigrid_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),*multigrid_hierarchy(0),
                                                                    channel_base,channel_multigrid,level,copy_to);
    }

    void Initialize_Right_Hand_Side(T Base_struct_type::* channel) const
    {Copy_Channel_Values(channel,b_channel);}

    void Smooth(const int mg_level,const int boundary_smoothing_iterations,const int interior_smoothing_iterations) const
    {
        Log::cout<<"Boundary Smoothing: "<<std::endl;
        Multigrid_Smoother<Multigrid_struct_type,T,d>::Jacobi_Iteration(*multigrid_hierarchy(mg_level),multigrid_hierarchy(mg_level)->Boundary_Blocks(mg_level),
                                                                        mg_level,x_channel,b_channel,temp_channel,boundary_smoothing_iterations,(unsigned)MG_Boundary,FICKS,dt,diff_coeff,Fc,tau);
        Log::cout<<"Interior Smoothing: "<<std::endl;                                                            
        Multigrid_Smoother<Multigrid_struct_type,T,d>::Jacobi_Iteration(*multigrid_hierarchy(mg_level),multigrid_hierarchy(mg_level)->Blocks(mg_level),
                                                                        mg_level,x_channel,b_channel,temp_channel,interior_smoothing_iterations,(unsigned)Cell_Type_Interior,FICKS,dt,diff_coeff,Fc,tau);
        Log::cout<<"Boundary Smoothing: "<<std::endl;
        Multigrid_Smoother<Multigrid_struct_type,T,d>::Jacobi_Iteration(*multigrid_hierarchy(mg_level),multigrid_hierarchy(mg_level)->Boundary_Blocks(mg_level),
                                                                        mg_level,x_channel,b_channel,temp_channel,boundary_smoothing_iterations,(unsigned)MG_Boundary,FICKS,dt,diff_coeff,Fc,tau);
    }

    void V_Cycle(const int boundary_smoothing_iterations,const int interior_smoothing_iterations,const int bottom_smoothing_iterations) const
    {
        const int mg_levels=(int)multigrid_hierarchy.size();
        // mg_levels: 2
        // downstroke
        // i: 0
        for(int i=0;i<mg_levels-1;++i){
            // smooth
            // Smooth(i,boundary_smoothing_iterations,interior_smoothing_iterations);
            // compute residual of level 0
            Compute_Residual(i);    // residual vector is stored in temp_channel
            // restrict residual to level 1
            Multigrid_Refinement<Multigrid_struct_type,T,d>::Restrict(*multigrid_hierarchy(i),*multigrid_hierarchy(i+1),temp_channel,b_channel,Vector<int,2>({i,i+1}));
            // clear x of level 1
            for(int level=0;level<multigrid_hierarchy(i+1)->Levels();++level)
                SPGrid::Clear<Multigrid_struct_type,T,d>(multigrid_hierarchy(i+1)->Allocator(level),multigrid_hierarchy(i+1)->Blocks(level),x_channel);}

        // exact solve
        // mg_levels-1=1: coarsest
        // temp_channel is residual vector restricted from fine mesh
        // Ax0 = b - r
        // Ax' = r
        // x = x0 + x' ==> Ax = b
        Log::cout<<"Before: "<<std::endl;
        Log::cout<<"rhs norm: "<<Convergence_Norm(mg_levels-1,b_channel)<<std::endl;
        Log::cout<<"residual norm: "<<Convergence_Norm(mg_levels-1,temp_channel)<<std::endl;
        Log::cout<<"x norm: "<<Convergence_Norm(mg_levels-1,x_channel)<<std::endl;
        // Solve Ax'=r on level 1
        Multigrid_Smoother<Multigrid_struct_type,T,d>::Exact_Solve(*multigrid_hierarchy(mg_levels-1),x_channel,b_channel,
                                                                   temp_channel,bottom_smoothing_iterations,(unsigned)Cell_Type_Interior,FICKS,dt,diff_coeff,Fc,tau);
        Log::cout<<"After: "<<std::endl;
        Log::cout<<"level 0 residual norm: "<<Convergence_Norm(0,temp_channel)<<std::endl;
        Log::cout<<"level "<<mg_levels-1<<" residual norm: "<<Convergence_Norm(mg_levels-1,temp_channel)<<std::endl;
        Log::cout<<"level 0 x norm: "<<Convergence_Norm(0,x_channel)<<std::endl;
        Log::cout<<"level "<<mg_levels-1<<" x norm: "<<Convergence_Norm(mg_levels-1,x_channel)<<std::endl;
        Log::cout<<"level 0 rhs norm: "<<Convergence_Norm(0,b_channel)<<std::endl;
        Log::cout<<"level "<<mg_levels-1<<" rhs norm: "<<Convergence_Norm(mg_levels-1,b_channel)<<std::endl;
        
        // upstroke
        for(int i=mg_levels-2;i>=0;--i){
            // prolongate
            // (i,i+1) = (0,1)
            // fine: temp_channel
            // coarse: x_channel
            Multigrid_Refinement<Multigrid_struct_type,T,d>::Prolongate(*multigrid_hierarchy(i),*multigrid_hierarchy(i+1),temp_channel,x_channel,Vector<int,2>({i,i+1}));
            // add correction
            for(int level=0;level<multigrid_hierarchy(i)->Levels();++level)
                SPGrid::Masked_Add<Multigrid_struct_type,T,d>(multigrid_hierarchy(i)->Allocator(level),multigrid_hierarchy(i)->Blocks(level),
                                                              x_channel,temp_channel,x_channel,(unsigned)Cell_Type_Interior);
            // propagate ghost values
            // Grid_Hierarchy_Projection<Multigrid_struct_type,T,d>::Propagate_Ghost_Values(*multigrid_hierarchy(i),x_channel);
            // smooth
            // Smooth(i,boundary_smoothing_iterations,interior_smoothing_iterations);
        }
    }

    void Compute_Residual(const int mg_level) const
    {
        Multigrid_Smoother<Multigrid_struct_type,T,d>::Compute_Residual(*multigrid_hierarchy(mg_level),x_channel,b_channel,temp_channel,(unsigned)Cell_Type_Interior,FICKS,dt,diff_coeff,Fc,tau);
    }

    T Convergence_Norm(const int mg_level,T Multigrid_struct_type::* channel) const
    {
        T max_value=(T)0.;

        for(int level=0;level<multigrid_hierarchy(mg_level)->Levels();++level)
            Convergence_Norm_Helper<Multigrid_struct_type,T,d>(multigrid_hierarchy(mg_level)->Allocator(level),multigrid_hierarchy(mg_level)->Blocks(level),
                                                               channel,max_value,(unsigned)Cell_Type_Interior);

        return max_value;
    }
};
}
#endif

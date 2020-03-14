//!#####################################################################
//! \file Restriction_Stencil_Helper.h
//!#####################################################################
// Class Restriction_Stencil_Helper
//######################################################################
#ifndef __Restriction_Stencil_Helper__
#define __Restriction_Stencil_Helper__

#include <nova/SPGrid/Core/SPGrid_Mask.h>
#include <nova/SPGrid/Core/SPGrid_Utilities.h>
#include <nova/Tools/Utilities/Pthread_Queue.h>

namespace Nova{
template<class Struct_type,class T,int d> class Restriction_Stencil_Helper;
//######################################################################
// Class Restriction_Stencil_Thread_Helper
//######################################################################
template<class Struct_type,class T,int d>
struct Restriction_Stencil_Thread_Helper: public Pthread_Queue::Task
{
    Restriction_Stencil_Helper<Struct_type,T,d>* const obj;
    const int index_start,index_end;

    Restriction_Stencil_Thread_Helper(Restriction_Stencil_Helper<Struct_type,T,d>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input)
    {}

    void Run()
    {obj->Run_Index_Range(index_start,index_end);}
};
//######################################################################
// Class Restriction_Stencil_Helper - 2D
//######################################################################
template<class Struct_type,class T>
class Restriction_Stencil_Helper<Struct_type,T,2>
{
    enum{d=2};
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

    T* const coarse_data;
    const T* const fine_data;
    const unsigned* const coarse_mask;
    const uint64_t* const b;        // block offset stream
    const int size;                 // number of blocks to process
    const T scale;
    const unsigned active_flag_mask;

    enum{
        block_xsize = 1u << Flag_array_mask::block_xbits,
        block_ysize = 1u << Flag_array_mask::block_ybits,
        og_xsize = block_xsize+2,
        og_ysize = block_ysize+2,
        xmin = 1,
        ymin = 1,
        // Inclusive!!! give mins and maxs for actual block within shadow grid
        xmax = og_xsize-2,
        ymax = og_ysize-2,
        fine_og_xsize = 2*block_xsize+2,
        fine_og_ysize = 2*block_ysize+2
    };

  public:
    explicit Restriction_Stencil_Helper(T* const coarse_data_input,const T* const fine_data_input,const unsigned* const coarse_mask_input,
                                        const uint64_t* const b_input,const int size_input,const T scale_input,const unsigned active_flag_mask_input)
        :coarse_data(coarse_data_input),fine_data(fine_data_input),coarse_mask(coarse_mask_input),b(b_input),size(size_input),scale(scale_input),
        active_flag_mask(active_flag_mask_input)
    {}

    void Run()
    {Run_Index_Range(0,size-1);}

    static void ComputeFineShadowGrid(uint64_t* offset_grid_ptr,const uint64_t fine_packed_offset)
    {
        typedef uint64_t (&offset_grid_type)[fine_og_xsize][fine_og_ysize];
        offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
        uint64_t current_offset=fine_packed_offset;
        for(int i=0;i<fine_og_xsize;i++){o_grid[i][0]=current_offset;
            current_offset=Flag_array_mask::template Packed_OffsetXdim<1>(current_offset);}
        for(int i=0;i<fine_og_xsize;i++){current_offset=o_grid[i][0];
        for(int j=1;j<fine_og_ysize;j++){current_offset=Flag_array_mask::template Packed_OffsetYdim<1>(current_offset);
            o_grid[i][j]=current_offset;}}
    }

    void Run_Parallel(const int number_of_partitions)
    {
        for(int partition=0;partition<number_of_partitions;partition++){
            // calculate indices of current partition
            int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
            int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
            // create helper object
            Restriction_Stencil_Thread_Helper<Struct_type,T,d>* task=new Restriction_Stencil_Thread_Helper<Struct_type,T,d>(this,first_index_of_partition,last_index_of_partition);
            // enqueue
            pthread_queue->Queue(task);}
        // wait for all tasks to complete
        pthread_queue->Wait();
    }

    void Run_Index_Range(const int index_start,const int index_end)
    {
        // shadow grid declaration
        uint64_t* offset_grid_ptr=(uint64_t*)malloc((fine_og_xsize)*(fine_og_ysize)*sizeof(uint64_t));
        typedef uint64_t (&offset_grid_type)[fine_og_xsize][fine_og_ysize];
        offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
        const uint64_t fine_data_base_addr=reinterpret_cast<uint64_t>(fine_data);
        // iterating through all block indices
        for(int index=index_start;index<=index_end;index++){
            T* const coarse_data_ptr=reinterpret_cast<T*>((uint64_t)coarse_data+b[index]);
            const T* const fine_data_ptr=reinterpret_cast<T*>((uint64_t)fine_data+b[index]);
            const unsigned* const coarse_mask_ptr=reinterpret_cast<unsigned*>((uint64_t)coarse_mask+b[index]);
            const uint64_t packed_offset=(uint64_t)coarse_data_ptr-(uint64_t)coarse_data;
            uint64_t fine_packed_offset=Flag_array_mask::UpsampleOffset(packed_offset);
            fine_packed_offset=Flag_array_mask::template Packed_Offset<-1,-1>(fine_packed_offset);
            ComputeFineShadowGrid(offset_grid_ptr,fine_packed_offset);
            // iterate within block
            int cur_index=0;
            for(int i=xmin;i<=xmax;i++) for(int j=ymin;j<=ymax;j++){
                const unsigned& coarse_mask_value=coarse_mask_ptr[cur_index];
                if(coarse_mask_value&active_flag_mask){T value=(T)0.;
                    const int ii(2*(i-1)),jj(2*(j-1));
                    
                    value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0]));
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1]));
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2]));
                    value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3]));
                    
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0]));
                    value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1]));
                    value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2]));
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3]));
                    
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0]));
                    value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1]));
                    value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2]));
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3]));
                    
                    value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0]));
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1]));
                    value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2]));
                    value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3]));
                    
                    coarse_data_ptr[cur_index]+=(value*scale);}
                cur_index++;}}
        free(offset_grid_ptr);
    }
//######################################################################
};

//######################################################################
// Class Restriction_Stencil_Helper - 3D
//######################################################################
template<class Struct_type,class T>
class Restriction_Stencil_Helper<Struct_type,T,3>
{
    enum{d=3};
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

    T* const coarse_data;
    const T* const fine_data;
    const unsigned* const coarse_mask;
    const uint64_t* const b;        // block offset stream
    const int size;                 // number of blocks to process
    const T scale;
    const unsigned active_flag_mask;

    enum{
        block_xsize = 1u << Flag_array_mask::block_xbits,
        block_ysize = 1u << Flag_array_mask::block_ybits,
        block_zsize = 1u << Flag_array_mask::block_zbits,
        og_xsize = block_xsize+2,
        og_ysize = block_ysize+2,
        og_zsize = block_zsize+2,
        xmin = 1,
        ymin = 1,
        zmin = 1,
        // Inclusive!!! give mins and maxs for actual block within shadow grid
        xmax = og_xsize-2,
        ymax = og_ysize-2,
        zmax = og_zsize-2,
        fine_og_xsize = 2*block_xsize+2,
        fine_og_ysize = 2*block_ysize+2,
        fine_og_zsize = 2*block_zsize+2
    };

  public:
    explicit Restriction_Stencil_Helper(T* const coarse_data_input,const T* const fine_data_input,const unsigned* const coarse_mask_input,
                                        const uint64_t* const b_input,const int size_input,const T scale_input,const unsigned active_flag_mask_input)
        :coarse_data(coarse_data_input),fine_data(fine_data_input),coarse_mask(coarse_mask_input),b(b_input),size(size_input),scale(scale_input),
        active_flag_mask(active_flag_mask_input)
    {}

    void Run()
    {Run_Index_Range(0,size-1);}

    static void ComputeFineShadowGrid(uint64_t* offset_grid_ptr,const uint64_t fine_packed_offset)
    {
        typedef uint64_t (&offset_grid_type)[fine_og_xsize][fine_og_ysize][fine_og_zsize];
        offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
        uint64_t current_offset=fine_packed_offset;
        for(int i=0;i<fine_og_xsize;i++){o_grid[i][0][0]=current_offset;
            current_offset=Flag_array_mask::template Packed_OffsetXdim<1>(current_offset);}
        for(int i=0;i<fine_og_xsize;i++){current_offset=o_grid[i][0][0];
        for(int j=1;j<fine_og_ysize;j++){current_offset=Flag_array_mask::template Packed_OffsetYdim<1>(current_offset);
            o_grid[i][j][0]=current_offset;}}
        for(int i=0;i<fine_og_xsize;i++)
        for(int j=0;j<fine_og_ysize;j++){current_offset=o_grid[i][j][0];
        for(int k=1;k<fine_og_zsize;k++){current_offset=Flag_array_mask::template Packed_OffsetZdim<1>(current_offset);
            o_grid[i][j][k]=current_offset;}}
    }

    void Run_Parallel(const int number_of_partitions)
    {
        for(int partition=0;partition<number_of_partitions;partition++){
            // calculate indices of current partition
            int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
            int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
            // create helper object
            Restriction_Stencil_Thread_Helper<Struct_type,T,d>* task=new Restriction_Stencil_Thread_Helper<Struct_type,T,d>(this,first_index_of_partition,last_index_of_partition);
            // enqueue
            pthread_queue->Queue(task);}
        // wait for all tasks to complete
        pthread_queue->Wait();
    }

    void Run_Index_Range(const int index_start,const int index_end)
    {
        // shadow grid declaration
        uint64_t* offset_grid_ptr=(uint64_t*)malloc((fine_og_xsize)*(fine_og_ysize)*(fine_og_zsize)*sizeof(uint64_t));
        typedef uint64_t (&offset_grid_type)[fine_og_xsize][fine_og_ysize][fine_og_zsize];
        offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
        const uint64_t fine_data_base_addr=reinterpret_cast<uint64_t>(fine_data);
        // iterating through all block indices
        for(int index=index_start;index<=index_end;index++){
            T* const coarse_data_ptr=reinterpret_cast<T*>((uint64_t)coarse_data+b[index]);
            const T* const fine_data_ptr=reinterpret_cast<T*>((uint64_t)fine_data+b[index]);
            const unsigned* const coarse_mask_ptr=reinterpret_cast<unsigned*>((uint64_t)coarse_mask+b[index]);
            const uint64_t packed_offset=(uint64_t)coarse_data_ptr-(uint64_t)coarse_data;
            uint64_t fine_packed_offset=Flag_array_mask::UpsampleOffset(packed_offset);
            fine_packed_offset=Flag_array_mask::template Packed_Offset<-1,-1,-1>(fine_packed_offset);
            ComputeFineShadowGrid(offset_grid_ptr,fine_packed_offset);
            // iterate within block
            int cur_index=0;
            for(int i=xmin;i<=xmax;i++) for(int j=ymin;j<=ymax;j++) for(int k=zmin;k<=zmax;k++){
                const unsigned& coarse_mask_value=coarse_mask_ptr[cur_index];
                if(coarse_mask_value&active_flag_mask){T value=(T)0.;
                    const int ii(2*(i-1)),jj(2*(j-1)),kk(2*(k-1));
                    
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+0]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+1]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+2]));
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+3]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+3]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+3]));
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+0]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+1]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+2]));
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+3]));
                    
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+3]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+0]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+1]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+2]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+3]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+0]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+1]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+2]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+3]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+3]));
                    
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+3]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+0]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+1]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+2]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+3]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+0]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+1]));
                    value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+2]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+3]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+3]));
                    
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+0]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+1]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+2]));
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+3]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+3]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+0]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+1]));
                    value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+2]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+3]));
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+0]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+1]));
                    value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+2]));
                    value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+3]));
                    
                    coarse_data_ptr[cur_index]+=(value*scale);}
                cur_index++;}}
        free(offset_grid_ptr);
    }
//######################################################################
};
}
#endif

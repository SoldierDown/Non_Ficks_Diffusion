//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Hierarchy/Visualization/Grid_Hierarchy_Visualization.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Log/Log.h>
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Random_Numbers/Random_Numbers.h>
#include <nova/Tools/Utilities/Constants.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Pthread_Queue.h>
#include "../MPM_Data.h"
#include "../Multigrid_Solver/Multigrid_Data.h"
#include <omp.h>

using namespace Nova;
using namespace SPGrid;

extern Pthread_Queue* pthread_queue;

namespace Nova{
int number_of_threads=0;
}

int main(int argc,char** argv)
{
    enum {d=2};
    typedef float T;typedef Vector<T,d> TV;
    typedef Vector<int,d> T_INDEX;

    using Struct_type                           = MPM_Data<T>;
    using Multigrid_struct_type                 = Multigrid_Data<T>;
    using Allocator_type                        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper                       = Grid_Topology_Helper<Flag_array_mask>;
    using Channel_Vector                        = Vector<T Struct_type::*,d>;
    using Hierarchy                             = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Lookup                      = Grid_Hierarchy_Lookup<Struct_type,T,d>;
    using Hierarchy_Rasterizer                  = Hierarchical_Rasterizer<Struct_type,T,d>;
    using Hierarchy_Visualization               = Grid_Hierarchy_Visualization<Struct_type,T>;
    enum {number_of_faces_per_cell              = Topology_Helper::number_of_faces_per_cell};

    Log::Initialize_Logging();

    Log::Finish_Logging();
    return 0;
}

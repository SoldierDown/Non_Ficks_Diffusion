//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include <nova/Tools/Utilities/Pthread_Queue.h>
#include "../Smoke_Driver.h"
#include "Standard_Tests/Standard_Tests.h"
using namespace Nova;

namespace Nova{
int number_of_threads=0;
}

extern Pthread_Queue* pthread_queue;

int main(int argc,char** argv)
{
    enum {d=2};
    typedef float T;typedef Vector<T,d> TV;
    typedef Vector<int,d> T_INDEX;

    Smoke_Example<T,d> *example=new Standard_Tests<T,d>();
    example->Parse(argc,argv);

    if(number_of_threads) pthread_queue=new Pthread_Queue(number_of_threads);

    File_Utilities::Create_Directory(example->output_directory);
    File_Utilities::Create_Directory(example->output_directory+"/common");
    File_Utilities::Create_Directory(example->output_directory+"/density_data");
    Log::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    Smoke_Driver<T,d> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}

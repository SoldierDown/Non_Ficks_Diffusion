//!#####################################################################
//! \file main.cpp
//!#####################################################################
#include "../MPM_Driver.h"
#include "Standard_Tests/Standard_Tests.h"
using namespace Nova;

int main(int argc,char** argv)
{
    enum {d=2};
    typedef float T;

    MPM_Example<T,d> *example=new Standard_Tests<T,d>();
    example->Parse(argc,argv);
    example->bbox=Range<T,d>(example->domain.max_corner,example->domain.min_corner);
    File_Utilities::Create_Directory(example->output_directory);
    File_Utilities::Create_Directory(example->output_directory+"/common");
    Log::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    MPM_Driver<T,d> driver(*example);
    driver.Execute_Main_Program();
       
    delete example;

    return 0;
}

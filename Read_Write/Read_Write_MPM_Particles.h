//!#####################################################################
//! \file Read_Write_MPM_Particles.h
//!#####################################################################
// Class Read_Write_MPM_Particles
//###################################################################### 
#ifndef __Read_Write_MPM_Particles__
#define __Read_Write_MPM_Particles__

#include <nova/Tools/Arrays/Array.h>
#include <nova/Tools/Read_Write/Utilities/Read_Write.h>
#include "../MPM_Particle.h"

namespace Nova{
template<class T,int d>
class Read_Write<Array<MPM_Particle<T,d>>>
{
    using TV                        = Vector<T,d>;
    using T_Particle                = MPM_Particle<T,d>;

  public:
    static void Read(std::istream& input,Array<T_Particle>& object)
    {
        size_t number_of_points;
        Read_Write<size_t>::Read(input,number_of_points);
        object.resize(number_of_points);
        for(size_t i=0;i<number_of_points;++i){
            Read_Write<TV>::Read(input,object(i).X);
            Read_Write<TV>::Read(input,object(i).V);
            Read_Write<T>::Read(input,object(i).mass);
            Read_Write<T>::Read(input,object(i).volume);}
    }

    static void Write(std::ostream& output,const Array<T_Particle>& object)
    {
        int valid_particles=0;
        for(size_t i=0;i<object.size();++i) if(object(i).valid) valid_particles++;
        Read_Write<size_t>::Write(output,valid_particles);
        for(size_t i=0;i<object.size();++i){
            if(object(i).valid){
                Read_Write<TV>::Write(output,object(i).X);
                Read_Write<TV>::Write(output,object(i).V);
                Read_Write<T>::Write(output,object(i).mass);
                Read_Write<T>::Write(output,object(i).volume);
            }
}
    }
};
}
#endif

//!#####################################################################
//! \file Particle_Indicator_Helper.h
//!#####################################################################
// Class Particle_Indicator_Helper
//######################################################################
#ifndef __Particle_Indicator_Helper__
#define __Particle_Indicator_Helper__

#include <iostream>
#include <stdlib.h>
#include <string>
#include "MPM_Particle.h"
#include <nova/Tools/Vectors/Vector.h>

using namespace std;
namespace Nova{
template<class T,int d>
class Particle_Indicator_Helper
{
    using T_Particle                    = MPM_Particle<T,d>;
  public:
    Particle_Indicator_Helper(const Array<T_Particle>& particles,const std::string& filename)
    {Run(particles,filename);}

    void Run(const Array<T_Particle>& particles,const std::string& filename) const
    {
        FILE* fp=fopen(filename.c_str(),"w");
        fprintf(fp, "ID\tVALID\tFLUID\n");
        for(int i=0;i<particles.size();++i) fprintf(fp, "%d\t%d\t%d\n",i,particles(i).valid,particles(i).eos);
        fclose(fp);
    }

};
}
#endif

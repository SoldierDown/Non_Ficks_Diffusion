//!#####################################################################
//! \file Volume2Mesh_Helper.h
//!#####################################################################
// Class Volume2Mesh_Helper
//######################################################################
#ifndef __Volume2Mesh_Helper__
#define __Volume2Mesh_Helper__
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Dynamics/Rigid_Bodies/Rigid_Body.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include <vector>
#include <assert.h>
#include <openvdb/openvdb.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/math/Vec4.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <iostream>
using namespace std;
using namespace openvdb;
namespace Nova{
template<class T>
class Volume2Mesh_Helper
{
    public: 
    Parse_Args *parse_args;
    std::string output_directory;
    int first_frame, last_frame;
    int step;
    int thickness;
    T threshold;
    T density_factor;
  public:
    Volume2Mesh_Helper()
        :parse_args(nullptr)
    {}

    ~Volume2Mesh_Helper()
    {}

    void Convert_Frame(const int current_frame)
    {
		openvdb::initialize();
        string vdb_filename=output_directory+"/"+std::to_string(current_frame)+".vdb";
        openvdb::io::File file(vdb_filename);
        file.open(); openvdb::GridBase::Ptr base_grid=file.readGrid("density"); file.close();
        openvdb::FloatGrid::Ptr mygrid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
        vector<math::Vec3<T>> points; vector<math::Vec3<uint32_t>> triangles; vector<math::Vec4<uint32_t>> quads;
        tools::volumeToMesh(*mygrid,points,triangles,quads,threshold);
        Save_To_OBJ(current_frame,points,triangles,quads);
    }

    void Register_Options()
    {
        if(!parse_args) return;
        parse_args->Add_String_Argument("-o","","output directory");
        parse_args->Add_Integer_Argument("-first_frame",0,"first frame");
        parse_args->Add_Integer_Argument("-last_frame",0,"last frame");
        parse_args->Add_Integer_Argument("-step",1,"step");
        parse_args->Add_Double_Argument("-th",.5,"density threshold");
    }

    void Parse_Options()
    {
        if(!parse_args) return;
        output_directory=parse_args->Get_String_Value("-o");
        first_frame=parse_args->Get_Integer_Value("-first_frame");
        last_frame=parse_args->Get_Integer_Value("-last_frame");
        step=parse_args->Get_Integer_Value("-step");
        threshold=parse_args->Get_Double_Value("-th");        
    }

    void Parse(int argc,char* argv[])
    {
        parse_args=new Parse_Args;
        Register_Options();    
        parse_args->Parse(argc,argv);
        Parse_Options();
    }

    void Save_To_OBJ(const int current_frame,vector<math::Vec3<T>>& points,vector<math::Vec3<uint32_t>>& triangles,vector<math::Vec4<uint32_t>>& quads)
    {
        std::cout<<"saving to obj"<<std::endl;
        std::cout<<"points: "<<points.size()<<", triangles: "<<triangles.size()<<", quads: "<<quads.size()<<std::endl;
        string obj_filename=output_directory+"/obj_data/"+std::to_string(current_frame)+".obj";
		// const char* file_name = name_s.c_str();
		FILE* fp=fopen(obj_filename.c_str(), "w");
		for (int i=0;i<points.size();++i)
		    fprintf(fp, "v %.6f %.6f %.6f\n",points[i](0),points[i](1),points[i](2));
        for(int i=0;i<triangles.size();++i)
		    fprintf(fp, "f %d %d %d\n",triangles[i](0)+1,triangles[i](1)+1,triangles[i](2)+1);
		for (int i=0;i<quads.size();++i)
		    fprintf(fp, "f %d %d %d %d\n",quads[i](0)+1,quads[i](1)+1,quads[i](2)+1,quads[i](3)+1);
		fclose(fp);
    }
};
}
#endif

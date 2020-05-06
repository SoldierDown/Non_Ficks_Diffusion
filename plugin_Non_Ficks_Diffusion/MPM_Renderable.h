//!#####################################################################
//! \file MPM_Renderable.h
//!#####################################################################
// Class MPM_Renderable
//######################################################################
#ifndef __MPM_Renderable_h__
#define __MPM_Renderable_h__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

#include "ViewportManager.h"
#include "plugins/Simulation/Simulation_Renderable.h"

#include "Viewer_Data.h"
#include "../MPM_Particle.h"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

namespace Nova{
template<class T,int d>
class MPM_Renderable: public Simulation_Renderable<T,d>
{
    using TV                                = Vector<T,d>;
    using Base                              = Simulation_Renderable<T,d>;
    using T_INDEX                           = Vector<int,d>;
    using T_Particle                        = MPM_Particle<T,d>;
    using Index_type                        = std::array<SPGrid::ucoord_t,d>;
    using Struct_type                       = Viewer_Data<T>;
    using Hierarchy                         = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type                    = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                   = typename Allocator_type::template Array_mask<unsigned>;
    using T_Range_Iterator                  = Range_Iterator<d,T_INDEX>;

    using Base::directory_name;using Base::active_frame;using Base::_app;

    struct Voxel
    {
        glm::vec3 location;
        float dx;

        Voxel() {}
        
        Voxel(const TV& location_input,const T dx_input)
        {
            for(int v=0;v<d;++v) location[v]=location_input[v];
            if(d==2) location[2]=1;
            dx=dx_input;
        }
    };

    Hierarchy *hierarchy;
    std::vector<Voxel> voxels;
    unsigned elements_per_block;
    int levels;
    Index_type block_size;
    Grid<T,d> grid;
    Array<T_Particle> particles;
    std::vector<glm::vec3> solid_locations,fluid_locations;
    unsigned int VASO,VAFO,VAVO;
    unsigned int VBSO,VBFO,VBVO;
    std::string frame_title;
    bool selected,draw_voxels,draw_solid,draw_fluid;

  public:
    MPM_Renderable(ApplicationFactory& app,const std::string& directory_name,int max_frame)
        :Base(app,directory_name,max_frame),hierarchy(nullptr),elements_per_block(0),levels(0),selected(false),draw_voxels(false),draw_solid(true),draw_fluid(true)
    {
        std::istream *input=(d==3)?File_Utilities::Safe_Open_Input(directory_name+"/common/hierarchy.struct3d"):File_Utilities::Safe_Open_Input(directory_name+"/common/hierarchy.struct2d");
        Read_Write<int>::Read(*input,levels);
        std::cout<<"Levels: "<<levels<<std::endl;
        Read_Write<unsigned>::Read(*input,elements_per_block);
        std::cout<<"# of elements/block: "<<elements_per_block<<std::endl;
        for(int v=0;v<d;++v) Read_Write<unsigned>::Read(*input,block_size[v]);
        delete input;

        File_Utilities::Read_From_File(directory_name+"/common/fine_grid",grid);
        app.GetIOService().On("DISPLAY-VOXELS",[&](IOEvent& event){this->Display_Voxels();});
        app.GetIOService().On("DISPLAY-SOLID",[&](IOEvent& event){this->Display_Solid();});
        app.GetIOService().On("DISPLAY-FLUID",[&](IOEvent& event){this->Display_Fluid();});
    }

    virtual ~MPM_Renderable()
    {if(hierarchy!=nullptr) delete hierarchy;}

    void Display_Voxels()
    {draw_voxels=!draw_voxels;}

    void Display_Solid()
    {draw_solid=!draw_solid;}

    void Display_Fluid()
    {draw_fluid=!draw_fluid;}

    void Clear_Buffers()
    {
        if(VASO) glDeleteVertexArrays(1,&VASO);
        if(VAFO) glDeleteVertexArrays(1,&VAFO);
        if(VAVO) glDeleteVertexArrays(1,&VAVO);
        if(VBSO) glDeleteBuffers(1,&VBSO);
        if(VBFO) glDeleteBuffers(1,&VBFO);
        if(VBVO) glDeleteBuffers(1,&VBVO);

        voxels.clear();
        particles.Clear();
        solid_locations.clear();
        fluid_locations.clear();
    }

    virtual void Initialize_Buffers()
    {
        glGenVertexArrays(1,&VASO);
        glGenBuffers(1,&VBSO);
        glBindVertexArray(VASO);
        glBindBuffer(GL_ARRAY_BUFFER,VBSO);
        glBufferData(GL_ARRAY_BUFFER,solid_locations.size()*sizeof(glm::vec3),&solid_locations[0],GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(glm::vec3),(GLvoid*)0);
        glBindVertexArray(0);

        glGenVertexArrays(1,&VAFO);
        glGenBuffers(1,&VBFO);
        glBindVertexArray(VAFO);
        glBindBuffer(GL_ARRAY_BUFFER,VBFO);
        glBufferData(GL_ARRAY_BUFFER,fluid_locations.size()*sizeof(glm::vec3),&fluid_locations[0],GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(glm::vec3),(GLvoid*)0);
        glBindVertexArray(0);

        glGenVertexArrays(1,&VAVO);
        glGenBuffers(1,&VBVO);
        glBindVertexArray(VAVO);
        glBindBuffer(GL_ARRAY_BUFFER,VBVO);
        glBufferData(GL_ARRAY_BUFFER,voxels.size()*sizeof(Voxel),&voxels[0],GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(Voxel),(GLvoid*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1,1,GL_FLOAT,GL_FALSE,sizeof(Voxel),(GLvoid*)(sizeof(glm::vec3)));
        glBindVertexArray(0);
    }

    virtual void Load_Active_Frame() override
    {
        Clear_Buffers();

        if(hierarchy!=nullptr) delete hierarchy;
        hierarchy=new Hierarchy(grid,levels);

        File_Utilities::Read_From_File(directory_name+"/"+std::to_string(active_frame)+"/particles",particles);
        std::istream *id_file=File_Utilities::Safe_Open_Input(directory_name+"/particle_indicator/"+std::to_string(active_frame)+".txt",false);
        std::string line;std::getline(*id_file,line);
        std::vector<bool> is_fluid(particles.size(),false);
        for(size_t i=0;i<particles.size();++i){size_t id;int fluid;
            *id_file>>id>>fluid>>fluid;
            if(fluid) is_fluid[i]=true;}
        delete id_file;

        for(size_t i=0;i<particles.size();++i){
            if(particles[i].density>(T)0.){glm::vec3 X;
                for(int v=0;v<d;++v) X[v]=particles[i].X[v];
                if(d==2) X[2]=1;
                if(!is_fluid[i]) solid_locations.push_back(X);
                else fluid_locations.push_back(X);}}

#if 0
        std::stringstream ribs;ribs<<"particles/hydrogel_";
        if(active_frame<10) ribs<<"0"<<active_frame;
        else ribs<<active_frame;
        ribs<<".rib";
        std::fstream outFile(ribs.str(),std::ios::out);
        for(size_t i=0;i<locations.size();++i){
            outFile<<"                AttributeBegin"<<std::endl;
            outFile<<"                        Bxdf \"PxrSurface\" \"particles\" \"string __materialid\" [\"particles_SG\"]"<<std::endl;
            outFile<<"                        Translate "<<locations[i][0]<<" "<<locations[i][1]<<" "<<locations[i][2]<<std::endl;
            outFile<<"                        Sphere 0.001 -0.001 0.001 360"<<std::endl;
            outFile<<"                AttributeEnd"<<std::endl;}
        outFile.close();
#endif

        std::stringstream ss;ss<<directory_name<<"/"<<active_frame;
        std::istream* input1=File_Utilities::Safe_Open_Input(ss.str()+"/flags");
        std::istream* input2=File_Utilities::Safe_Open_Input(ss.str()+"/block_offsets");

        for(int level=0;level<levels;++level){unsigned number_of_blocks=0;
            Read_Write<unsigned>::Read(*input2,number_of_blocks);

            for(unsigned block=0;block<number_of_blocks;++block){T_INDEX base_index;
                Read_Write<T_INDEX>::Read(*input2,base_index);
                T_Range_Iterator range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);

                for(unsigned e=0;e<elements_per_block;++e){unsigned flag;T phi;
                    const T_INDEX index=base_index+range_iterator.Index();
                    uint64_t offset=Flag_array_mask::Linear_Offset(index._data);
                    Read_Write<unsigned>::Read(*input1,flag);
                    hierarchy->template Set<unsigned>(level,&Struct_type::flags).Mask(offset,flag);

                    TV location=hierarchy->Lattice(level).Center(index);
                    glm::vec3 X;for(int v=0;v<d;++v) X[v]=location[v];
                    if(d==2) X[2]=1;

                    if(flag&(Cell_Type_Interior|Cell_Type_Dirichlet)) voxels.push_back(Voxel(location,hierarchy->Lattice(level).dX[0]));
                    range_iterator.Next();}}}

        std::cout<<"# of voxels: "<<voxels.size()<<std::endl;

        delete input2;delete input1;
        hierarchy->Update_Block_Offsets();

        Initialize_Buffers();
    }

    virtual void draw() override
    {
        Base::draw();

        glm::mat4 projection,view,model;
        view = _app.GetWorld().Get_ViewMatrix();
        model = _app.GetWorld().Get_ModelMatrix();
        projection = _app.GetWorld().Get_ProjectionMatrix();
        auto slicePlanes = _app.GetWorld().Get_Slice_Planes();

        if(solid_locations.size()>0 && draw_solid){
            auto shader = _app.GetShaderManager().GetShader("Points");
            shader->SetMatrix4("projection",projection);
            shader->SetMatrix4("view",view);
            shader->SetMatrix4("model",model);
            shader->SetVector4f("slice_plane0",slicePlanes[0]);
            shader->SetVector4f("slice_plane1",slicePlanes[1]);
            glBindVertexArray(VASO);
            glEnable(GL_PROGRAM_POINT_SIZE);
            shader->SetVector3f("point_color",glm::vec3(254./256.,110./256.,0.));
            glDrawArrays(GL_POINTS,0,solid_locations.size());
            glBindVertexArray(0);}

        if(fluid_locations.size()>0 && draw_fluid){
            auto shader = _app.GetShaderManager().GetShader("Points");
            shader->SetMatrix4("projection",projection);
            shader->SetMatrix4("view",view);
            shader->SetMatrix4("model",model);
            shader->SetVector4f("slice_plane0",slicePlanes[0]);
            shader->SetVector4f("slice_plane1",slicePlanes[1]);
            glBindVertexArray(VAFO);
            glEnable(GL_PROGRAM_POINT_SIZE);
            shader->SetVector3f("point_color",glm::vec3(153./256.,204./256.,255./256.));
            glDrawArrays(GL_POINTS,0,fluid_locations.size());
            glBindVertexArray(0);}

        if(draw_voxels){auto shader = _app.GetShaderManager().GetShader("Grid");
            shader->SetMatrix4("projection",projection);
            shader->SetMatrix4("view",view);
            shader->SetMatrix4("model",model);
            shader->SetVector4f("slice_plane0",slicePlanes[0]);
            shader->SetVector4f("slice_plane1",slicePlanes[1]);
            shader->SetVector3f("grid_color", glm::vec3( .9, .9, .9 ));
            shader->SetFloat("size_multiplier", 1.0f);
            glBindVertexArray(VAVO);
            glDrawArrays(GL_POINTS,0,voxels.size());
            glBindVertexArray(0);}

        if( _app.GetWorld().GetViewportManager().IsPrimaryViewport() )
            _app.GetTextRenderingService() << frame_title;
    }
};
}
#endif

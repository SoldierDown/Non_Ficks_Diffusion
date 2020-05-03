//!#####################################################################
//! \file Autonomous_Navigation_Renderable.h
//!#####################################################################
// Class Autonomous_Navigation_Renderable
//######################################################################
#ifndef __Autonomous_Navigation_Renderable__
#define __Autonomous_Navigation_Renderable__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Lookup.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Dynamics/Rigid_Bodies/Rigid_Body.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

#include "ViewportManager.h"
#include "plugins/Simulation/Simulation_Renderable.h"

#include "Viewer_Data.h"
#include "../Poisson_Solver/Grid_Hierarchy_Projection.h"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

namespace Nova{
template<class T,int d>
class Autonomous_Navigation_Renderable: public Simulation_Renderable<T,d>
{
    using TV                                = Vector<T,d>;
    using BOUNDARY_ELEMENT_INDEX            = Vector<int,d>;
    using Base                              = Simulation_Renderable<T,d>;
    using T_INDEX                           = Vector<int,d>;
    using Index_type                        = std::array<SPGrid::ucoord_t,d>;
    using Struct_type                       = Viewer_Data<T>;
    using Channel_Vector                    = Vector<T Struct_type::*,d>;
    using Hierarchy                         = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Lookup                  = Grid_Hierarchy_Lookup<Struct_type,T,d>;
    using Allocator_type                    = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask                   = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper                   = Grid_Topology_Helper<Flag_array_mask>;
    using T_Range_Iterator                  = Range_Iterator<d,T_INDEX>;
    using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;

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

    struct Color_Voxel
    {
        glm::vec3 location,color;

        Color_Voxel() {}

        Color_Voxel(const TV& location_input,const glm::vec3& color_input)
        {
            for(int v=0;v<d;++v) location[v]=location_input[v];
            if(d==2) location[2]=1;
            color=color_input;
        }
    };

    Hierarchy *hierarchy;
    std::vector<Voxel> voxels;
    std::vector<Color_Voxel> density_voxels;
    Voxel selected_voxel;
    unsigned elements_per_block;
    int selected_voxel_level,levels;
    bool draw_density,draw_velocity,draw_surface,draw_voxels;
    Index_type block_size;
    Grid<T,d> grid;
    T_INDEX selected_voxel_index;
    unsigned int VAO,VBO;
    unsigned int VASO,VBSO;
    unsigned int VADO,VBDO;
    unsigned int VABO,VBBO,EBBO;
    std::string frame_title;

    T Struct_type::* density_channel;
    T Struct_type::* divergence_channel;
    Channel_Vector face_velocity_channels;
    Rigid_Body<T,d> body;

    glm::vec3 _slice_start;
    float _dx;
    glm::vec3 _max_slice;
    glm::vec4 _bounding_sphere;

  public:
    Autonomous_Navigation_Renderable(ApplicationFactory& app,const std::string& directory_name,int max_frame)
        :Base(app,directory_name,max_frame),hierarchy(nullptr),elements_per_block(0),levels(0),draw_density(false),
        draw_velocity(false),draw_voxels(true)
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

        density_channel                         = &Struct_type::ch0;
        face_velocity_channels(0)               = &Struct_type::ch1;
        face_velocity_channels(1)               = &Struct_type::ch2;
        if(d==3) face_velocity_channels(2)      = &Struct_type::ch3;
        divergence_channel                      = &Struct_type::ch4;
    }

    virtual ~Autonomous_Navigation_Renderable()
    {if(hierarchy!=nullptr) delete hierarchy;}

    void Display_Voxels()
    {draw_voxels=!draw_voxels;}

    void Clear_Buffers()
    {
        if(VAO) glDeleteVertexArrays(1,&VAO);
        if(VASO) glDeleteVertexArrays(1,&VASO);
        if(VADO) glDeleteVertexArrays(1,&VADO);
        if(VBO) glDeleteBuffers(1,&VBO);
        if(VBSO) glDeleteBuffers(1,&VBSO);
        if(VBDO) glDeleteBuffers(1,&VBDO);
        if(VABO) glDeleteVertexArrays(1,&VABO);
        if(VBBO) glDeleteBuffers(1,&VBBO);
        if(EBBO) glDeleteBuffers(1,&EBBO);

        voxels.clear();
        density_voxels.clear();
        selected_voxel_level=-1;
    }

    virtual void Initialize_Buffers()
    {
        glGenVertexArrays(1,&VAO);
        glGenBuffers(1,&VBO);
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER,VBO);
        glBufferData(GL_ARRAY_BUFFER,voxels.size()*sizeof(Voxel),&voxels[0],GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(Voxel),(GLvoid*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1,1,GL_FLOAT,GL_FALSE,sizeof(Voxel),(GLvoid*)(sizeof(glm::vec3)));
        glBindVertexArray(0);

        glGenVertexArrays(1,&VASO);
        glGenBuffers(1,&VBSO);
        glBindVertexArray(VASO);
        glBindBuffer(GL_ARRAY_BUFFER,VBSO);
        glBufferData(GL_ARRAY_BUFFER,sizeof(Voxel),&selected_voxel,GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(Voxel),(GLvoid*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1,1,GL_FLOAT,GL_FALSE,sizeof(Voxel),(GLvoid*)(sizeof(glm::vec3)));
        glBindVertexArray(0);

        glGenVertexArrays(1,&VADO);
        glGenBuffers(1,&VBDO);
        glBindVertexArray(VADO);
        glBindBuffer(GL_ARRAY_BUFFER,VBDO);
        glBufferData(GL_ARRAY_BUFFER,density_voxels.size()*sizeof(Color_Voxel),&density_voxels[0],GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(Color_Voxel),(GLvoid*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(Color_Voxel),(GLvoid*)(sizeof(glm::vec3)));
        glBindVertexArray(0);

#if 0
        glGenVertexArrays(1,&VABO);
        glGenBuffers(1,&VBBO);
        glGenBuffers(1,&EBBO);
        glBindVertexArray(VABO);
        glBindBuffer(GL_ARRAY_BUFFER,VBBO);
        glBufferData(GL_ARRAY_BUFFER,body.object->points.size()*sizeof(TV),&body.object->points[0],GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0,d,GL_FLOAT,GL_FALSE,sizeof(TV),(GLvoid*)0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,EBBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,body.object->elements.size()*sizeof(BOUNDARY_ELEMENT_INDEX),&body.object->elements[0],GL_STATIC_DRAW);
        glBindVertexArray(0);
#endif
    }

    void Draw_Surface(Shader& shader)
    {
        glBindVertexArray(VABO);
        glDrawElements(GL_TRIANGLES,body.object->elements.size()*3,GL_UNSIGNED_INT,0);
        glBindVertexArray(0);
    }

    virtual void Load_Active_Frame() override
    {
        Clear_Buffers();

        if(hierarchy!=nullptr) delete hierarchy;
        hierarchy=new Hierarchy(grid,levels);

        std::stringstream ss;ss<<directory_name<<"/"<<active_frame;
        std::istream* input1=File_Utilities::Safe_Open_Input(ss.str()+"/flags");
        std::istream* input2=File_Utilities::Safe_Open_Input(ss.str()+"/block_offsets");

        // check for density
        draw_density=File_Utilities::File_Exists(ss.str()+"/spgrid_density",false);
        std::istream* input3=draw_density?File_Utilities::Safe_Open_Input(ss.str()+"/spgrid_density"):nullptr;

        // check for velocity
        draw_velocity=File_Utilities::File_Exists(ss.str()+"/spgrid_u",false);
        Vector<std::istream*,d> velocity_input(nullptr);
        if(draw_velocity){velocity_input(0)=File_Utilities::Safe_Open_Input(ss.str()+"/spgrid_u");
            velocity_input(1)=File_Utilities::Safe_Open_Input(ss.str()+"/spgrid_v");
            if(d==3) velocity_input(2)=File_Utilities::Safe_Open_Input(ss.str()+"/spgrid_w");}

        draw_surface=File_Utilities::File_Exists(ss.str()+"/body.rgd",false);
        if(draw_surface){File_Utilities::Read_From_File(ss.str()+"/body.rgd",body);
            // transform all points to world space
            for(int i=0;i<(int)body.object->points.size();++i){
                TV X=body.World_Space_Point(body.object->points[i]);
                body.object->points[i]=X;}}

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

                    if(draw_density){T density;
                        Read_Write<T>::Read(*input3,density);
                        hierarchy->Allocator(level).template Get_Array<Struct_type,T>(density_channel)(offset)=density;
                        if(flag&Cell_Type_Interior && density>(T).01) density_voxels.push_back(Color_Voxel(location,glm::vec3(1,1,1)*std::max(density,(T)0.)));}

                    if(draw_velocity) for(int axis=0;axis<d;++axis){T velocity;
                        Read_Write<T>::Read(*velocity_input(axis),velocity);
                        TV face_location=hierarchy->Lattice(level).Face(axis,index);
                        hierarchy->Allocator(level).template Get_Array<Struct_type,T>(face_velocity_channels(axis))(offset)=velocity;}

                    if(flag&(Cell_Type_Interior|Cell_Type_Dirichlet)) voxels.push_back(Voxel(location,hierarchy->Lattice(level).dX[0]));
                    range_iterator.Next();}}}

        delete input2;delete input1;
        if(input3!=nullptr) delete input3;
        if(draw_velocity) for(int v=0;v<d;++v) if(velocity_input(v)!=nullptr) delete velocity_input(v);
        hierarchy->Update_Block_Offsets();

        if(draw_velocity){
            for(int level=0;level<levels;++level)
                SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),divergence_channel);
            Hierarchy_Projection::Compute_Divergence(*hierarchy,face_velocity_channels,divergence_channel);}

        Initialize_Buffers();

        _slice_start.x = grid.domain.min_corner[0]+ grid.dX[0]/2.0f;
        _slice_start.y = grid.domain.min_corner[1]+ grid.dX[1]/2.0f;
        _slice_start.z = (d==3)?grid.domain.min_corner[2]+ grid.dX[2]/2.0f:1;
        _dx = grid.dX[0];
        _max_slice.x = grid.counts[0]-1; // We slice on cells, so one less than the node count
        _max_slice.y = grid.counts[1]-1;
        _max_slice.z = (d==3)?grid.counts[2]-1:1;
        
        glm::vec3 sphere_center;
        sphere_center.x = grid.domain.min_corner[0];
        sphere_center.y = grid.domain.min_corner[1];
        sphere_center.z = (d==3)?grid.domain.min_corner[2]:1;
        sphere_center = glm::mix(sphere_center,_max_slice*_dx,0.5f);
        float sphere_radius = glm::distance(sphere_center,glm::vec3(_max_slice*_dx));
        _bounding_sphere = glm::vec4(sphere_center.x,sphere_center.y,sphere_center.z,sphere_radius);
    }

    virtual void draw() override
    {
        Base::draw();

        glm::mat4 projection,view,model;
        view = _app.GetWorld().Get_ViewMatrix();
        model = _app.GetWorld().Get_ModelMatrix();
        projection = _app.GetWorld().Get_ProjectionMatrix();
        auto slicePlanes = _app.GetWorld().Get_Slice_Planes();
        glm::vec3 cameraPos = _app.GetWorld().GetCameraPosition();
        std::vector<glm::vec3> lights = _app.GetWorld().GetSceneLightPositions();

        if(d==2 && draw_surface){
            auto shader = _app.GetShaderManager().GetShader("BasicColored");
            shader->SetMatrix4("projection",projection);
            shader->SetMatrix4("view",view);
            shader->SetMatrix4("model",model);
            shader->SetVector4f("slice_plane0",slicePlanes[0]);
            shader->SetVector4f("slice_plane1",slicePlanes[1]);
            shader->SetVector3f("basecolor",glm::vec3(228/255.0f,26/255.0f,28/255.0f));
            shader->SetInteger("enable_slice",1);

            glBindVertexArray(VABO);
            glDrawElements(GL_LINES,body.object->elements.size()*d,GL_UNSIGNED_INT,0);
            glBindVertexArray(0);}

        if(d==3 && draw_surface){
            if(_app.GetWorld().SolidMode()){
                auto _shader = _app.GetShaderManager().GetShader("BasicMeshShader");
                _shader->SetMatrix4("projection",projection);
                _shader->SetMatrix4("view",view);
                _shader->SetMatrix4("model",model);
                _shader->SetVector4f("slice_plane0",slicePlanes[0]);
                _shader->SetVector4f("slice_plane1",slicePlanes[1]);
                _shader->SetInteger("shaded",_app.GetWorld().LightingMode()?1:2);
                _shader->SetVector3f("cameraPos",cameraPos);
                _shader->SetVector3f("defaultChannels.diffuse",glm::vec3(77/255.0f,175/255.0f,74/255.0f));
                _shader->SetFloat("defaultChannels.shininess",10.0f);
                _shader->SetInteger("activeLights",std::min((int)lights.size(),4));         // no more than 4 lights.
                for(auto l:lights) _shader->SetVector3f("lights[0].position",l);
                Draw_Surface(*(_shader.get()));}

            if(_app.GetWorld().WireframeMode()){
                auto _shader = _app.GetShaderManager().GetShader("BasicColored");
                _shader->SetMatrix4("projection",projection);
                _shader->SetMatrix4("view",view);
                _shader->SetMatrix4("model",model);
                _shader->SetVector4f("slice_plane0",slicePlanes[0]);
                _shader->SetVector4f("slice_plane1",slicePlanes[1]);
                _shader->SetVector3f("basecolor",glm::vec3(0.8,0.9,0.8));
                _shader->SetInteger("enable_slice",1);
                
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
                glLineWidth(1);
                glEnable(GL_POLYGON_OFFSET_LINE);
                glPolygonOffset(-1,-1);
                Draw_Surface(*(_shader.get()));
                glDisable(GL_POLYGON_OFFSET_LINE);
                glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);}}

        if(draw_density){
                auto shader = _app.GetShaderManager().GetShader("Density");
                shader->SetMatrix4("projection",projection);
                shader->SetMatrix4("view",view);
                shader->SetMatrix4("model",model);
                shader->SetVector4f("slice_plane0",slicePlanes[0]);
                shader->SetVector4f("slice_plane1",slicePlanes[1]);
                glBindVertexArray(VADO);
                glEnable(GL_PROGRAM_POINT_SIZE);
                glDrawArrays(GL_POINTS,0,density_voxels.size());
                glBindVertexArray(0);}

        if(draw_voxels){auto shader = _app.GetShaderManager().GetShader("Grid");
            shader->SetMatrix4("projection",projection);
            shader->SetMatrix4("view",view);
            shader->SetMatrix4("model",model);
            shader->SetVector4f("slice_plane0",slicePlanes[0]);
            shader->SetVector4f("slice_plane1",slicePlanes[1]);
            shader->SetVector3f("grid_color", glm::vec3( .9, .9, .9 ));
            shader->SetFloat("size_multiplier", 1.0f);
            glBindVertexArray(VAO);
            glDrawArrays(GL_POINTS,0,voxels.size());
            glBindVertexArray(0); 

            if(selected_voxel_level != -1){
                glBindVertexArray(VASO);
                shader->SetVector3f("grid_color",glm::vec3(.9,.9,0));
                shader->SetFloat("size_multiplier",.98f);
                glDrawArrays(GL_POINTS,0,1);
                glBindVertexArray(0);}}

        if(_app.GetWorld().GetViewportManager().IsPrimaryViewport())
            _app.GetTextRenderingService() << frame_title;
    }

    virtual bool selectable()
    {
        unsigned int interactionViewport = _app.GetWorld().GetViewportManager().InteractionViewportId();
        if(_app.GetWorld().SliceMode(interactionViewport) || d==2) return true;
        return false;
    }

    virtual float hit_test(glm::vec3 start_point,glm::vec3 end_point)
    {
        unsigned int interactionViewport = _app.GetWorld().GetViewportManager().InteractionViewportId();
        auto slicePlanes = _app.GetWorld().Get_Slice_Planes(interactionViewport);
        glm::vec3 normal1=(d==3)?glm::vec3(slicePlanes[0]):glm::vec3(0,0,1);
        glm::vec3 normal2=(d==3)?glm::vec3(slicePlanes[1]):glm::vec3(0,0,-1);
        glm::vec3 e = glm::normalize(end_point-start_point);

        if((fabs(glm::dot(normal1, e)) < 1e-6) || (fabs(glm::dot(normal2, e)) < 1e-6 )) return 0.0f;

        float t1,t2;
        if(d==3)
        {
            t1 = -(glm::dot(normal1,start_point) + slicePlanes[0][3])/glm::dot(normal1,e);
            t2 = -(glm::dot(normal2,start_point) + slicePlanes[1][3])/glm::dot(normal2,e);
        }
        else
        {
            t1 = -(glm::dot(normal1,start_point) + -(1-_dx*.5))/glm::dot(normal1,e);
            t2 = -(glm::dot(normal2,start_point) + 1-_dx*.5)/glm::dot(normal2,e);
        }

        glm::vec3 intersection=start_point+(t1+t2)*0.5f*e;
        TV X;for(int v=0;v<d;++v) X[v]=intersection[v];
        bool valid=grid.domain.Inside(X);
        if(!valid) return -1.0f;

        if(d==2)
        {
            t1 = -(glm::dot(normal1,start_point) + -(1-grid.dX[0]*.5))/glm::dot(normal1,e);
            t2 = -(glm::dot(normal2,start_point) + 1-grid.dX[0]*.5)/glm::dot(normal2,e);
        }
        
        return (t1+t2)*0.5f;
    }

    virtual glm::vec4 bounding_sphere()
    {return _bounding_sphere;}

    virtual void assign_selection(glm::vec3 start_point,glm::vec3 end_point,glm::vec3 intersection)
    {
        uint64_t offset;TV X,weights;
        for(int v=0;v<d;++v) X[v]=intersection[v];

        if(grid.domain.Inside(X) && Hierarchy_Lookup::Cell_Lookup(*hierarchy,X,offset,selected_voxel_level,weights)){
            selected_voxel_index=T_INDEX(Flag_array_mask::LinearToCoord(offset));
            selected_voxel=Voxel(hierarchy->Lattice(selected_voxel_level).Center(selected_voxel_index),hierarchy->Lattice(selected_voxel_level).dX[0]);
            uint64_t selected_voxel_offset=Flag_array_mask::Linear_Offset(selected_voxel_index._data);

            std::stringstream ss;ss<<std::endl<<"Selected Cell: "<<std::endl;
            ss<<"\tIndex: "<<selected_voxel_index<<", Level: "<<selected_voxel_level<<std::endl;
            unsigned flag=hierarchy->Allocator(selected_voxel_level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags)(selected_voxel_offset);

            if(flag&Cell_Type_Dirichlet) ss<<"\tDirichlet"<<std::endl;
            else if(flag&Cell_Type_Interior) ss<<"\tInterior"<<std::endl;
            if(flag&MG_Boundary) ss<<"\tBoundary"<<std::endl;

            if(draw_density){T density=hierarchy->Allocator(selected_voxel_level).template Get_Const_Array<Struct_type,T>(density_channel)(selected_voxel_offset);
                ss<<"\tDensity: "<<density<<std::endl;}

            if(draw_velocity){for(int axis=0;axis<d;++axis){const unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
                ss<<"Axis "<<axis<<": Left Face: ";
                if(hierarchy->template Set<unsigned>(selected_voxel_level,&Struct_type::flags).Is_Set(selected_voxel_offset,face_valid_mask))
                    ss<<hierarchy->Allocator(selected_voxel_level).template Get_Const_Array<Struct_type,T>(face_velocity_channels(axis))(selected_voxel_offset);
                else ss<<"N/A";
                ss<<", Right Face: ";
                const uint64_t other_offset=Flag_array_mask::Packed_Add(selected_voxel_offset,Topology_Helper::Axis_Vector_Offset(axis));
                if(hierarchy->template Set<unsigned>(selected_voxel_level,&Struct_type::flags).Is_Set(other_offset,face_valid_mask))
                    ss<<hierarchy->Allocator(selected_voxel_level).template Get_Const_Array<Struct_type,T>(face_velocity_channels(axis))(other_offset)<<std::endl;
                else ss<<"N/A"<<std::endl;}

                ss<<"Divergence: "<<hierarchy->Channel(selected_voxel_level,divergence_channel)(selected_voxel_offset)<<std::endl;}

            frame_title=ss.str();

            glBindVertexArray(VASO);
            glBindBuffer(GL_ARRAY_BUFFER,VBSO);
            glBufferData(GL_ARRAY_BUFFER,sizeof(Voxel),&selected_voxel,GL_STATIC_DRAW);
            glBindVertexArray(0);}
        else unassign_selection();
    }

    virtual void unassign_selection()
    {
        selected_voxel_level = -1;
        frame_title="";
    }

    virtual bool slice_provider()
    {return true;}

    virtual glm::vec3 slice_start()
    {return _slice_start;}

    virtual float slice_dx()
    {return _dx;}
    
    virtual glm::vec3 max_slice()
    {return _max_slice;}
};
}
#endif

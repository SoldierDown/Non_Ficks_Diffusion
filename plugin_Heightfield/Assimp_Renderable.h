//!#####################################################################
//! \file Assimp_Renderable.h
//!#####################################################################
// Class Assimp_Renderable
//######################################################################
#ifndef __Assimp_Renderable_h__
#define __Assimp_Renderable_h__

#include "ViewportManager.h"
#include "plugins/Assimp/AssimpRenderable_Model.h"
#include "plugins/Simulation/Simulation_Renderable.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <fstream>
#include <iostream>

namespace Nova{
template<class T,int d>
class Assimp_Renderable: public Simulation_Renderable<T,d>
{
    using BASE  = Simulation_Renderable<T,d>;

    using BASE::directory_name;using BASE::active_frame;using BASE::_app;

    std::unique_ptr<AssimpRenderable_Model> _model;
    bool selected;

  public:
    Assimp_Renderable(ApplicationFactory& app,const std::string& directory_name,int max_frame)
        :BASE(app,directory_name,max_frame),selected(false)
    {
        _model=std::unique_ptr<AssimpRenderable_Model>(new AssimpRenderable_Model(app));
    }

    virtual ~Assimp_Renderable() {}

    virtual void Load_Active_Frame() override
    {
        std::cout<<"Loading model for frame: "<<active_frame<<std::endl;
        if(_model!=nullptr) _model.reset(new AssimpRenderable_Model(_app));
        std::stringstream ss;ss<<directory_name<<"/"<<active_frame<<"/heightfield.obj";
        _model->Load_Model(ss.str());
    }

    virtual void draw() override
    {
        BASE::draw();

        glm::mat4 projection,view,model;
        view = _app.GetWorld().Get_ViewMatrix();
        model = _app.GetWorld().Get_ModelMatrix();
        projection = _app.GetWorld().Get_ProjectionMatrix();
        auto slicePlanes = _app.GetWorld().Get_Slice_Planes();
        glm::vec3 cameraPos = _app.GetWorld().GetCameraPosition();
        std::vector<glm::vec3> lights = _app.GetWorld().GetSceneLightPositions();
        
        if(_app.GetWorld().SolidMode())
        {
            auto _shader = _app.GetShaderManager().GetShader("BasicMeshShader");
            _shader->SetMatrix4("projection",projection);
            _shader->SetMatrix4("view",view);
            _shader->SetMatrix4("model",model);
            _shader->SetVector4f("slice_plane0",slicePlanes[0]);
            _shader->SetVector4f("slice_plane1",slicePlanes[1]);
            _shader->SetInteger("selected",selected?1:0);
            _shader->SetInteger("shaded",_app.GetWorld().LightingMode()?1:0);
            _shader->SetVector3f("cameraPos",cameraPos);
            _shader->SetInteger("activeLights", std::min( (int)lights.size(),4));       // no more than 4 lights.
            for(auto l : lights) _shader->SetVector3f("lights[0].position",l);
            _model->Draw(*(_shader.get()));
        }

        if(_app.GetWorld().WireframeMode())
        {
            auto _shader = _app.GetShaderManager().GetShader("BasicColored");
            _shader->SetMatrix4("projection",projection);
            _shader->SetMatrix4("view",view);
            _shader->SetMatrix4("model",model);
            _shader->SetVector4f("slice_plane0",slicePlanes[0]);
            _shader->SetVector4f("slice_plane1",slicePlanes[1]);
            _shader->SetVector3f("basecolor",glm::vec3(0.8,0.8,0.9));
            _shader->SetInteger("enable_slice",1);
            
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
            glLineWidth(1);
            glEnable(GL_POLYGON_OFFSET_LINE);
            glPolygonOffset(-1,-1);
            _model->Draw(*(_shader.get()));
            glDisable(GL_POLYGON_OFFSET_LINE);
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        }
    }

    virtual bool selectable()
    {return true;}

    virtual float hit_test(glm::vec3 start_point,glm::vec3 end_point)
    {return _model->TestIntersection(start_point,end_point);}

    virtual glm::vec4 bounding_sphere()
    {return _model->BoundingSphere();}

    virtual void assign_selection(glm::vec3 start_point,glm::vec3 end_point,glm::vec3 intersection)
    {selected=true;}

    virtual void unassign_selection()
    {selected=false;}
};
}
#endif

#include <igl/opengl/glfw/Viewer.h>
#include <thread>
#include "PhysicsHook.h"
#include "GooHook.h"
#include <igl/unproject.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

static PhysicsHook *hook = NULL;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
        hook->run();
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    hook->reset();
}

bool drawCallback(igl::opengl::glfw::Viewer &viewer)
{
    if (!hook)
        return false;

    hook->render(viewer);
    return false;
}


bool keyCallback(igl::opengl::glfw::Viewer &viewer, unsigned int key, int modifiers)
{
	if (key == ' ')
	{
		toggleSimulation();
		return true;
	}
	Eigen::Vector4f look4 = viewer.core.view.inverse() * Eigen::Vector4f(0, 0, 1.0, 0.0);
	Eigen::Vector4f left4 = viewer.core.view.inverse() * Eigen::Vector4f(1.0, 0.0, 0.0, 0.0);
	Eigen::Vector3f look(look4[0], look4[1], look4[2]);
	Eigen::Vector3f left(left4[0], left4[1], left4[2]);
	if (key == 'w')
	{
		viewer.core.camera_base_translation += look;
	}
	if (key == 's')
	{
		viewer.core.camera_base_translation -= look;
	}
	if (key == 'a')
	{
		viewer.core.camera_base_translation += left;
	}
	if (key == 'd')
	{
		viewer.core.camera_base_translation -= left;
	}
	return false;
}

bool mouseCallback(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    Eigen::Vector3f pos(viewer.down_mouse_x, viewer.down_mouse_y, 0);
    Eigen::Matrix4f model = viewer.core.view;
    Eigen::Vector3f unproj = igl::unproject(pos, model, viewer.core.proj, viewer.core.viewport);
    hook->mouseClicked(unproj[0], -unproj[1], button);

    return true;
}

bool mouseScroll(igl::opengl::glfw::Viewer& viewer, float delta)
{
    return true;
}


bool drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Control", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Run/Pause Sim", ImVec2(-1, 0)))
        {
            toggleSimulation();
        }
        if (ImGui::Button("Reset Sim", ImVec2(-1, 0)))
        {
            resetSimulation();
        }
		if (ImGui::Button("Test Force and Gradient", ImVec2(-1, 0)))
		{
            auto gooHook = static_cast<GooHook*>(hook);
            std::cout<<"Test Force"<<std::endl;
            gooHook->testForceDifferential();
            //std::cout<<std::endl<<"Test Step Projection"<<std::endl;
            //gooHook->testStepProjection();
            //std::cout<<std::endl<<"Test Rigid rod constraints and gradient"<<std::endl;
            //gooHook->testConstriantAndGradient();
            //std::cout<<std::endl<<"Test Lagrange Multiple"<<std::endl;
            //gooHook->testLagrangeMultiple();

		}
    }

	
    hook->drawGUI(menu);
    return false;
}

int main(int argc, char *argv[])
{
  igl::opengl::glfw::Viewer viewer;

  hook = new GooHook();
  hook->reset();
  viewer.core.orthographic = true;
  viewer.core.camera_zoom = 4.0;
  viewer.data().show_lines = false;
  viewer.data().set_face_based(false);
  viewer.core.is_animating = true;
  viewer.callback_key_pressed = keyCallback;
  viewer.callback_pre_draw = drawCallback;
  viewer.callback_mouse_down = mouseCallback;
  viewer.callback_mouse_scroll = mouseScroll;

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  menu.callback_draw_viewer_menu = [&]() {drawGUI(menu); };
  viewer.launch();
}

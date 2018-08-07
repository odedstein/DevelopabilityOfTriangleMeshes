/*
 
 2018, Oded Stein, Eitan Grinspun and Keenan Crane
 
 This file is part of the code for "Developability of Triangle Meshes".
 
 The code for "Developability of Triangle Meshes" is free software: you can
 redistribute it and/or modify it under the terms of the GNU General Public
 License as published by the Free Software Foundation, either version 2 of the
 License, or (at your option) any later version.
 
 The code for "Developability of Triangle Meshes" is distributed in the hope
 that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with the code for "Developability of Triangle Meshes". If not,
 see <https://www.gnu.org/licenses/>.
 
 */


#ifndef IGL_OVIEWER_H
#define IGL_OVIEWER_H

#ifndef IGL_OPENGL_4
#define IGL_OPENGL_4
#endif


#ifdef _WIN32
#  include <windows.h>
#  undef max
#  undef min
#endif

#include <chrono>
#include <thread>

#ifdef __APPLE__
#   include <OpenGL/gl3.h>
#   define __gl_h_ /* Prevent inclusion of the old gl.h */
#else
#   include <GL/gl.h>
#endif

#include <igl/get_seconds.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>


#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>

namespace igl
{
    
    class OViewer
    {
    public:
        IGL_INLINE OViewer();
        IGL_INLINE OViewer(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
        IGL_INLINE OViewer(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& E);
        
        IGL_INLINE void set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
        IGL_INLINE void set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& E);
        IGL_INLINE int launch();
        
        IGL_INLINE void set_colors(const Eigen::Vector3d& color);
        IGL_INLINE void set_colors(const Eigen::MatrixXd& color);
        IGL_INLINE void set_colors(const Eigen::MatrixXd& color_amb, const Eigen::MatrixXd& color_spec, const Eigen::MatrixXd& color_dif);
        IGL_INLINE void set_linetransparency(const Eigen::VectorXd& transp);
        IGL_INLINE void set_face_based_coloring(const bool& colorMode);
        IGL_INLINE bool get_face_based_coloring();
        
        IGL_INLINE const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& get_V() const;
        IGL_INLINE void modify_V(const Eigen::MatrixXd& V);
        IGL_INLINE void modify_F(const Eigen::MatrixXi& F, Eigen::MatrixXi E = Eigen::MatrixXi());
        
        IGL_INLINE int get_width() const;
        IGL_INLINE int get_height() const;
        IGL_INLINE void set_width(const int& w);
        IGL_INLINE void set_height(const int& h);
        IGL_INLINE void set_window_size(const int& w, const int& h);
        
        enum ProjectionType {
            Orthographic,
            Frustum
        };
        
        IGL_INLINE ProjectionType get_projection_type() const;
        IGL_INLINE void set_projection_type(const ProjectionType& type);
        
        enum DisplayMode {
            Fill,
            FillWireframe,
            Wireframe,
            Texture,
            TextureWireframe
        };
        
        IGL_INLINE DisplayMode get_display_mode() const;
        IGL_INLINE void set_display_mode(const DisplayMode& mode);
        
        IGL_INLINE bool get_dragging_enabled() const;
        IGL_INLINE void set_dragging_enabled(const bool& mode);
        
        IGL_INLINE Eigen::Vector3f get_background_color() const;
        IGL_INLINE void set_background_color(const Eigen::Vector3d& color);
        IGL_INLINE Eigen::Vector3f get_wire_color() const;
        IGL_INLINE void set_wire_color(const Eigen::Vector3d& color);
        IGL_INLINE Eigen::Vector3f get_points_outline_color() const;
        IGL_INLINE void set_points_outline_color(const Eigen::Vector3d& color);
        IGL_INLINE Eigen::Vector3f get_points_fill_color() const;
        IGL_INLINE void set_points_fill_color(const Eigen::Vector3d& color);
        IGL_INLINE Eigen::Vector3f get_points_important_color() const;
        IGL_INLINE void set_points_important_color(const Eigen::Vector3d& color);
        IGL_INLINE Eigen::Vector3f get_points_selected_color() const;
        IGL_INLINE void set_points_selected_color(const Eigen::Vector3d& color);
        
        IGL_INLINE void set_point_size(const float& size);
        IGL_INLINE void set_point_size(const float& size, const float& importantSize, const float& selectedSize);
        
        IGL_INLINE float get_camera_zoom() const;
        IGL_INLINE void set_camera_zoom(const float& zoom);
        IGL_INLINE Eigen::Quaternion<float> get_camera_rotation() const;
        IGL_INLINE void set_camera_rotation(const Eigen::Quaternion<float>& rotation);
        IGL_INLINE Eigen::Vector3f get_camera_translation() const;
        IGL_INLINE void set_camera_translation(const Eigen::Vector3f& translation);
        
        IGL_INLINE bool get_rotation_enabled() const;
        IGL_INLINE void set_rotation_enabled(const bool& enabled);
        
        IGL_INLINE Eigen::Vector3f get_light_position() const;
        IGL_INLINE void set_light_position(const Eigen::Vector3f& position);
        IGL_INLINE float get_lighting_factor() const;
        IGL_INLINE void set_lighting_factor(const float& factor);
        IGL_INLINE float get_shininess() const;
        IGL_INLINE void set_shininess(const float& shiny);
        IGL_INLINE bool get_cel_shading() const;
        IGL_INLINE void set_cel_shading(const bool& cel);
        
        IGL_INLINE void mark_points(const Eigen::VectorXi& points);
        IGL_INLINE const Eigen::VectorXi& get_marked_points() const;
        IGL_INLINE void set_important_points(const Eigen::VectorXi& points);
        IGL_INLINE const Eigen::VectorXi& get_important_points() const;
        IGL_INLINE const Eigen::VectorXi& get_selected_points() const;
        IGL_INLINE const Eigen::MatrixXd& get_marked_points_translations() const;
        
        IGL_INLINE const Eigen::Matrix4f& get_model() const;
        IGL_INLINE const Eigen::Matrix4f& get_view() const;
        IGL_INLINE const Eigen::Matrix4f& get_proj() const;
        
        IGL_INLINE void stop_dragging();
        IGL_INLINE bool get_currently_dragging() const;
        
        IGL_INLINE void set_plottitles(const std::vector<std::string>& titles);
        IGL_INLINE void set_plotlines(const std::vector<Eigen::VectorXd>& lines);
        IGL_INLINE void set_plots(const std::vector<std::string>& titles, const std::vector<Eigen::VectorXd>& lines);
        IGL_INLINE void set_plotcolors(const std::vector<Eigen::Vector3d>& colors);
        IGL_INLINE void set_plotting_enabled(const bool& mode);
        IGL_INLINE bool get_plotting_enabled();
        IGL_INLINE void set_plotting_logplot(const bool& mode);
        IGL_INLINE bool get_plotting_logplot();
        
        IGL_INLINE void update_projection_matrix();
        IGL_INLINE void update_view_matrix();
        IGL_INLINE void update_model_matrix();
        IGL_INLINE void update_lighting();
        IGL_INLINE void update_colors();
        IGL_INLINE void update_linetransp();
        IGL_INLINE void update_texture();
        IGL_INLINE void update_points();
        IGL_INLINE void update_plots();
        
        IGL_INLINE void exit();
        
        //Callback, called once per frame. The bool signifies if any modification to the mesh is happening.
        void (*predraw_callback)(OViewer&, bool) = NULL;
        //Callback, called after internal keyboard handling is done. Arguments as in GLFW function
        void (*keyboard_callback)(OViewer&, int, int, int, int) = NULL;
        
        double highdpi = 1.;
        IGL_INLINE void handle_selection(float x, float y);
        IGL_INLINE void snap();
        IGL_INLINE void add_translation_to_selected(const Eigen::Vector3d& diff);
        
        
        //Navigation variables
        bool CTRL_pressed = false;
        bool ALT_pressed = false;
        bool SHIFT_pressed = false;
        bool LMB_pressed = false;
        bool RMB_pressed = false;
        bool MMB_pressed = false;
        bool MOD_happened = false;
        
        double mouse_x, mouse_y, down_mouse_x, down_mouse_y, down_mouse_z=0.;
        Eigen::Quaternion<float> down_rot;
        Eigen::Vector3f down_translation;
        Eigen::MatrixXd marked_down_translation;
        
        
    private:
        GLFWwindow* window;
        
        double start_t, last_t;
        
        Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> m_V, m_V_mod, m_VN, m_VN_mod, m_facedVN_mod, m_amb, m_spec, m_dif;
        Eigen::VectorXf m_linetransp;
        bool linetranspEnabled = false;
        Eigen::Matrix<GLuint, Eigen::Dynamic, 3, Eigen::RowMajor> m_F, m_F_mod;
        Eigen::Matrix<GLuint, Eigen::Dynamic, 2, Eigen::RowMajor> m_E, m_E_mod;
        
        Eigen::VectorXi markedPoints;
        Eigen::VectorXi pointIsImportant;
        Eigen::VectorXi pointIsSelected;
        Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> pointsCoords;
        
        bool launched = false;
        bool meshSet = false;
        bool meshBuffersGenerated = false, colorsBuffersGenerated = false, linetranspBufferGenerated = false, pointsBuffersGenerated = false, plotBuffersGenerated = false;
        
        IGL_INLINE bool showLines() const;
        IGL_INLINE bool showFill() const;
        IGL_INLINE bool showTexture() const;
        
        int windowWidth = 800;
        int windowHeight = 600;
        
        float near = 1.;
        float far = 100.;
        float cameraAngle = 35./360.*M_PI;
        Eigen::Vector3f camera_eye;
        Eigen::Vector3f camera_center;
        Eigen::Vector3f camera_up;
        
        Eigen::Affine3f initialShift;
        
        Eigen::Vector3f lightPosition;
        float lightingFactor = 0.8;
        float shininess = 35.;
        bool celShading = false;
        
        Eigen::Matrix4f model, view, proj;
        
        float cameraZoom = 1.;
        Eigen::Quaternion<float> cameraRotation;
        Eigen::Vector3f cameraTranslation;
        bool rotationEnabled = true;
        
        GLuint fill_prog_id, texture_prog_id, lines_prog_id, points_prog_id, plot_prog_id;
        GLuint VAO;
        GLuint VBO, edgebasedVBO, NBO, ambBO, specBO, difBO, linetranspBO, FBO, EBO, PBO, importantPBO, selectedPBO, plotLBO, plotBBO;
        
        Eigen::Vector3f backgroundColor;
        Eigen::Vector3d defaultColor; bool defaultColorUsed = true;
        Eigen::Vector3f wireColor;
        Eigen::Vector3f pointsOutlineColor;
        Eigen::Vector3f pointsFillColor;
        Eigen::Vector3f pointsImportantColor;
        Eigen::Vector3f pointsSelectedColor;
        float pointSize=20., importantPointSize=20., selectedPointSize=27.;
        
        bool faceBasedColoring = false;
        bool draggingEnabled = false;
        bool currentlyDragging = false;
        
        std::vector<std::string> plotTitles;
        std::vector<Eigen::VectorXd> plotLines;
        std::vector<Eigen::Vector3f> plotColors;
        Eigen::Vector3f plotBGColor;
        bool plottingEnabled = false;
        bool plottingLogPlot = false;
        
        Eigen::MatrixXd markedPointsTranslations;
        
        ProjectionType projectionType = Frustum;
        DisplayMode displayMode = FillWireframe;
        
        IGL_INLINE void defaults();
        
        IGL_INLINE void update_mesh_data(const bool& updateFacesEdges = true);
        
        IGL_INLINE GLuint start_program(const std::string& vertexShader, const std::string& fragmentShader, const std::string& geometryShader = std::string());
        IGL_INLINE void drawing();
        IGL_INLINE void draw_plots();
        
        
#ifdef IGL_VIEWER_WITH_NANOGUI
        
        //Text renderer
        IGL_INLINE void render_text(const Eigen::Vector2f& position, const std::string& text);
        IGL_INLINE void render_text(const Eigen::Vector2f& position, const std::string& text, const Eigen::Vector3f& color);
        struct NVGcontext;
        NVGcontext *ctx;
        bool textRendererInitialized = false;
#endif
        
    };
    
}


#ifndef IGL_STATIC_LIBRARY
#  include "OViewer.cpp"
#endif

#endif

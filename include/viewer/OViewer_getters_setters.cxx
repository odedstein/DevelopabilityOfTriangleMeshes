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


//All getters/setters in the OVIewer

IGL_INLINE void OViewer::set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    Eigen::MatrixXi E;
    igl::edges(F, E);
    set_mesh(V, F, E);
}

IGL_INLINE void OViewer::set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& E)
{
    assert(V.cols()==3 && "Every vertex must have 3 coordinates.\n");
    m_V = V.cast<float>();
    m_V_mod = m_V;
    
    //Center and scale
    float scale = 1.;
    Eigen::Vector3f transl;
    if(V.rows() != 0) {
        Eigen::MatrixXd barycenters;
        igl::barycenter(V, F, barycenters);
        Eigen::Vector3d minPoint = barycenters.colwise().minCoeff();
        Eigen::Vector3d maxPoint = barycenters.colwise().maxCoeff();
        Eigen::Vector3d center = 0.5*(minPoint + maxPoint);
        transl = -center.cast<float>();
        scale = 2. / (maxPoint-minPoint).array().abs().maxCoeff();
    }
    initialShift = Eigen::Affine3f::Identity();
    initialShift.scale(scale);
    initialShift.translate(transl);
    
    assert(F.cols()==3 && "Every face must have 3 vertices.\n");
    m_F = F.cast<GLuint>();
    m_F_mod = m_F;
    
    per_face_normals(m_V, m_F_mod, m_facedVN_mod);
    per_vertex_normals(m_V, m_F_mod, m_facedVN_mod, m_VN);
    m_VN_mod = m_VN;
    
    assert(E.cols()==2 && "Every edge must have 2 vertices.\n");
    m_E = E.cast<GLuint>();
    m_E_mod = m_E;
    
    meshSet = true;
    
    if(launched) {
        update_mesh_data(true);
        if(defaultColorUsed) {
            set_colors(defaultColor);
            defaultColorUsed = true;
        }
        linetranspEnabled = false;
    }
}


IGL_INLINE void OViewer::set_colors(const Eigen::Vector3d& color)
{
    Eigen::MatrixXd colors(faceBasedColoring ? m_F_mod.rows() : m_V_mod.rows(), 3);
    for(int i=0; i<colors.rows(); ++i)
        colors.row(i) = color;
    
    set_colors(colors);
}

IGL_INLINE void OViewer::set_colors(const Eigen::MatrixXd& color)
{
    set_colors(0.4*color, 0.3 + 0.2*(color.array()-0.3), color);
}

IGL_INLINE void OViewer::set_colors(const Eigen::MatrixXd& color_amb, const Eigen::MatrixXd& color_spec, const Eigen::MatrixXd& color_dif)
{
    assert(color_amb.cols() == 3 && color_spec.cols() == 3 && color_dif.cols() == 3 && "Color matrices must have 3 cols.\n");
    assert(color_amb.rows() == color_spec.rows() && color_spec.rows() == color_dif.rows() && "Color matrices must all have the same number of rows.\n");
    assert( (color_amb.rows() == m_V_mod.rows() || color_amb.rows() == m_F_mod.rows() || color_amb.rows() == 3*m_F_mod.rows()) && "Color matrices must either have the same number of rows as vertices, or as faces in the mesh.\n");
    
    m_amb = color_amb.cast<float>();
    m_spec = color_spec.cast<float>();
    m_dif = color_dif.cast<float>();
    
    if(color_amb.rows() == m_F_mod.rows() || color_amb.rows() == 3*m_F_mod.rows())
        faceBasedColoring = true;
    else if(color_amb.rows() == m_V.rows())
        faceBasedColoring = false;
    
    if(colorsBuffersGenerated)
        update_colors();
    
    defaultColorUsed = false;
}

IGL_INLINE void OViewer::set_linetransparency(const Eigen::VectorXd& transp)
{
    m_linetransp = transp.cast<float>();
    
    bool linetranspWasEnabled = linetranspEnabled;
    linetranspEnabled = true;
    if(launched) {
        if(!linetranspWasEnabled)
            update_mesh_data();
        if(linetranspBufferGenerated)
            update_linetransp();
    }
}


IGL_INLINE void OViewer::set_face_based_coloring(const bool& colorMode)
{
    faceBasedColoring = colorMode;
    
    if(colorsBuffersGenerated)
        update_colors();
}

IGL_INLINE bool OViewer::get_face_based_coloring()
{
    return faceBasedColoring;
}

IGL_INLINE const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& OViewer::get_V() const
{
    return m_V_mod;
}

IGL_INLINE void OViewer::modify_V(const Eigen::MatrixXd& V)
{
    m_V_mod = V.cast<float>();
    per_face_normals(m_V_mod, m_F_mod, m_facedVN_mod);
    per_vertex_normals(m_V_mod, m_F_mod, m_facedVN_mod, m_VN_mod);
    update_mesh_data(false);
    
    if(defaultColorUsed) {
        set_colors(defaultColor);
        defaultColorUsed = true;
    }
}


IGL_INLINE void OViewer::modify_F(const Eigen::MatrixXi& F, Eigen::MatrixXi E)
{
    m_F_mod = F.cast<GLuint>();
    per_face_normals(m_V_mod, m_F_mod, m_facedVN_mod);
    per_vertex_normals(m_V_mod, m_F_mod, m_facedVN_mod, m_VN_mod);
    
    if(E.rows()==0)
        igl::edges(F, E);
    m_E_mod = E.cast<GLuint>();
    
    update_mesh_data(true);
    linetranspEnabled = false;
}


IGL_INLINE int OViewer::get_width() const
{
    return windowWidth;
}

IGL_INLINE int OViewer::get_height() const
{
    return windowHeight;
}

IGL_INLINE void OViewer::set_width(const int& w)
{
    assert(w>0 && "Width must be positive.\n");
    windowWidth = w;
}

IGL_INLINE void OViewer::set_height(const int& h)
{
    assert(h>0 && "Height must be positive.\n");
    windowHeight = h;
}

IGL_INLINE void OViewer::set_window_size(const int& w, const int& h)
{
    set_width(w);
    set_height(h);
}

IGL_INLINE OViewer::ProjectionType OViewer::get_projection_type() const
{
    return projectionType;
}

IGL_INLINE void OViewer::set_projection_type(const ProjectionType& type)
{
    projectionType = type;
}

IGL_INLINE OViewer::DisplayMode OViewer::get_display_mode() const
{
    return displayMode;
}

IGL_INLINE void OViewer::set_display_mode(const DisplayMode& mode)
{
    displayMode = mode;
}


IGL_INLINE bool OViewer::get_dragging_enabled() const
{
    return draggingEnabled;
}

IGL_INLINE void OViewer::set_dragging_enabled(const bool& mode)
{
    draggingEnabled = mode;
}


IGL_INLINE Eigen::Vector3f OViewer::get_background_color() const
{
    return backgroundColor;
}

IGL_INLINE void OViewer::set_background_color(const Eigen::Vector3d& color)
{
    backgroundColor = color.cast<float>();
}

IGL_INLINE Eigen::Vector3f OViewer::get_wire_color() const
{
    return wireColor;
}

IGL_INLINE void OViewer::set_wire_color(const Eigen::Vector3d& color)
{
    wireColor = color.cast<float>();
}

IGL_INLINE Eigen::Vector3f OViewer::get_points_outline_color() const
{
    return pointsOutlineColor;
}

IGL_INLINE void OViewer::set_points_outline_color(const Eigen::Vector3d& color)
{
    pointsOutlineColor = color.cast<float>();
}

IGL_INLINE Eigen::Vector3f OViewer::get_points_fill_color() const
{
    return pointsFillColor;
}

IGL_INLINE void OViewer::set_points_fill_color(const Eigen::Vector3d& color)
{
    pointsFillColor = color.cast<float>();
}

IGL_INLINE Eigen::Vector3f OViewer::get_points_important_color() const
{
    return pointsImportantColor;
}

IGL_INLINE void OViewer::set_points_important_color(const Eigen::Vector3d& color)
{
    pointsImportantColor = color.cast<float>();
}

IGL_INLINE Eigen::Vector3f OViewer::get_points_selected_color() const
{
    return pointsSelectedColor;
}

IGL_INLINE void OViewer::set_points_selected_color(const Eigen::Vector3d& color)
{
    pointsSelectedColor = color.cast<float>();
}


IGL_INLINE void OViewer::set_point_size(const float& size)
{
    set_point_size(size, size, 4./3.*size);
}

IGL_INLINE void OViewer::set_point_size(const float& size, const float& importantSize, const float& selectedSize)
{
    pointSize = size;
    importantPointSize = importantSize;
    selectedPointSize = selectedSize;
}


IGL_INLINE float OViewer::get_camera_zoom() const
{
    return cameraZoom;
}

IGL_INLINE void OViewer::set_camera_zoom(const float& zoom)
{
    cameraZoom = zoom;
}

IGL_INLINE Eigen::Quaternion<float> OViewer::get_camera_rotation() const
{
    return cameraRotation;
}

IGL_INLINE void OViewer::set_camera_rotation(const Eigen::Quaternion<float>& rotation)
{
    cameraRotation = rotation;
}

IGL_INLINE Eigen::Vector3f OViewer::get_camera_translation() const
{
    return cameraTranslation;
}

IGL_INLINE void OViewer::set_camera_translation(const Eigen::Vector3f& translation)
{
    cameraTranslation = translation;
}


IGL_INLINE bool OViewer::get_rotation_enabled() const
{
    return rotationEnabled;
}

IGL_INLINE void OViewer::set_rotation_enabled(const bool& enabled)
{
    rotationEnabled = enabled;
}


IGL_INLINE Eigen::Vector3f OViewer::get_light_position() const
{
    return lightPosition;
}

IGL_INLINE void OViewer::set_light_position(const Eigen::Vector3f& position)
{
    lightPosition = position;
    if(colorsBuffersGenerated)
        update_lighting();
}

IGL_INLINE float OViewer::get_lighting_factor() const
{
    return lightingFactor;
}

IGL_INLINE void OViewer::set_lighting_factor(const float& factor)
{
    lightingFactor = factor;
    if(colorsBuffersGenerated)
        update_lighting();
}

IGL_INLINE float OViewer::get_shininess() const
{
    return shininess;
}

IGL_INLINE void OViewer::set_shininess(const float& shiny)
{
    shininess = shiny;
    if(colorsBuffersGenerated)
        update_lighting();
}

IGL_INLINE bool OViewer::get_cel_shading() const
{
    return celShading;
}

IGL_INLINE void OViewer::set_cel_shading(const bool& cel)
{
    celShading = cel;
    if(colorsBuffersGenerated)
        update_lighting();
}


IGL_INLINE void OViewer::mark_points(const Eigen::VectorXi& points)
{
    markedPoints = points;
    pointIsImportant = Eigen::VectorXi::Constant(points.size(), false);
    pointIsSelected = Eigen::VectorXi::Constant(points.size(), false);
    markedPointsTranslations = Eigen::MatrixXd::Constant(points.size(), 3, 0.);
}

IGL_INLINE const Eigen::VectorXi& OViewer::get_marked_points() const
{
    return markedPoints;
}

IGL_INLINE void OViewer::set_important_points(const Eigen::VectorXi& points)
{
    pointIsImportant = points;
}

IGL_INLINE const Eigen::VectorXi& OViewer::get_important_points() const
{
    return pointIsImportant;
}

IGL_INLINE const Eigen::VectorXi& OViewer::get_selected_points() const
{
    return pointIsSelected;
}

IGL_INLINE const Eigen::MatrixXd& OViewer::get_marked_points_translations() const
{
    return markedPointsTranslations;
}


IGL_INLINE const Eigen::Matrix4f& OViewer::get_model() const
{
    return model;
}

IGL_INLINE const Eigen::Matrix4f& OViewer::get_view() const
{
    return view;
}

IGL_INLINE const Eigen::Matrix4f& OViewer::get_proj() const
{
    return proj;
}


IGL_INLINE void OViewer::stop_dragging()
{
    currentlyDragging = false;
}

IGL_INLINE bool OViewer::get_currently_dragging() const
{
    return currentlyDragging;
}


IGL_INLINE bool OViewer::showLines() const
{
    return displayMode == Wireframe || displayMode == FillWireframe || displayMode == TextureWireframe;
}

IGL_INLINE bool OViewer::showFill() const
{
    return displayMode == Fill || displayMode == FillWireframe;
}

IGL_INLINE bool OViewer::showTexture() const
{
    return displayMode == Texture || displayMode == TextureWireframe;
}


IGL_INLINE void OViewer::set_plottitles(const std::vector<std::string>& titles)
{
    plotTitles = titles;
    if(plotBuffersGenerated)
        update_plots();
}

IGL_INLINE void OViewer::set_plotlines(const std::vector<Eigen::VectorXd>& lines)
{
    plotLines = lines;
    if(plotBuffersGenerated)
        update_plots();
}

IGL_INLINE void OViewer::set_plots(const std::vector<std::string>& titles, const std::vector<Eigen::VectorXd>& lines)
{
    plotTitles = titles;
    plotLines = lines;
    if(plotBuffersGenerated)
        update_plots();
}

IGL_INLINE void OViewer::set_plotcolors(const std::vector<Eigen::Vector3d>& colors)
{
    plotColors = std::vector<Eigen::Vector3f>(colors.size());
    for(int i=0; i<plotColors.size(); ++i) {
        plotColors[i] = colors[i].cast<float>();
    }
    if(plotBuffersGenerated)
        update_plots();
}

IGL_INLINE void OViewer::set_plotting_enabled(const bool& mode)
{
    plottingEnabled = mode;
    if(plotBuffersGenerated)
        update_plots();
}

IGL_INLINE bool OViewer::get_plotting_enabled()
{
    return plottingEnabled;
}

IGL_INLINE void OViewer::set_plotting_logplot(const bool& mode)
{
    plottingLogPlot = mode;
    if(plotBuffersGenerated)
        update_plots();
}

IGL_INLINE bool OViewer::get_plotting_logplot()
{
    return plottingLogPlot;
}


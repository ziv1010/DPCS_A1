// Renderer.h
#ifndef RENDERER_H
#define RENDERER_H

#include "Model_3D.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>
#include <vector>

class Renderer
{
public:
    Renderer();
    ~Renderer();

    // Render the model based on model type (2D or 3D)
    void renderModel(const Model_3D& model, int modelType);

    // Save the rendered image to a file
    void saveImage(const std::string& filename);

private:
    // Initialization and cleanup functions
    void initialize();
    void cleanup();

    // Setup framebuffer for offscreen rendering (for 2D views)
    void setupFramebuffer();

    // Render all views for a 3D model (for 3D to 2D case)
    void renderViews(const Model_3D& model);

    // Render a single 2D view
    void renderView(const Model_2D& view);

    // **New methods for interactive 3D rendering**
    void renderInteractive(const Model_3D& model);
    void drawModel();

    // OpenGL and rendering parameters
    GLuint fbo;              // Framebuffer Object
    GLuint texture;          // Texture attached to FBO
    int total_width;         // Total image width
    int total_height;        // Total image height
    int view_width;          // Width of each view
    int view_height;         // Height of each view
    GLFWwindow* window;      // GLFW window for context
    std::vector<unsigned char> pixels; // Pixel data for the rendered image
    bool initialized;        // Flag to check if OpenGL is initialized

    // **Camera and interaction parameters**
    float camera_distance;
    float camera_angle_x;
    float camera_angle_y;
    bool left_mouse_button_pressed;
    double last_mouse_x;
    double last_mouse_y;

    // **Model reference**
    const Model_3D* model;

    // **Input handling**
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
    static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

 // **New methods and variables**
    void calculateModelBounds(const Model_3D& model, float& centerX, float& centerY, float& centerZ, float& maxExtent);
    void drawAxes();
    void render2DView(const Model_2D& view);


// Model center and extent
    float model_center_x;
    float model_center_y;
    float model_center_z;
    float model_max_extent;


    // **Static instance pointer for callbacks**
    static Renderer* instance;
};

#endif // RENDERER_H

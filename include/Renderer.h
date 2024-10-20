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

    // Setup framebuffer for offscreen rendering
    void setupFramebuffer();

    // Render all views for a 3D model
    void renderViews(const Model_3D& model);

    // Render a single 2D view
    void renderView(const Model_2D& view);

    // OpenGL and rendering parameters
    GLuint fbo;              // Framebuffer Object
    GLuint texture;          // Texture attached to FBO
    int total_width;         // Total image width
    int total_height;        // Total image height
    int view_width;          // Width of each view
    int view_height;         // Height of each view
    GLFWwindow* window;      // Hidden GLFW window for context
    std::vector<unsigned char> pixels; // Pixel data for the rendered image
    bool initialized;        // Flag to check if OpenGL is initialized
};

#endif // RENDERER_H

#include "Renderer.h"
#include <iostream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


Renderer::Renderer()
{
    fbo = 0;
    texture = 0;
    total_width = 0;
    total_height = 0;
    view_width = 800;
    view_height = 800;
    window = nullptr;
    initialized = false;
}

Renderer::~Renderer()
{
    cleanup();
}

void Renderer::initialize()
{
    if (!initialized)
    {
        // Initialize GLFW
        if (!glfwInit())
        {
            std::cerr << "Failed to initialize GLFW\n";
            return;
        }

        // Create hidden window
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        window = glfwCreateWindow(800, 800, "Offscreen", NULL, NULL);
        if (!window)
        {
            std::cerr << "Failed to create hidden window\n";
            glfwTerminate();
            return;
        }

        glfwMakeContextCurrent(window);

        // Initialize GLEW
        glewExperimental = GL_TRUE; // Ensure GLEW uses modern techniques for managing OpenGL functionality
        if (glewInit() != GLEW_OK)
        {
            std::cerr << "Failed to initialize GLEW\n";
            glfwDestroyWindow(window);
            glfwTerminate();
            return;
        }

        initialized = true;
    }
}

void Renderer::cleanup()
{
    if (texture != 0)
    {
        glDeleteTextures(1, &texture);
        texture = 0;
    }
    if (fbo != 0)
    {
        glDeleteFramebuffers(1, &fbo);
        fbo = 0;
    }
    if (window != nullptr)
    {
        glfwDestroyWindow(window);
        window = nullptr;
    }
    glfwTerminate();
    initialized = false;
}

void Renderer::setupFramebuffer()
{
    // Calculate total dimensions for three views side by side
    total_width = view_width * 3;
    total_height = view_height;

    // Create Framebuffer Object (FBO)
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    // Create Texture to render to
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, total_width, total_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Attach texture to FBO
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

    // Check FBO status
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        std::cerr << "Failed to create framebuffer\n";
        cleanup();
        return;
    }
}

void Renderer::renderModel(const Model_3D& model, int modelType)
{
    initialize();
    if (!initialized)
    {
        return;
    }

    setupFramebuffer();

    // Set the viewport for the entire framebuffer
    glViewport(0, 0, total_width, total_height);

    // Clear the framebuffer
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // White background
    glClear(GL_COLOR_BUFFER_BIT);

    // Determine which model to draw based on modelType
    if (modelType == 3)
    {
        // Draw all 2D projections
        renderViews(model);
    }
    else if (modelType == 2)
    {
        // Handle 3D model rendering (not implemented in this example)
        std::cerr << "3D model rendering is not implemented in Renderer.\n";
    }
    else
    {
        std::cerr << "Invalid modelType provided to Renderer.\n";
    }

    // Read pixels from framebuffer
    pixels.resize(total_width * total_height * 3);
    glReadPixels(0, 0, total_width, total_height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

    // Unbind the framebuffer
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Renderer::renderViews(const Model_3D& model)
{
    Model_2D views[3] = { model.topView, model.frontView, model.sideView };

    for (int i = 0; i < 3; ++i)
    {
        // Set the viewport for the current view
        glViewport(i * view_width, 0, view_width, view_height);

        // Set up projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)view_width / (float)view_height;
        glOrtho(-1.5f * aspect, 1.5f * aspect, -1.5f, 1.5f, -1.0f, 1.0f);

        // Set up modelview matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Render the current view
        renderView(views[i]);
    }
}

void Renderer::renderView(const Model_2D& view)
{
    // Find max coordinate for scaling (same as before)
    float max_coord = 1.0f;
    for (const auto& vertex : view.v)
    {
        if (view.direction == 0) // Top View
        {
            max_coord = std::max(max_coord, std::abs(vertex.x));
            max_coord = std::max(max_coord, std::abs(vertex.y));
        }
        else if (view.direction == 1) // Front View
        {
            max_coord = std::max(max_coord, std::abs(vertex.x));
            max_coord = std::max(max_coord, std::abs(vertex.z));
        }
        else if (view.direction == 2) // Side View
        {
            max_coord = std::max(max_coord, std::abs(vertex.y));
            max_coord = std::max(max_coord, std::abs(vertex.z));
        }
    }
    max_coord *= 1.1f; // Add a margin

    // Set the line stipple pattern once
    glLineStipple(1, 0x00FF); // Dashed line pattern

    // Draw edges
    glColor3f(0.0f, 0.0f, 0.0f); // Black color for lines
    for (const auto& edge : view.e)
    {
        if (!edge.isTrue)
            continue;

        // Enable or disable line stippling based on edge.hidden
        if (edge.hidden)
        {
            glEnable(GL_LINE_STIPPLE);
        }
        else
        {
            glDisable(GL_LINE_STIPPLE);
        }

        glBegin(GL_LINES);
        float ax = 0.0f, ay = 0.0f, bx = 0.0f, by = 0.0f;

        // Adjust coordinate selection based on the view's direction
        if (view.direction == 0) // Top View
        {
            ax = edge.a.x / max_coord;
            ay = edge.a.y / max_coord;
            bx = edge.b.x / max_coord;
            by = edge.b.y / max_coord;
        }
        else if (view.direction == 1) // Front View
        {
            ax = edge.a.x / max_coord;
            ay = edge.a.z / max_coord;
            bx = edge.b.x / max_coord;
            by = edge.b.z / max_coord;
        }
        else if (view.direction == 2) // Side View
        {
            ax = edge.a.y / max_coord;
            ay = edge.a.z / max_coord;
            bx = edge.b.y / max_coord;
            by = edge.b.z / max_coord;
        }
        else
        {
            std::cerr << "Error: Invalid view direction: " << view.direction << std::endl;
            continue;
        }

        glVertex2f(ax, ay);
        glVertex2f(bx, by);
        glEnd();
    }

    // Disable line stippling after drawing all edges
    glDisable(GL_LINE_STIPPLE);
}


void Renderer::saveImage(const std::string& filename)
{
    if (pixels.empty())
    {
        std::cerr << "No image data to save. Did you call renderModel()?\n";
        return;
    }

    // OpenGL images are bottom-up, so we need to flip the image vertically
    std::vector<unsigned char> flipped_pixels(total_width * total_height * 3);
    for (int y = 0; y < total_height; ++y)
    {
        memcpy(&flipped_pixels[y * total_width * 3], &pixels[(total_height - 1 - y) * total_width * 3], total_width * 3);
    }

    // Save the image using stb_image_write
    if (stbi_write_png(filename.c_str(), total_width, total_height, 3, flipped_pixels.data(), total_width * 3))
    {
        std::cout << "Image saved to " << filename << "\n";
    }
    else
    {
        std::cerr << "Failed to save image\n";
    }
}

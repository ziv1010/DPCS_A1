// Renderer.cpp
#include "Renderer.h"
#include <iostream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cfloat> // For FLT_MAX
#include <algorithm> // For std::min and std::max

// Initialize the static instance pointer
Renderer* Renderer::instance = nullptr;

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

    // Initialize camera and interaction parameters
    camera_distance = 5.0f;
    camera_angle_x = 0.0f;
    camera_angle_y = 0.0f;
    left_mouse_button_pressed = false;
    last_mouse_x = 0.0;
    last_mouse_y = 0.0;
    model = nullptr;
}

Renderer::~Renderer()
{
    cleanup();
}

void Renderer::calculateModelBounds(const Model_3D& model, float& centerX, float& centerY, float& centerZ, float& maxExtent)
{
    float minX = FLT_MAX, minY = FLT_MAX, minZ = FLT_MAX;
    float maxX = -FLT_MAX, maxY = -FLT_MAX, maxZ = -FLT_MAX;

    for (const auto& vertex : model.v)
    {
        if (!vertex.isTrue)
            continue;

        minX = std::min(minX, vertex.x);
        minY = std::min(minY, vertex.y);
        minZ = std::min(minZ, vertex.z);

        maxX = std::max(maxX, vertex.x);
        maxY = std::max(maxY, vertex.y);
        maxZ = std::max(maxZ, vertex.z);
    }

    centerX = (minX + maxX) / 2.0f;
    centerY = (minY + maxY) / 2.0f;
    centerZ = (minZ + maxZ) / 2.0f;

    float extentX = maxX - minX;
    float extentY = maxY - minY;
    float extentZ = maxZ - minZ;

    maxExtent = std::max({ extentX, extentY, extentZ });
}


void Renderer::initialize()
{

    if (!initialized)
    {
        // Set the instance pointer for callbacks
        instance = this;

        // Initialize GLFW
        if (!glfwInit())
        {
            std::cerr << "Failed to initialize GLFW\n";
            return;
        }

        // Create window (visible for 3D rendering)
        window = glfwCreateWindow(view_width, view_height, "3D Model Viewer", NULL, NULL);
        if (!window)
        {
            std::cerr << "Failed to create GLFW window\n";
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

        // Set callbacks
        glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
        glfwSetCursorPosCallback(window, cursor_position_callback);
        glfwSetMouseButtonCallback(window, mouse_button_callback);
        glfwSetScrollCallback(window, scroll_callback);
        glfwSetKeyCallback(window, key_callback);
        // Enable depth testing
        glEnable(GL_DEPTH_TEST);

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

void Renderer::renderModel(const Model_3D& model, int modelType)
{
    initialize();
    if (!initialized)
    {
        return;
    }

    if (modelType == 3)
    {
        // Set up framebuffer for offscreen rendering
        setupFramebuffer();

        // Set the viewport for the entire framebuffer
        glViewport(0, 0, total_width, total_height);

        // Clear the framebuffer
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // White background
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw all 2D projections
        renderViews(model);

        // Read pixels from framebuffer
        pixels.resize(total_width * total_height * 3);
        glReadPixels(0, 0, total_width, total_height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

        // Unbind the framebuffer
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
    else if (modelType == 2)
    {
        // Render the 3D model interactively
        renderInteractive(model);
    }
    else
    {
        std::cerr << "Invalid modelType provided to Renderer.\n";
    }
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
void Renderer::renderInteractive(const Model_3D& model)
{
    this->model = &model;

    // Calculate model bounds
    float centerX, centerY, centerZ, maxExtent;
    calculateModelBounds(model, centerX, centerY, centerZ, maxExtent);

    // Store model center and extent for use in other methods
    model_center_x = centerX;
    model_center_y = centerY;
    model_center_z = centerZ;
    model_max_extent = maxExtent;

    // Adjust camera distance based on model size
    camera_distance = maxExtent * 1.5f;

    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        // Handle events
        glfwPollEvents();

        // Clear buffers
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f); // Light gray background
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Set up the projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect_ratio = static_cast<float>(view_width) / static_cast<float>(view_height);
        gluPerspective(45.0, aspect_ratio, 0.1, model_max_extent * 20.0f);

        // Set up the modelview matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Apply camera transformations
        glTranslatef(0.0f, 0.0f, -camera_distance);
        glRotatef(camera_angle_x, 1.0f, 0.0f, 0.0f);
        glRotatef(camera_angle_y, 0.0f, 1.0f, 0.0f);

        // Translate model to center it
        glTranslatef(-model_center_x, -model_center_y, -model_center_z);

        // Draw the model
        drawModel();

        // Swap buffers
        glfwSwapBuffers(window);
    }

    // Clean up after rendering
    cleanup();
}


void Renderer::drawModel()
{
    if (!model)
        return;

    // Draw coordinate axes
    drawAxes();

    // Draw surfaces (if available)
    glColor3f(0.8f, 0.8f, 0.8f); // Light gray color for surfaces
    for (const auto& surface : model->s)
    {
        glBegin(GL_POLYGON);
        for (const auto& edge : surface.boundary)
        {
            glVertex3f(edge.a.x, edge.a.y, edge.a.z);
        }
        glEnd();
    }

    // Draw edges
    glColor3f(0.0f, 0.0f, 0.0f); // Black color for edges
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    for (const auto& edge : model->e)
    {
        if (!edge.isTrue)
            continue;

        glVertex3f(edge.a.x, edge.a.y, edge.a.z);
        glVertex3f(edge.b.x, edge.b.y, edge.b.z);
    }
    glEnd();
}


// Callback functions
void Renderer::framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    if (instance)
    {
        instance->view_width = width;
        instance->view_height = height;
        glViewport(0, 0, width, height);
    }
}

void Renderer::cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (instance && instance->left_mouse_button_pressed)
    {
        float dx = static_cast<float>(xpos - instance->last_mouse_x);
        float dy = static_cast<float>(ypos - instance->last_mouse_y);

        instance->camera_angle_y += dx * 0.5f;
        instance->camera_angle_x += dy * 0.5f;

        // Limit the vertical rotation angle to prevent flipping
        if (instance->camera_angle_x > 89.0f)
            instance->camera_angle_x = 89.0f;
        if (instance->camera_angle_x < -89.0f)
            instance->camera_angle_x = -89.0f;
    }

    if (instance)
    {
        instance->last_mouse_x = xpos;
        instance->last_mouse_y = ypos;
    }
}

void Renderer::mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (instance)
    {
        if (button == GLFW_MOUSE_BUTTON_LEFT)
        {
            if (action == GLFW_PRESS)
            {
                instance->left_mouse_button_pressed = true;
                glfwGetCursorPos(window, &instance->last_mouse_x, &instance->last_mouse_y);
            }
            else if (action == GLFW_RELEASE)
            {
                instance->left_mouse_button_pressed = false;
            }
        }
    }
}

void Renderer::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    if (instance)
    {
        float zoomSensitivity = 0.1f;
        instance->camera_distance -= static_cast<float>(yoffset) * zoomSensitivity * instance->model_max_extent;

        float minDistance = instance->model_max_extent * 0.1f; // Minimum zoom in
        float maxDistance = instance->model_max_extent * 10.0f; // Maximum zoom out

        if (instance->camera_distance < minDistance)
            instance->camera_distance = minDistance;
        if (instance->camera_distance > maxDistance)
            instance->camera_distance = maxDistance;
    }
}

void Renderer::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (instance)
    {
        if (action == GLFW_PRESS || action == GLFW_REPEAT)
        {
            float zoomStep = instance->model_max_extent * 0.05f;

            if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_KP_ADD) // '+' key
            {
                instance->camera_distance -= zoomStep;
                if (instance->camera_distance < instance->model_max_extent * 0.1f)
                    instance->camera_distance = instance->model_max_extent * 0.1f;
            }
            else if (key == GLFW_KEY_MINUS || key == GLFW_KEY_KP_SUBTRACT) // '-' key
            {
                instance->camera_distance += zoomStep;
                if (instance->camera_distance > instance->model_max_extent * 10.0f)
                    instance->camera_distance = instance->model_max_extent * 10.0f;
            }
        }
    }
}

void Renderer::drawAxes()
{
    float axisLength = model_max_extent * 1.2f;

    glLineWidth(2.0f);

    // X-axis
    glColor3f(1.0f, 0.0f, 0.0f); // Red
    glBegin(GL_LINES);
    glVertex3f(-axisLength, 0.0f, 0.0f);
    glVertex3f(axisLength, 0.0f, 0.0f);
    glEnd();

    // Y-axis
    glColor3f(0.0f, 1.0f, 0.0f); // Green
    glBegin(GL_LINES);
    glVertex3f(0.0f, -axisLength, 0.0f);
    glVertex3f(0.0f, axisLength, 0.0f);
    glEnd();

    // Z-axis
    glColor3f(0.0f, 0.0f, 1.0f); // Blue
    glBegin(GL_LINES);
    glVertex3f(0.0f, 0.0f, -axisLength);
    glVertex3f(0.0f, 0.0f, axisLength);
    glEnd();
}




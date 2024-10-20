// main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Model_3D.h"
#include "Model_2D.h"
#include "Renderer.h" // Include the Renderer header

// Global Variables

// Model Instances
Model_3D Model;
Model_2D Top;
Model_2D Front;
Model_2D Side;

// File Paths
std::string file_path;
std::string save_file_path;

// Model Type
int modelType = 0;  // 2 for 2D input, 3 for 3D input

// Function Declarations
void generate_model();
void save_output();

int main(int argc, char** argv)
{
    // Input File Path
    if (argc >= 2)
    {
        file_path = argv[1];
    }
    else
    {
        std::cout << "Enter the input file path: ";
        std::getline(std::cin, file_path);
    }

    // Generate Model from Input File
    generate_model();
    if (modelType == 0)
    {
        std::cerr << "Model type not set. Exiting.\n";
        return -1;
    }

    // Save Output
    save_output();

    // Render and save image using Renderer
    Renderer renderer;
    renderer.renderModel(Model, modelType);

    if (modelType == 3)
    {
        std::string output_image_path = "output.png"; // Or prompt the user for the output image path
        renderer.saveImage(output_image_path);
    }

    return 0;
}

void generate_model()
{
    // Read the input file
    std::ifstream inFile(file_path);
    if (!inFile.is_open())
    {
        std::cerr << "Failed to open input file: " << file_path << "\n";
        return;
    }

    int n;
    inFile >> n;
    modelType = n;  // Set the global modelType variable

    if (n == 2)
    {
        // File is 2D
        std::cout << "File is 2D\n";

        // Read Top View Details
        int dir;
        inFile >> dir;
        if (dir != 0)
        {
            std::cerr << "Wrong input for direction in Top View!\n";
            return;
        }
        Top.direction = dir;

        int nv;
        inFile >> nv;
        for (int i = 0; i < nv; ++i)
        {
            int a;
            float b, c;
            inFile >> a >> b >> c;
            Vertex t(b, c, 0.0f, a);
            Top.v.push_back(t);
        }

        int ne;
        inFile >> ne;
        for (int i = 0; i < ne; ++i)
        {
            int a, b, c, d;
            inFile >> a >> b >> c >> d;
            Edge t(Top.v[b - 1], Top.v[c - 1], a);
            t.hidden = (d == 1);
            Top.e.push_back(t);
        }

        // Read Front View Details
        inFile >> dir;
        if (dir != 1)
        {
            std::cerr << "Wrong input for direction in Front View!\n";
            return;
        }
        Front.direction = dir;

        inFile >> nv;
        for (int i = 0; i < nv; ++i)
        {
            int a;
            float b, c;
            inFile >> a >> b >> c;
            Vertex t(b, 0.0f, c, a);
            Front.v.push_back(t);
        }

        inFile >> ne;
        for (int i = 0; i < ne; ++i)
        {
            int a, b, c, d;
            inFile >> a >> b >> c >> d;
            Edge t(Front.v[b - 1], Front.v[c - 1], a);
            t.hidden = (d == 1);
            Front.e.push_back(t);
        }

        // Read Side View Details
        inFile >> dir;
        if (dir != 2)
        {
            std::cerr << "Wrong input for direction in Side View!\n";
            return;
        }
        Side.direction = dir;

        inFile >> nv;
        for (int i = 0; i < nv; ++i)
        {
            int a;
            float b, c;
            inFile >> a >> b >> c;
            Vertex t(0.0f, b, c, a);
            Side.v.push_back(t);
        }

        inFile >> ne;
        for (int i = 0; i < ne; ++i)
        {
            int a, b, c, d;
            inFile >> a >> b >> c >> d;
            Edge t(Side.v[b - 1], Side.v[c - 1], a);
            t.hidden = (d == 1);
            Side.e.push_back(t);
        }

        // Assign Views to Model
        Model.frontView = Front;
        Model.topView = Top;
        Model.sideView = Side;

        // Form Wireframe from 2D Views
        Model.form_wireframe();

        // Add this line to form 3D faces
        Model.form_3dfaces();
    }
    else if (n == 3)
    {
        // File is 3D
        std::cout << "File is 3D\n";

        // Read 3D Model Details
        int nv;
        inFile >> nv;
        for (int i = 0; i < nv; ++i)
        {
            int a;
            float b, c, d;
            inFile >> a >> b >> c >> d;
            Vertex t(b, c, d, a);
            Model.v.push_back(t);
        }

        int ne;
        inFile >> ne;
        for (int i = 0; i < ne; ++i)
        {
            int a, b, c;
            inFile >> a >> b >> c;
            Edge t(Model.v[b - 1], Model.v[c - 1], a);
            Model.e.push_back(t);
        }

        int ns;
        inFile >> ns;
        for (int i = 0; i < ns; ++i)
        {
            int c;
            inFile >> c;
            int a;
            inFile >> a;
            std::vector<Edge> t;
            for (int j = 0; j < a; ++j)
            {
                int b;
                inFile >> b;
                t.push_back(Model.e[b - 1]);
            }
            Surface t1(t, c);
            Model.s.push_back(t1);
        }

        // Generate 2D Projections from 3D Model
        Model.topView = Model.basic_proj(0);
        Model.frontView = Model.basic_proj(1);
        Model.sideView = Model.basic_proj(2);

        Model.calc_2d(0);
        Model.calc_2d(1);
        Model.calc_2d(2);
    }
    else
    {
        std::cerr << "Invalid file format!\n";
        return;
    }

    inFile.close();
}

void save_output()
{
    std::cout << "Enter the output file path: ";
    std::getline(std::cin, save_file_path);

    std::ofstream outFile(save_file_path);
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open output file: " << save_file_path << "\n";
        return;
    }

    if (modelType == 3)
    {
        Model.save2Dmodels(outFile);
    }
    else if (modelType == 2)
    {
        Model.save3Dmodel(outFile);
    }
    else
    {
        std::cerr << "Invalid model type!\n";
    }

    outFile.close();
}



/*void render_and_save_image()
{
    // Initialize GLFW
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW\n";
        return;
    }

    // Create hidden window
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    GLFWwindow* window = glfwCreateWindow(800, 800, "Offscreen", NULL, NULL);
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

    // Create Framebuffer Object (FBO)
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    // Create Texture to render to
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    const int width = 800;
    const int height = 800;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Attach texture to FBO
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

    // Check FBO status
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        std::cerr << "Failed to create framebuffer\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return;
    }

    // Set the viewport
    glViewport(0, 0, width, height);

    // Set up projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = (float)width / (float)height;
    glOrtho(-1.5 * aspect, 1.5 * aspect, -1.5, 1.5, -1.0, 1.0);

    // Set up modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Clear the framebuffer
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // White background
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw the model
    glColor3f(0.0f, 0.0f, 0.0f); // Black color for lines

    // Determine which model to draw based on modelType
    if (modelType == 3)
    {
        // Draw 2D projections
        // Example: Draw top view
        Model_2D* view = &Model.topView;

        // Find max coordinate for scaling
        float max_coord = 1.0f;
        for (const auto& vertex : view->v)
        {
            max_coord = std::max(max_coord, std::abs(vertex.x));
            max_coord = std::max(max_coord, std::abs(vertex.y));
        }
        max_coord *= 1.1f; // Add a margin

        // Draw edges
        for (const auto& edge : view->e)
        {
            if (!edge.isTrue)
                continue;

            if (edge.hidden)
            {
                glLineStipple(1, 0x00FF);
                glEnable(GL_LINE_STIPPLE);
            }
            else
            {
                glDisable(GL_LINE_STIPPLE);
            }

            glBegin(GL_LINES);
            glVertex2f(edge.a.x / max_coord, edge.a.y / max_coord);
            glVertex2f(edge.b.x / max_coord, edge.b.y / max_coord);
            glEnd();
        }

        glDisable(GL_LINE_STIPPLE);
    }
    else if (modelType == 2)
    {
        // Draw 3D model
        // Project the 3D model onto 2D plane for rendering
        // Simple orthographic projection onto XY plane

        // Find max coordinate for scaling
        float max_coord = 1.0f;
        for (const auto& vertex : Model.v)
        {
            if (!vertex.isTrue)
                continue;
            max_coord = std::max(max_coord, std::abs(vertex.x));
            max_coord = std::max(max_coord, std::abs(vertex.y));
        }
        max_coord *= 1.1f; // Add a margin

        // Draw edges
        for (const auto& edge : Model.e)
        {
            if (!edge.isTrue)
                continue;

            if (edge.hidden)
            {
                glLineStipple(1, 0x00FF);
                glEnable(GL_LINE_STIPPLE);
            }
            else
            {
                glDisable(GL_LINE_STIPPLE);
            }

            glBegin(GL_LINES);
            glVertex2f(edge.a.x / max_coord, edge.a.y / max_coord);
            glVertex2f(edge.b.x / max_coord, edge.b.y / max_coord);
            glEnd();
        }

        glDisable(GL_LINE_STIPPLE);
    }

    // Read pixels from framebuffer
    std::vector<unsigned char> pixels(width * height * 3);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

    // OpenGL images are bottom-up, so we need to flip the image vertically
    std::vector<unsigned char> flipped_pixels(width * height * 3);
    for (int y = 0; y < height; ++y)
    {
        memcpy(&flipped_pixels[y * width * 3], &pixels[(height - 1 - y) * width * 3], width * 3);
    }

    // Save the image using stb_image_write
    const char* output_image = "output.png";
    if (stbi_write_png(output_image, width, height, 3, flipped_pixels.data(), width * 3))
    {
        std::cout << "Image saved to " << output_image << "\n";
    }
    else
    {
        std::cerr << "Failed to save image\n";
    }

    // Clean up
    glDeleteTextures(1, &texture);
    glDeleteFramebuffers(1, &fbo);

    glfwDestroyWindow(window);
    glfwTerminate();
}

*/


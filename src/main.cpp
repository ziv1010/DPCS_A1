#include <iostream>
#include "Object3D.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>

int main() {
    try {
        Object3D obj;
        obj.loadFromFile("cube.txt");

        std::cout << "Cube created with " << obj.getFace(0).id + 1 << " face(s)." << std::endl;

        // Initialize GLFW
        if (!glfwInit()) {
            std::cerr << "Failed to initialize GLFW" << std::endl;
            return -1;
        }

        // Create a windowed mode window and its OpenGL context
        GLFWwindow* window = glfwCreateWindow(800, 600, "Object3D Viewer", NULL, NULL);
        if (!window) {
            glfwTerminate();
            std::cerr << "Failed to create GLFW window" << std::endl;
            return -1;
        }

        // Make the window's context current
        glfwMakeContextCurrent(window);

        // Initialize GLEW
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) {
            std::cerr << "Failed to initialize GLEW" << std::endl;
            return -1;
        }

        // Initialize OpenGL buffers
        obj.initializeOpenGL();

        // Main loop
        while (!glfwWindowShouldClose(window)) {
            // Render here
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Draw the object
            obj.draw();

            // Swap front and back buffers
            glfwSwapBuffers(window);

            // Poll for and process events
            glfwPollEvents();
        }

        glfwDestroyWindow(window);
        glfwTerminate();
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return -1;
    }

    return 0;
}
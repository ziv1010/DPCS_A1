#include <iostream>
#include "Object3D.h"
#include "Transformation.h"  // Include Transformation header
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

        // Set up a larger orthographic projection to view the transformations better
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-2.0, 2.0, -2.0, 2.0, -1.0, 1.0);  // Modify the range to zoom out a bit
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Initialize OpenGL buffers
        obj.initializeOpenGL();

        // Main loop
        while (!glfwWindowShouldClose(window)) {
            // Clear the screen
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Apply transformations
            glPushMatrix();  // Save the current matrix state
            glTranslatef(0.5f, 0.5f, 0.0f);  // Translate by 0.5 units on X and Y
            glScalef(1.5f, 1.5f, 1.5f);      // Scale by 1.5x
            glRotatef(45.0f, 0.0f, 0.0f, 1.0f);  // Rotate by 45 degrees around the Z-axis

            // Draw the object
            obj.draw();

            glPopMatrix();  // Restore the original matrix state

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
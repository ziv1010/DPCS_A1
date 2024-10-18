#include "OrthographicView.h"
#include "Object3D.h"

void OrthographicView::projectOntoXY(const Object3D& obj) {
    // Use the getter method to access vertices
    const auto& vertices = obj.getVertices();

    // Render 2D projection of obj on the XY plane using OpenGL
    glBegin(GL_LINES);
    for (const auto& [id, vertex] : vertices) {
        glVertex2f(vertex.x, vertex.y);  // Discard the z-coordinate
    }
    glEnd();
}

void OrthographicView::projectOntoXZ(const Object3D& obj) {
    const auto& vertices = obj.getVertices();

    glBegin(GL_LINES);
    for (const auto& [id, vertex] : vertices) {
        glVertex2f(vertex.x, vertex.z);  // Discard the y-coordinate
    }
    glEnd();
}

void OrthographicView::projectOntoYZ(const Object3D& obj) {
    const auto& vertices = obj.getVertices();

    glBegin(GL_LINES);
    for (const auto& [id, vertex] : vertices) {
        glVertex2f(vertex.y, vertex.z);  // Discard the x-coordinate
    }
    glEnd();
}
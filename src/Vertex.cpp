// Vertex.cpp

#include "Vertex.h"
#include <cmath>    // For sqrt and pow
#include <glm/glm.hpp>  // For GLM vector operations

#include "Vertex.h"

// Default constructor
Vertex::Vertex() : x(0.0f), y(0.0f), z(0.0f), vNo(0), isTrue(true) {}

// Parameterized constructor
Vertex::Vertex(float x_, float y_, float z_, int n_) : x(x_), y(y_), z(z_), vNo(n_), isTrue(true) {}

// Copy constructor
Vertex::Vertex(const Vertex& v) : x(v.x), y(v.y), z(v.z), vNo(v.vNo), isTrue(v.isTrue) {}

// Get the XY Projection of the Vertex for Top View
Vertex Vertex::getXY() const {
    // Generates a new Vertex object by projecting the current vertex onto the XY-plane.
    // The Z-coordinate is eliminated (set to zero) to represent the top view.
    return Vertex(this->x, this->y, 0.0f, this->vNo);
}

// Get the YZ Projection of the Vertex for Side View
Vertex Vertex::getYZ() const {
    // Generates a new Vertex object by projecting the current vertex onto the YZ-plane.
    // The X-coordinate is eliminated (set to zero) to represent the side view.
    return Vertex(0.0f, this->y, this->z, this->vNo);
}

// Get the XZ Projection of the Vertex for Right Side View
Vertex Vertex::getXZ() const {
    // Generates a new Vertex object by projecting the current vertex onto the XZ-plane.
    // The Y-coordinate is eliminated (set to zero) to represent the right side view.
    return Vertex(this->x, 0.0f, this->z, this->vNo);
}

// Compare if Two Vertices have the Same Coordinates
bool Vertex::same(const Vertex& a) const {
    const float EPSILON = 1e-6f;
    return (std::abs(this->x - a.x) < EPSILON) &&
           (std::abs(this->y - a.y) < EPSILON) &&
           (std::abs(this->z - a.z) < EPSILON);
}

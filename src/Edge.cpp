// Edge.cpp

#include "Edge.h"
#include <cmath>    // For sqrt and pow
#include <iostream> // For printing messages

// Include GLM headers
#include <glm/glm.hpp>          // For glm::vec3
#include <glm/gtx/norm.hpp>     // For glm::length, glm::length2
#include <glm/geometric.hpp>    // For glm::dot

// Default Constructor
Edge::Edge() {
    // Initializes the vertices to default-constructed Vertices
    this->a = Vertex();
    this->b = Vertex();
    this->eno = -1;         // Indicates an undefined or invalid edge number
    this->hidden = false;   // Marks the edge as visible by default
    this->isTrue = true;    // Marks the edge as active/valid
}

// Parameterized Constructor with Vertices and Edge Number
Edge::Edge(const Vertex& a, const Vertex& b, int n) {
    this->a = a;            // Assigns the first vertex
    this->b = b;            // Assigns the second vertex
    this->eno = n;          // Assigns the provided edge number
    this->hidden = false;   // Marks the edge as visible by default
    this->isTrue = true;    // Marks the edge as active/valid

    const float EPSILON = 1e-6f; // Or use Model_2D::EPS

    if (std::abs(a.x - b.x) < EPSILON &&
    std::abs(a.y - b.y) < EPSILON &&
    std::abs(a.z - b.z) < EPSILON) {
    std::cout << "Same vertices " << a.vNo << " and " << b.vNo << ". No edge added!!" << std::endl;
    this->isTrue = false; // Marks the edge as invalid if both vertices are the same
}
}

// Retrieve the Neighboring Vertex of a Given Vertex in the Edge
Vertex Edge::getNeighbour(const Vertex& d) const {
    // Given a vertex 'd', this function returns the other vertex connected by the edge.
    if (this->a.vNo == d.vNo) {
        return this->b;  // If 'd' is the first vertex, return the second vertex
    } else {
        return this->a;  // Otherwise, return the first vertex
    }
}

// Calculate and Return the Midpoint Vertex of the Edge
Vertex Edge::midpoint() const {
    // Computes the geometric midpoint of the edge by averaging the coordinates of the two vertices.
    float midX = (this->a.x + this->b.x) / 2.0f;
    float midY = (this->a.y + this->b.y) / 2.0f;
    float midZ = (this->a.z + this->b.z) / 2.0f;

    // The vertex number is set to -1, indicating that it's a new, derived vertex.
    return Vertex(midX, midY, midZ, -1);
}

// Return the Current Edge Object
Edge Edge::getEdge() const {
    // Provides access to the current Edge object.
    return *this;
}

// Determine if Another Edge Overlaps with the Current Edge
bool Edge::overlap(const Edge& e2) const {
    // Calculate the direction vectors of both edges
    glm::vec3 dir1(this->a.x - this->b.x, this->a.y - this->b.y, this->a.z - this->b.z);
    glm::vec3 dir2(e2.a.x - e2.b.x, e2.a.y - e2.b.y, e2.a.z - e2.b.z);

    // Calculate the dot product of the direction vectors
    float res = glm::dot(dir1, dir2);

    // Calculate the magnitudes of both edge vectors
    float a1 = glm::length(dir1);
    float a2 = glm::length(dir2);

    // Normalize the dot product to obtain the cosine of the angle between the edges
    if (a1 == 0.0f || a2 == 0.0f) {
        return false; // Avoid division by zero
    }

    res = fabs(res / (a1 * a2));

    // Check if the edges are parallel (cosine near 1) and share at least one common vertex
    if ((res >= 0.99f && res <= 1.01f) &&
        (e2.a.same(this->a) || e2.a.same(this->b) ||
         e2.b.same(this->a) || e2.b.same(this->b))) {
        return true; // Edges overlap
    }

    return false;
}
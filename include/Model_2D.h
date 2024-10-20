// Model_2D.h

#ifndef MODEL_2D_H
#define MODEL_2D_H

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <glm/glm.hpp> // For GLM vector operations
#include <fstream>     // For file input

#include "Vertex.h"
#include "Edge.h"
#include "Surface.h"

class Model_2D {
public:
    // Projection Directions
    enum Direction {
        TOP_VIEW = 0,   // Top view (x, y)
        FRONT_VIEW = 1, // Front view (x, z)
        SIDE_VIEW = 2   // Side view (y, z)
    };

    // Constants
    static constexpr float EPS = 0.0001f; // Epsilon for floating-point comparisons
    static constexpr float INF = 1e9f;    // Infinity value for ray casting

    // Attributes
    int direction;               // Projection direction
    std::vector<Vertex> v;       // List of vertices
    std::vector<Edge> e;         // List of edges
    std::vector<Surface> s;      // List of surfaces (if applicable)

    // Constructors
    Model_2D();

    // Methods
    void add_vertices(const std::vector<Vertex>& vt);
    void add_edges(const std::vector<Edge>& et);
    void add_vertex(const Vertex& vt);
    const std::vector<Vertex>& ret_vertices() const;
    const std::vector<Edge>& ret_edges() const;
    Vertex* is_intersect(int p, int q);
    bool onSegment(const glm::vec2& p, const glm::vec2& q, const glm::vec2& r) const;
    int orientation(const glm::vec2& p, const glm::vec2& q, const glm::vec2& r) const;
    bool doIntersect(const glm::vec2& p1, const glm::vec2& q1,
                     const glm::vec2& p2, const glm::vec2& q2) const;
    bool is_inside(const Vertex& t, int q) const;

    // New method to read from file
    void read_from_file(std::ifstream& inFile);
};

#endif // MODEL_2D_H
#ifndef FACE_H
#define FACE_H

#include <vector>
#include "edge.h"

class Face {
public:
    std::vector<Edge> edges;

    Face(const std::vector<Edge>& edges) : edges(edges) {}

    // Add edge to the face
    void addEdge(const Edge& edge) {
        edges.push_back(edge);
    }
};

#endif // FACE_H
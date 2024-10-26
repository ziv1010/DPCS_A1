#ifndef OBJECT3D_H
#define OBJECT3D_H

#include <vector>
#include "vertex.h"
#include "edge.h"
#include "face.h"

class Object3D {
public:
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;

    // Add a vertex
    void addVertex(const Vertex& v) {
        vertices.push_back(v);
    }

    // Add an edge
    void addEdge(const Edge& e) {
        edges.push_back(e);
    }

    // Add a face
    void addFace(const Face& f) {
        faces.push_back(f);
    }
};

#endif // OBJECT3D_H
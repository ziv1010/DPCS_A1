// add_elements.h
#ifndef ADD_ELEMENTS_H
#define ADD_ELEMENTS_H

#include "model_2d_core.h"
#include "Vertex.h"
#include "Edge.h"
#include <vector>

// Adding elements functionalities for Model2D
class Model2DAddElements : public Model2DCore {
public:
    // Function to add a list of vertices to the model
    // Input: const std::vector<Vertex>& newVertices
    void addVertices(const std::vector<Vertex>& newVertices) {
        for (const auto& vertex : newVertices) {
            vertices.push_back(vertex);
        }
    }

    // Function to add a list of edges to the model
    // Input: const std::vector<Edge>& newEdges
    void addEdges(const std::vector<Edge>& newEdges) {
        for (const auto& edge : newEdges) {
            edges.push_back(edge);
        }
    }

    // Function to add a single vertex to the model
    // Input: const Vertex& newVertex
    void addVertex(const Vertex& newVertex) {
        vertices.push_back(newVertex);
    }
};

#endif // ADD_ELEMENTS_H
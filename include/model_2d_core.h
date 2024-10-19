// model_2d_core.h
#ifndef MODEL_2D_CORE_H
#define MODEL_2D_CORE_H

#include <vector>
#include "Vertex.h"
#include "Edge.h"
#include "Surface.h"

// Core functionalities and member variables of Model2D
// In model_2d_core.h
class Model2DCore {
public:
    // Member variables
    int direction;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Surface> surfaces;

    // Default constructor
    Model2DCore()
        : direction(-1), vertices(), edges(), surfaces() {}

    // Function to add a single vertex to the model
    void addVertex(const Vertex& newVertex) {
        vertices.push_back(newVertex);
    }
    // Member variables
    int direction;                                 // Direction of the model
    std::vector<Vertex> vertices;                  // List of vertices in the model
    std::vector<Edge> edges;                       // List of edges in the model
    std::vector<Surface> surfaces;                 // List of surfaces in the model

    // Default constructor
    Model2DCore()
        : direction(0), vertices(), edges(), surfaces() {}

    // Function to retrieve the list of vertices
    const std::vector<Vertex>& getVertices() const {
        return vertices;
    }

    // Function to retrieve the list of edges
    const std::vector<Edge>& getEdges() const {
        return edges;
    }

    // Function to return the model object itself
    Model2DCore* getModel() {
        return this;
    }
};

#endif // MODEL_2D_CORE_H
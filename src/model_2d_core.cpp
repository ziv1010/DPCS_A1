// model_2d_core.cpp
#include "model_2d_core.h"

// Constructor to initialize the model
Model2DCore::Model2DCore()
    : direction(-1), vertices(), edges(), surfaces() {
    // Initializes the direction and empty lists for vertices, edges, and surfaces
}

// Function to retrieve the list of vertices
const std::vector<Vertex>& Model2DCore::getVertices() const {
    return vertices;
}

// Function to retrieve the list of edges
const std::vector<Edge>& Model2DCore::getEdges() const {
    return edges;
}

// Function to return the model object itself
Model2DCore* Model2DCore::getModel() {
    return this;
}
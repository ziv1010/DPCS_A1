// add_elements.cpp
#include "add_elements.h"

// Function to add a list of vertices to the model
void Model2DAddElements::addVertices(const std::vector<Vertex>& newVertices) {
    vertices.insert(vertices.end(), newVertices.begin(), newVertices.end());
}

// Function to add a list of edges to the model
void Model2DAddElements::addEdges(const std::vector<Edge>& newEdges) {
    edges.insert(edges.end(), newEdges.begin(), newEdges.end());
}

// Function to add a single vertex to the model
void Model2DAddElements::addVertex(const Vertex& newVertex) {
    vertices.push_back(newVertex);
}

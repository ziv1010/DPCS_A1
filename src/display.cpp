// display.cpp
#include "display.h"
#include <iostream>

void Model2DDisplay::displayModel() const {
    displayVertices();
    displayEdges();
    displaySurfaces();
}

void Model2DDisplay::displayVertices() const {
    std::cout << "Vertices (Direction: " << direction << "):\n";
    for (const auto& vertex : vertices) {
        std::cout << "Vertex ID: " << vertex.vertexID
                  << " | Coordinates: (" << vertex.x << ", "
                  << vertex.y << ", " << vertex.z << ")\n";
    }
}

void Model2DDisplay::displayEdges() const {
    std::cout << "Edges (Direction: " << direction << "):\n";
    for (const auto& edge : edges) {
        std::cout << "Edge ID: " << edge.edgeID
                  << " | From Vertex ID: " << (edge.A ? edge.A->vertexID : -1)
                  << " To Vertex ID: " << (edge.B ? edge.B->vertexID : -1)
                  << " | Hidden: " << (edge.isHidden ? "Yes" : "No")
                  << " | Real: " << (edge.isReal ? "Yes" : "No") << "\n";
    }
}

void Model2DDisplay::displaySurfaces() const {
    std::cout << "Surfaces (Direction: " << direction << "):\n";
    for (const auto& surface : surfaces) {
        std::cout << "Surface ID: " << surface.surfaceID << " | Boundary Edges: ";
        for (const auto& edge : surface.boundaryEdges) {
            std::cout << edge->edgeID << " ";
        }
        std::cout << "\nCoefficients: ("
                  << surface.coefficients[0] << ", "
                  << surface.coefficients[1] << ", "
                  << surface.coefficients[2] << ", "
                  << surface.coefficients[3] << ")\n";
    }
}
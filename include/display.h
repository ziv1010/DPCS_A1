// display.h
#ifndef DISPLAY_H
#define DISPLAY_H

#include "model_2d_core.h"
#include "Vertex.h"
#include "Edge.h"
#include "Surface.h"
#include <iostream>

// Display functionalities for Model2D
class Model2DDisplay : virtual public Model2DCore {
public:
    // Function to display the model's structure (vertices, edges, surfaces)
    void displayModel() const {
        displayVertices();
        displayEdges();
        displaySurfaces();
    }

    // Function to display the vertices in the model
    void displayVertices() const {
        std::cout << "Vertices:\n";
        for (const auto& vertex : vertices) {
            std::cout << "Vertex ID: " << vertex.vertexID
                      << " | Coordinates: (" << vertex.x << ", "
                      << vertex.y << ", " << vertex.z << ")\n";
        }
    }

    // Function to display the edges in the model
    void displayEdges() const {
        std::cout << "Edges:\n";
        for (const auto& edge : edges) {
            std::cout << "Edge ID: " << edge.edgeID
                      << " | From Vertex ID: " << (edge.A ? edge.A->vertexID : -1)
                      << " To Vertex ID: " << (edge.B ? edge.B->vertexID : -1) << "\n";
        }
    }

    // Function to display the surfaces in the model
    void displaySurfaces() const {
        std::cout << "Surfaces:\n";
        for (const auto& surface : surfaces) {
            std::cout << "Surface ID: " << surface.surfaceID << " | Boundary Edges: ";
            for (const auto& edge : surface.boundaryEdges) {
                std::cout << edge->edgeID << " ";
            }
            std::cout << "\n";
        }
    }
};

#endif // DISPLAY_H
#ifndef EDGE_H
#define EDGE_H

#include "Vertex.h" // Ensure that Vertex.h is included or forward declared

// Edge class represents an edge between two vertices in a graph or model
class Edge {
public:
    // Member variables

    Vertex* A;         // Pointer to vertex A (one end of the edge)
    Vertex* B;         // Pointer to vertex B (other end of the edge)
    bool isHidden;     // Flag to indicate if the edge is hidden
    bool isReal;       // Flag to indicate if the edge is valid (exists in the model)
    int edgeID;        // Unique identifier for the edge

    // Default constructor
    // Initializes member variables to default values
    Edge()
        : A(nullptr), B(nullptr), isHidden(false), isReal(false), edgeID(0) {}

    // Parameterized constructor
    // Initializes the edge with given vertices and edge ID
    Edge(Vertex* vertexA, Vertex* vertexB, int id)
        : A(vertexA), B(vertexB), isHidden(false), isReal(true), edgeID(id) {}

    // Get the neighboring vertex opposite to the current vertex
    // Returns a pointer to the neighboring vertex or nullptr if the current vertex is not part of the edge
    Vertex* getNeighbor(Vertex* currentVertex) const {
        if (currentVertex == A) {
            return B;
        } else if (currentVertex == B) {
            return A;
        } else {
            return nullptr; // Current vertex is not part of this edge
        }
    }

    // Calculate and return the midpoint of the edge
    // Returns a new Vertex object representing the midpoint
    Vertex* midpoint() const {
        Vertex* newVertex = new Vertex();
        newVertex->x = (A->x + B->x) / 2.0;
        newVertex->y = (A->y + B->y) / 2.0;
        newVertex->z = (A->z + B->z) / 2.0;
        return newVertex;
    }

    // Check if this edge overlaps with another edge
    // Returns true if both ends of the edges match, false otherwise
    bool overlap(const Edge& otherEdge) const {
        bool endsMatch = 
            ((A == otherEdge.A || A == otherEdge.B) &&
             (B == otherEdge.A || B == otherEdge.B));
        return endsMatch;
    }

    // Get the edge object itself
    // Returns a pointer to this edge
    Edge* getEdge() {
        return this;
    }
};

#endif // EDGE_H
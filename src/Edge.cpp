#include "Edge.h"
#include <glm/glm.hpp>
#include <cmath>

Edge::Edge()
    : A(nullptr), B(nullptr), isHidden(false), isReal(false), edgeID(-1) {
    // Initializes an edge with default values.
}

Edge::Edge(Vertex* vertexA, Vertex* vertexB, int id)
    : A(vertexA), B(vertexB), isHidden(false), isReal(true), edgeID(id) {
    // Initializes an edge with specified vertices and edge ID.

    // If both vertices have identical coordinates, mark the edge as invalid.
    if (A->same(*B)) {
        isReal = false;
    }
}

Vertex* Edge::getNeighbor(Vertex* currentVertex) const {
    // Returns the vertex at the opposite end of the edge from the given vertex.
    if (currentVertex == A) {
        return B;
    } else if (currentVertex == B) {
        return A;
    } else {
        return nullptr; // Current vertex is not part of this edge.
    }
}

Vertex* Edge::midpoint() const {
    // Returns a new vertex representing the midpoint of the edge.
    float midX = (A->x + B->x) / 2.0f;
    float midY = (A->y + B->y) / 2.0f;
    float midZ = (A->z + B->z) / 2.0f;
    return new Vertex(midX, midY, midZ, -1);
}

bool Edge::overlap(const Edge& otherEdge) const {
    // Determines if two edges overlap by checking parallelism and shared vertices.

    // Create direction vectors using GLM.
    glm::vec3 dir1(A->x - B->x, A->y - B->y, A->z - B->z);
    glm::vec3 dir2(otherEdge.A->x - otherEdge.B->x, otherEdge.A->y - otherEdge.B->y, otherEdge.A->z - otherEdge.B->z);

    // Calculate the cosine of the angle between the edges.
    float dotProduct = glm::dot(dir1, dir2);
    float magnitude1 = glm::length(dir1);
    float magnitude2 = glm::length(dir2);

    if (magnitude1 == 0.0f || magnitude2 == 0.0f) {
        return false; // One of the edges has zero length.
    }

    float cosineAngle = std::abs(dotProduct / (magnitude1 * magnitude2));

    // Check if edges are nearly parallel (cosine close to 1) and share a common vertex.
    if (cosineAngle > 0.99f && cosineAngle < 1.01f &&
        (A->same(*otherEdge.A) || A->same(*otherEdge.B) || B->same(*otherEdge.A) || B->same(*otherEdge.B))) {
        return true;
    }

    return false;
}

Edge* Edge::getEdge() {
    // Returns the current edge object.
    return this;
}
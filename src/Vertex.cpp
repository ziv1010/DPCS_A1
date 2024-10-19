#include "Vertex.h"
#include <glm/glm.hpp>

Vertex::Vertex()
    : x(0.0f), y(0.0f), z(0.0f), isReal(false), vertexID(-1) {
    // Initializes a vertex with default coordinates and invalid status.
}

Vertex::Vertex(float xCoord, float yCoord, float zCoord, int id)
    : x(xCoord), y(yCoord), z(zCoord), isReal(true), vertexID(id) {
    // Initializes a vertex with specified coordinates and vertex ID.
}

Vertex::Vertex(const Vertex& v)
    : x(v.x), y(v.y), z(v.z), isReal(v.isReal), vertexID(-1) {
    // Creates a copy of the given vertex with a default vertex ID.
}

Vertex Vertex::getXY() const {
    // Returns a new vertex with the XY coordinates for top view projection.
    return Vertex(x, y, 0.0f, vertexID);
}

Vertex Vertex::getYZ() const {
    // Returns a new vertex with the YZ coordinates for side view projection.
    return Vertex(0.0f, y, z, vertexID);
}

Vertex Vertex::getXZ() const {
    // Returns a new vertex with the XZ coordinates for front view projection.
    return Vertex(x, 0.0f, z, vertexID);
}

bool Vertex::same(const Vertex& a) const {
    // Checks if the current vertex has the same coordinates as another vertex.
    return (x == a.x) && (y == a.y) && (z == a.z);
}
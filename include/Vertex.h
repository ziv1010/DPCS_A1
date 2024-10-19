#ifndef VERTEX_H
#define VERTEX_H

// Vertex class represents a point in 3D space with coordinates and identification
class Vertex {
public:
    // Member variables
    float x, y, z;     // Coordinates of the vertex in 3D space
    bool isReal;       // Flag to indicate if the vertex is valid (exists in the model)
    int vertexID;      // Unique identifier for the vertex

    // Default constructor
    // Initializes member variables to default values
    Vertex()
        : x(0.0f), y(0.0f), z(0.0f), isReal(false), vertexID(0) {}

    // Parameterized constructor
    // Initializes the vertex with given coordinates and vertex ID
    Vertex(float xCoord, float yCoord, float zCoord, int id)
        : x(xCoord), y(yCoord), z(zCoord), isReal(true), vertexID(id) {}

    // Copy constructor
    // Creates a new vertex by copying another vertex
    Vertex(const Vertex& v)
        : x(v.x), y(v.y), z(v.z), isReal(v.isReal), vertexID(v.vertexID) {}

    // Return the projection of the vertex on the XY plane
    // Returns a new Vertex object with z-coordinate set to 0
    Vertex getXY() const {
        return Vertex(x, y, 0.0f, vertexID);
    }

    // Return the projection of the vertex on the YZ plane
    // Returns a new Vertex object with x-coordinate set to 0
    Vertex getYZ() const {
        return Vertex(0.0f, y, z, vertexID);
    }

    // Return the projection of the vertex on the XZ plane
    // Returns a new Vertex object with y-coordinate set to 0
    Vertex getXZ() const {
        return Vertex(x, 0.0f, z, vertexID);
    }

    // Check if two vertices are the same (have the same coordinates)
    // Returns true if coordinates match, false otherwise
    bool same(const Vertex& a) const {
        return (x == a.x) && (y == a.y) && (z == a.z);
    }
};

#endif // VERTEX_H
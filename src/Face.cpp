// Face.cpp

#include "Face.h"
#include <cmath>    // For sqrt and pow
#include <iostream> // For printing messages

// Default Constructor
Surface::Surface() {
    // Initialize the boundary as an empty vector
    this->boundary = std::vector<Edge>();

    this->trueSurface = false;   // Marks the surface as inactive or invalid
    coeff[0] = 0.0f;             // Initializes coefficient 'a' to 0
    coeff[1] = 0.0f;             // Initializes coefficient 'b' to 0
    coeff[2] = 0.0f;             // Initializes coefficient 'c' to 0
    coeff[3] = 0.0f;             // Initializes coefficient 'd' to 0
    this->sno = -1;              // Indicates an undefined or invalid surface number
}

// Parameterized Constructor with Boundary Edges and Surface Number
Surface::Surface(const std::vector<Edge>& a, int n) {
    if (a.size() < 3) {
        std::cout << "Number of edges can't be less than 3. Wrong input!" << std::endl;
        this->trueSurface = false;
        return; // Exits the constructor early due to insufficient edges
    }

    this->boundary = a;          // Assigns the provided list of edges as the boundary
    this->calc_coeff();          // Calculates the plane coefficients based on the boundary edges
    this->trueSurface = true;    // Marks the surface as active/valid
    this->sno = n;               // Assigns the provided surface number identifier
}

// Calculate the Coefficients of the Plane Equation for the Surface
void Surface::calc_coeff() {
    // Identify three distinct vertices to define the plane
    Vertex m = boundary[0].a;  // First vertex from the first edge
    Vertex n = boundary[0].b;  // Second vertex from the first edge
    Vertex o = boundary[1].a;  // Third vertex from the second edge

    // Ensure that the third vertex 'o' is not the same as 'm' or 'n'
    if ((o.x == m.x && o.y == m.y && o.z == m.z) ||
        (o.x == n.x && o.y == n.y && o.z == n.z)) {
        o = boundary[1].b; // If 'o' is the same as 'm' or 'n', choose the other vertex
    }

    // Calculate two vectors lying on the plane
    float x1 = m.x - n.x;  // X-component of the first vector
    float y1 = m.y - n.y;  // Y-component of the first vector
    float z1 = m.z - n.z;  // Z-component of the first vector

    float x2 = o.x - n.x;  // X-component of the second vector
    float y2 = o.y - n.y;  // Y-component of the second vector
    float z2 = o.z - n.z;  // Z-component of the second vector

    // Calculate the normal vector (a, b, c) of the plane using the cross product
    coeff[0] = (y1 * z2) - (y2 * z1);  // Coefficient 'a'
    coeff[1] = (x2 * z1) - (x1 * z2);  // Coefficient 'b'
    coeff[2] = (x1 * y2) - (x2 * y1);  // Coefficient 'c'

    // Calculate the coefficient 'd' using the plane equation ax + by + cz + d = 0
    coeff[3] = -1.0f * ((coeff[0] * n.x) + (coeff[1] * n.y) + (coeff[2] * n.z));

    // Normalize the coefficients to ensure that the normal vector has unit length
    float magnitude = sqrt((coeff[0] * coeff[0]) + (coeff[1] * coeff[1]) + (coeff[2] * coeff[2]));

    if (magnitude != 0.0f) {
        coeff[0] /= magnitude; // Normalize coefficient 'a'
        coeff[1] /= magnitude; // Normalize coefficient 'b'
        coeff[2] /= magnitude; // Normalize coefficient 'c'
        coeff[3] /= magnitude; // Normalize coefficient 'd'
    }
}

// Calculate the Projection of a Given Vertex onto the Surface Plane
float Surface::calcProj(const Vertex& v_test) const {
    // Evaluates the plane equation ax + by + cz + d for the given vertex coordinates.
    // The result indicates the relative position of the vertex concerning the surface plane.
    float result = (coeff[0] * v_test.x) + (coeff[1] * v_test.y) + (coeff[2] * v_test.z) + coeff[3];
    return result;
}

// Calculate the Dot Product between the Surface's Normal Vector and an Edge Vector
float Surface::dotProduct(const Edge& e_t) const {
    // Computes the dot product between the normal vector of the surface and the vector representing the edge.
    // This value can be used to determine the orientation or relative positioning of the edge with respect to the surface.
    float edgeVectorX = e_t.a.x - e_t.b.x; // X-component of the edge vector
    float edgeVectorY = e_t.a.y - e_t.b.y; // Y-component of the edge vector
    float edgeVectorZ = e_t.a.z - e_t.b.z; // Z-component of the edge vector

    float result = (coeff[0] * edgeVectorX) +
                   (coeff[1] * edgeVectorY) +
                   (coeff[2] * edgeVectorZ);
    return result;
}
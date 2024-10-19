#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include "Edge.h"    // Ensure that Edge.h is included since Surface depends on Edge

// Surface class represents a planar surface defined by its boundary edges and plane equation
class Surface {
public:
    // Member variables
    std::vector<Edge*> boundaryEdges;   // List of edges that form the boundary of the surface
    float coefficients[4];              // Coefficients of the plane equation (Ax + By + Cz + D = 0)
    bool isRealSurface;                 // Flag to indicate if the surface is valid/real
    int surfaceID;                      // Unique identifier for the surface

    // Default constructor
    // Initializes member variables to default values
    Surface()
        : boundaryEdges(), isRealSurface(false), surfaceID(0) {
        coefficients[0] = 0.0f; // A
        coefficients[1] = 0.0f; // B
        coefficients[2] = 0.0f; // C
        coefficients[3] = 0.0f; // D
    }

    // Parameterized constructor
    // Initializes the surface with given boundary edges and surface ID
    Surface(const std::vector<Edge*>& edges, int id)
        : boundaryEdges(edges), isRealSurface(true), surfaceID(id) {
        coefficients[0] = 0.0f; // Initialize coefficients; should be calculated later
        coefficients[1] = 0.0f;
        coefficients[2] = 0.0f;
        coefficients[3] = 0.0f;
    }

    // Get the list of edges that form the boundary of the surface
    // Returns a constant reference to the boundary edges vector
    const std::vector<Edge*>& getEdges() const {
        return boundaryEdges;
    }

    // Calculate the coefficients of the surface's plane equation (Ax + By + Cz + D = 0)
    void calculateCoefficients() {
        // Implement logic to calculate coefficients based on boundary edges or vertices
        // This is a placeholder implementation and should be replaced with actual logic

        // Example Placeholder:
        // Assume the first three vertices of the first three edges are non-collinear
        if (boundaryEdges.size() < 3) {
            // Not enough edges to define a plane
            coefficients[0] = coefficients[1] = coefficients[2] = coefficients[3] = 0.0f;
            return;
        }

        // Extract three distinct vertices
        Vertex* v1 = boundaryEdges[0]->A;
        Vertex* v2 = boundaryEdges[1]->A;
        Vertex* v3 = boundaryEdges[2]->A;

        // Vectors u and v
        float ux = v2->x - v1->x;
        float uy = v2->y - v1->y;
        float uz = v2->z - v1->z;

        float vx = v3->x - v1->x;
        float vy = v3->y - v1->y;
        float vz = v3->z - v1->z;

        // Cross product u x v to get normal vector (A, B, C)
        coefficients[0] = uy * vz - uz * vy; // A
        coefficients[1] = uz * vx - ux * vz; // B
        coefficients[2] = ux * vy - uy * vx; // C

        // Calculate D using point v1
        coefficients[3] = -(coefficients[0] * v1->x +
                           coefficients[1] * v1->y +
                           coefficients[2] * v1->z); // D
    }

    // Calculate the projection of a vertex onto the surface
    // Returns the signed distance from the vertex to the plane
    float calculateProjection(const Vertex& vTest) const {
        // Using the plane equation Ax + By + Cz + D
        return (coefficients[0] * vTest.x +
                coefficients[1] * vTest.y +
                coefficients[2] * vTest.z +
                coefficients[3]);
    }

    // Calculate the dot product between an edge and the surface's normal vector
    // Returns the dot product result
    float dotProduct(const Edge& edgeTest) const {
        // Calculate the edge vector (B - A)
        float edgeX = edgeTest.B->x - edgeTest.A->x;
        float edgeY = edgeTest.B->y - edgeTest.A->y;
        float edgeZ = edgeTest.B->z - edgeTest.A->z;

        // Dot product with the surface's normal vector (A, B, C)
        return (coefficients[0] * edgeX +
                coefficients[1] * edgeY +
                coefficients[2] * edgeZ);
    }
};

#endif // SURFACE_H


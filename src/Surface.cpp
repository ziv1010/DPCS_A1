#include "Surface.h"
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <iostream>

void Surface::calculateCoefficients() {
    // Calculates the plane coefficients (A, B, C, D) for the surface equation Ax + By + Cz + D = 0.

    if (boundaryEdges.size() < 3) {
        std::cerr << "Number of edges can't be less than 3. Wrong input!" << std::endl;
        isRealSurface = false;
        return;
    }

    // Extract three distinct vertices from the boundary edges.
    Vertex* m = boundaryEdges[0]->A;
    Vertex* n = boundaryEdges[0]->B;
    Vertex* o = boundaryEdges[1]->A;

    if (o->same(*m) || o->same(*n)) {
        o = boundaryEdges[1]->B;
    }

    // Calculate vectors using GLM.
    glm::vec3 u(m->x - n->x, m->y - n->y, m->z - n->z);
    glm::vec3 v(o->x - n->x, o->y - n->y, o->z - n->z);

    // Compute the normal vector (A, B, C) via cross product.
    glm::vec3 normal = glm::cross(u, v);

    if (glm::length2(normal) == 0.0f) {
        std::cerr << "Edges are colinear, cannot define a plane." << std::endl;
        isRealSurface = false;
        return;
    }

    // Normalize the normal vector.
    normal = glm::normalize(normal);

    // Assign the coefficients.
    coefficients[0] = normal.x; // A
    coefficients[1] = normal.y; // B
    coefficients[2] = normal.z; // C

    // Calculate D using point n.
    coefficients[3] = -glm::dot(normal, glm::vec3(n->x, n->y, n->z));

    isRealSurface = true;
}

float Surface::calculateProjection(const Vertex& vTest) const {
    // Calculates the signed distance from the vertex to the plane.
    glm::vec3 normal(coefficients[0], coefficients[1], coefficients[2]);
    glm::vec3 point(vTest.x, vTest.y, vTest.z);
    return glm::dot(normal, point) + coefficients[3];
}

float Surface::dotProduct(const Edge& edgeTest) const {
    // Calculates the dot product of the edge's direction vector with the surface's normal vector.
    glm::vec3 edgeVec(edgeTest.B->x - edgeTest.A->x,
                      edgeTest.B->y - edgeTest.A->y,
                      edgeTest.B->z - edgeTest.A->z);
    glm::vec3 normal(coefficients[0], coefficients[1], coefficients[2]);
    return glm::dot(normal, edgeVec);
}
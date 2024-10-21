// Face.h

#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include "Edge.h"

class Surface {
public:
    // Attributes
    std::vector<Edge> boundary;
    float coeff[4];  // Coefficients for the plane equation
    bool trueSurface;
    int sno;

    // Constructors
    Surface();
    Surface(const std::vector<Edge>& a, int n);

    // Methods
    std::vector<Edge> getEdges() const;
    void calc_coeff();
    float calcProj(const Vertex& v_test) const;
    float dotProduct(const Edge& et) const;
};

#endif // SURFACE_H
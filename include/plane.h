// plane.h

#ifndef PLANE_H
#define PLANE_H

#include <vector>
#include "edge3d.h"

class Plane {
public:
    float a, b, c, d; // Plane equation coefficients: ax + by + cz + d = 0
    std::vector<int> edgesOnPlane; // Indices of edges lying on this plane

    Plane(float a_, float b_, float c_, float d_)
        : a(a_), b(b_), c(c_), d(d_) {}

    // Method to compare two planes within a tolerance
    bool isSimilar(const Plane& other, float tolerance = 1e-5f) const {
        float diff = std::sqrt((a - other.a) * (a - other.a) +
                               (b - other.b) * (b - other.b) +
                               (c - other.c) * (c - other.c) +
                               (d - other.d) * (d - other.d));
        return diff <= tolerance;
    }
};

#endif // PLANE_H

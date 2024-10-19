// geometry.h
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "model_2d_core.h"
#include "Vertex.h"
#include <utility> // For std::pair

// Geometry helper functionalities for Model2D
class Model2DGeometry : public Model2DCore {
public:
    // Function to check if a point lies on a given segment
    // Input: std::pair<float, float> p (point),
    //        std::pair<float, float> q (start of segment),
    //        std::pair<float, float> r (end of segment)
    // Output: bool (True if point lies on segment, otherwise False)
    bool isPointOnSegment(const std::pair<float, float>& p,
                          const std::pair<float, float>& q,
                          const std::pair<float, float>& r) const {
        // Implement logic to check if point p lies on segment qr
        // Placeholder implementation: Always return false
        return false;
    }

    // Function to calculate the orientation of three points (p, q, r)
    // Output: int (0 if collinear, 1 if clockwise, 2 if counterclockwise)
    int calculateOrientation(const std::pair<float, float>& p,
                             const std::pair<float, float>& q,
                             const std::pair<float, float>& r) const {
        // Implement logic to determine the orientation of the points
        // Placeholder implementation: Always return 0 (collinear)
        return 0;
    }

    // Function to check if a point is inside the polygon
    // Input: Vertex point, int numVertices (number of vertices in polygon)
    // Output: bool (True if point is inside the polygon, otherwise False)
    bool isPointInsidePolygon(const Vertex& point, int numVertices) const {
        // Implement logic to check if the point is inside the polygon formed by vertices
        // Placeholder implementation: Always return false
        return false;
    }
};

#endif // GEOMETRY_H
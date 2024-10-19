#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "model_2d_core.h"
#include "Vertex.h"
#include <glm/glm.hpp>

// Geometry helper functionalities for Model2D
class Model2DGeometry : virtual public Model2DCore {
public:
    // Function to check if a point lies on a given segment
    bool isPointOnSegment(const glm::vec2& p,
                          const glm::vec2& q,
                          const glm::vec2& r) const;

    // Function to calculate the orientation of three points (p, q, r)
    int calculateOrientation(const glm::vec2& p,
                             const glm::vec2& q,
                             const glm::vec2& r) const;

    // Function to check if two segments intersect
    bool doIntersect(const glm::vec2& p1, const glm::vec2& q1,
                     const glm::vec2& p2, const glm::vec2& q2) const;

    // Function to check if a point is inside the polygon
    bool isPointInsidePolygon(const Vertex& point, int surfaceIndex) const;
};

#endif // GEOMETRY_H
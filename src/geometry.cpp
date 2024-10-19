// geometry.cpp
#include "geometry.h"
#include <glm/glm.hpp>
#include <limits>
#include <iostream>

bool Model2DGeometry::isPointOnSegment(const glm::vec2& p, const glm::vec2& q, const glm::vec2& r) const {
    return (q.x <= glm::max(p.x, r.x) && q.x >= glm::min(p.x, r.x) &&
            q.y <= glm::max(p.y, r.y) && q.y >= glm::min(p.y, r.y));
}

int Model2DGeometry::calculateOrientation(const glm::vec2& p, const glm::vec2& q, const glm::vec2& r) const {
    float val = (q.y - p.y) * (r.x - q.x) - 
                (q.x - p.x) * (r.y - q.y);
    if (val == 0.0f) return 0;  // Colinear
    return (val > 0.0f) ? 1 : 2; // Clockwise or Counterclockwise
}

bool Model2DGeometry::doIntersect(const glm::vec2& p1, const glm::vec2& q1,
                                  const glm::vec2& p2, const glm::vec2& q2) const {
    int o1 = calculateOrientation(p1, q1, p2);
    int o2 = calculateOrientation(p1, q1, q2);
    int o3 = calculateOrientation(p2, q2, p1);
    int o4 = calculateOrientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    if (o1 == 0 && isPointOnSegment(p1, p2, q1)) return true;
    if (o2 == 0 && isPointOnSegment(p1, q2, q1)) return true;
    if (o3 == 0 && isPointOnSegment(p2, p1, q2)) return true;
    if (o4 == 0 && isPointOnSegment(p2, q1, q2)) return true;

    return false; // No intersection
}

bool Model2DGeometry::isPointInsidePolygon(const Vertex& point, int surfaceIndex) const {
    if (surfaceIndex < 0 || surfaceIndex >= surfaces.size()) {
        std::cerr << "Invalid surface index!" << std::endl;
        return false;
    }

    const Surface& surface = surfaces[surfaceIndex];
    int n = static_cast<int>(surface.boundaryEdges.size());
    if (n < 3) {
        std::cerr << "Invalid surface: Less than 3 edges!" << std::endl;
        return false;
    }

    std::vector<glm::vec2> polygonPoints(n);

    // Project the point and polygon based on direction
    glm::vec2 p;
    switch (direction) {
        case 0: // Top view (x, y)
            p = glm::vec2(point.x, point.y);
            for (int i = 0; i < n; ++i) {
                polygonPoints[i] = glm::vec2(surface.boundaryEdges[i]->A->x, surface.boundaryEdges[i]->A->y);
            }
            break;
        case 1: // Front view (x, z)
            p = glm::vec2(point.x, point.z);
            for (int i = 0; i < n; ++i) {
                polygonPoints[i] = glm::vec2(surface.boundaryEdges[i]->A->x, surface.boundaryEdges[i]->A->z);
            }
            break;
        case 2: // Side view (y, z)
            p = glm::vec2(point.y, point.z);
            for (int i = 0; i < n; ++i) {
                polygonPoints[i] = glm::vec2(surface.boundaryEdges[i]->A->y, surface.boundaryEdges[i]->A->z);
            }
            break;
        default:
            return false;
    }

    // Ray Casting algorithm
    glm::vec2 extreme = glm::vec2(std::numeric_limits<float>::infinity(), p.y);
    int count = 0, i = 0;
    do {
        int next = (i + 1) % n;
        if (doIntersect(polygonPoints[i], polygonPoints[next], p, extreme)) {
            if (calculateOrientation(polygonPoints[i], p, polygonPoints[next]) == 0) {
                return isPointOnSegment(polygonPoints[i], p, polygonPoints[next]);
            }
            count++;
        }
        i = next;
    } while (i != 0);

    return (count % 2 == 1); // Point is inside if count is odd
}
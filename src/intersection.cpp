// intersection.cpp
#include "intersection.h"
#include "geometry.h" // Include this to use geometry functions
#include "Model2D.h"
#include <glm/glm.hpp>
#include <optional>

std::optional<Vertex> Model2DIntersection::checkIntersection(int edge1Index, int edge2Index) {
    if (edge1Index < 0 || edge1Index >= edges.size() ||
        edge2Index < 0 || edge2Index >= edges.size()) {
        return std::nullopt; // Invalid edge indices
    }

    Edge& edge1 = edges[edge1Index];
    Edge& edge2 = edges[edge2Index];

    glm::vec2 p1, p2, q1, q2;

    // Determine projection based on model direction
    switch (direction) {
        case 0: // Top view (x, y)
            p1 = glm::vec2(edge1.A->x, edge1.A->y);
            p2 = glm::vec2(edge1.B->x, edge1.B->y);
            q1 = glm::vec2(edge2.A->x, edge2.A->y);
            q2 = glm::vec2(edge2.B->x, edge2.B->y);
            break;
        case 1: // Front view (x, z)
            p1 = glm::vec2(edge1.A->x, edge1.A->z);
            p2 = glm::vec2(edge1.B->x, edge1.B->z);
            q1 = glm::vec2(edge2.A->x, edge2.A->z);
            q2 = glm::vec2(edge2.B->x, edge2.B->z);
            break;
        case 2: // Side view (y, z)
            p1 = glm::vec2(edge1.A->y, edge1.A->z);
            p2 = glm::vec2(edge1.B->y, edge1.B->z);
            q1 = glm::vec2(edge2.A->y, edge2.A->z);
            q2 = glm::vec2(edge2.B->y, edge2.B->z);
            break;
        default:
            return std::nullopt; // Undefined direction
    }

    // Check if the segments (p1,p2) and (q1,q2) intersect
    if (checkSegmentIntersection(p1, p2, q1, q2)) {
        // Calculate intersection point
        float denominator = (p1.x - p2.x) * (q1.y - q2.y) - (p1.y - p2.y) * (q1.x - q2.x);
        if (denominator == 0.0f) {
            return std::nullopt; // Lines are parallel
        }

        float x = ((p1.x * p2.y - p1.y * p2.x) * (q1.x - q2.x) - (p1.x - p2.x) * (q1.x * q2.y - q1.y * q2.x)) / denominator;
        float y = ((p1.x * p2.y - p1.y * p2.x) * (q1.y - q2.y) - (p1.y - p2.y) * (q1.x * q2.y - q1.y * q2.x)) / denominator;

        Vertex intersectionVertex;
        intersectionVertex.isReal = true;
        intersectionVertex.vertexID = static_cast<int>(vertices.size()) + 1;

        switch (direction) {
            case 0:
                intersectionVertex.x = x;
                intersectionVertex.y = y;
                intersectionVertex.z = 0.0f;
                break;
            case 1:
                intersectionVertex.x = x;
                intersectionVertex.y = 0.0f;
                intersectionVertex.z = y;
                break;
            case 2:
                intersectionVertex.x = 0.0f;
                intersectionVertex.y = x;
                intersectionVertex.z = y;
                break;
        }

        // Add intersection vertex to the model
        addVertex(intersectionVertex);

        // Split the edges at the intersection point
        // Create new edges by splitting edge1
        Edge newEdge1(edge1.A, &vertices.back(), static_cast<int>(edges.size()) + 1);
        Edge newEdge2(&vertices.back(), edge1.B, static_cast<int>(edges.size()) + 2);
        edges.push_back(newEdge1);
        edges.push_back(newEdge2);

        // Create new edges by splitting edge2
        Edge newEdge3(edge2.A, &vertices.back(), static_cast<int>(edges.size()) + 3);
        Edge newEdge4(&vertices.back(), edge2.B, static_cast<int>(edges.size()) + 4);
        edges.push_back(newEdge3);
        edges.push_back(newEdge4);

        // Mark original edges as invalid
        edge1.isReal = false;
        edge2.isReal = false;

        return intersectionVertex;
    }

    return std::nullopt; // No intersection
}

bool Model2DIntersection::checkSegmentIntersection(const glm::vec2& p1, const glm::vec2& q1,
                                                   const glm::vec2& p2, const glm::vec2& q2) {
    // Use geometryHelper to access Model2DGeometry functions
    int o1 = geometryHelper.calculateOrientation(p1, q1, p2);
    int o2 = geometryHelper.calculateOrientation(p1, q1, q2);
    int o3 = geometryHelper.calculateOrientation(p2, q2, p1);
    int o4 = geometryHelper.calculateOrientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases (colinear points)
    if (o1 == 0 && geometryHelper.isPointOnSegment(p1, p2, q1)) return true;
    if (o2 == 0 && geometryHelper.isPointOnSegment(p1, q2, q1)) return true;
    if (o3 == 0 && geometryHelper.isPointOnSegment(p2, p1, q2)) return true;
    if (o4 == 0 && geometryHelper.isPointOnSegment(p2, q1, q2)) return true;

    return false; // No intersection
}
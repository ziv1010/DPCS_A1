#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "model_2d_core.h"
#include "Vertex.h"
#include "Edge.h"
#include <optional>
#include <glm/glm.hpp>

// Intersection functionalities for Model2D
class Model2DIntersection : virtual public Model2DCore {
private:
    Model2DGeometry geometryHelper; // Add Model2DGeometry object
public:
    // Function to check if two edges intersect (returns a Vertex if they do)
    std::optional<Vertex> checkIntersection(int edge1Index, int edge2Index);

    // Function to check if two segments in 2D space intersect
    bool checkSegmentIntersection(const glm::vec2& p1, const glm::vec2& q1,
                                  const glm::vec2& p2, const glm::vec2& q2);
};

#endif // INTERSECTION_H
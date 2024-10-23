#ifndef WIREFRAME_H
#define WIREFRAME_H

#include <vector>
#include "vertex3d.h"    // Enhanced 3D vertex
#include "edge3d.h"      // 3D edge
#include "vertex2d.h"    // 2D vertex
#include "edge2d.h"      // 2D edge
#include "projection2d.h"

class Wireframe {
public:
    std::vector<Vertex3D> vertices3D;
    std::vector<Edge3D> edges3D;

    // Constructor
    Wireframe(const Projection2D& front, const Projection2D& top, const Projection2D& side);

    // Generate 3D vertices
    void generate3DVertices();

    // Generate 3D edges
    void generate3DEdges();

    // Remove redundant edges
    void removeRedundantEdges();

    // Remove pathological edges and vertices using PEVR
    void removePathologicalElements();

    // Save wireframe to file
    void saveToFile(const std::string& filename) const;

private:
    const Projection2D& frontView;
    const Projection2D& topView;
    const Projection2D& sideView;

    // Helper functions
    int findVertex3DIndex(float x, float y, float z) const;
    bool areCollinear(int edge1, int edge2) const;
    bool areCoplanar(const std::vector<int>& edgeIndices) const;
    void mergeEdges(int edge1, int edge2);
};

#endif // WIREFRAME_H

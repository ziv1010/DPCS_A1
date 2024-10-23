#ifndef WIREFRAME_H
#define WIREFRAME_H

#include <vector>
#include "vertex3d.h"    // 3D vertex
#include "edge3d.h"      // 3D edge
#include "vertex2d.h"    // 2D vertex
#include "edge2d.h"      // 2D edge
#include "projection2d.h"
#include "plane.h"       // New: Plane class
#include "vector3d.h"    // New: Vector3D class
#include <unordered_map>
#include "vector2d.h" 

// Forward declarations if necessary
class Wireframe {
public:
    std::vector<Vertex3D> vertices3D;
    std::vector<Edge3D> edges3D;
    std::vector<Plane> planes; // To store generated planes

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

    // Generate planar graphs
    void generatePlanarGraphs(); // New method

    // Save wireframe to file
    void saveToFile(const std::string& filename) const;

    void generateFaceLoops();

private:
    const Projection2D& frontView;
    const Projection2D& topView;
    const Projection2D& sideView;

    // Helper functions
    int findVertex3DIndex(float x, float y, float z) const;
    bool areCollinear(int edge1, int edge2) const;
    bool areCoplanar(const std::vector<int>& edgeIndices) const;
    void mergeEdges(int edge1, int edge2);

    // New Helper Functions for PGG
    bool isPlanarGraphClosed(const Plane& plane) const;
    void removePathologicalElementsAtVertex(int vertexIndex);
    void removeDanglingEdgesFromPlane(Plane& plane);
    void adjustIndicesAfterVertexRemoval(int removedVertexIndex);
    void adjustEdgeIndicesAfterEdgeRemoval(int removedEdgeIndex);

   // Helper functions for FLG
    std::vector<int> orderEdgesAtVertex(int vertexIdx, const std::vector<int>& connectedEdges, const Plane& plane);
    bool findBasicLoop(int startEdgeIdx, const Plane& plane,
                    const std::unordered_map<int, std::vector<int>>& orderedEdgesAtVertex,
                    std::vector<bool>& edgeVisited, std::vector<int>& basicLoop);
    std::vector<std::vector<bool>> computeInclusionMatrix(const Plane& plane);
    std::vector<std::vector<int>> formFaceLoops(const std::vector<std::vector<int>>& basicLoops,
                                                const std::vector<std::vector<bool>>& inclusionMatrix);
    bool isLoopIncludedInLoop(const std::vector<int>& innerLoop, const std::vector<int>& outerLoop, const Plane& plane);
    void projectLoopOntoPlane(const std::vector<int>& loop, const Plane& plane, std::vector<Vector2D>& loop2D);
    bool isPointInPolygon(const Vector2D& point, const std::vector<Vector2D>& polygon);

};

#endif // WIREFRAME_H

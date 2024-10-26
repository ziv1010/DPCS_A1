#ifndef ORTHOGRAPHIC_PROJECTIONS_H
#define ORTHOGRAPHIC_PROJECTIONS_H

#include "object3d.h"
#include <vector>
#include <utility>

struct Projection2D {
    std::vector<std::pair<float, float>> projectedVertices;
    std::vector<std::pair<int, int>> projectedEdges; // indices into projectedVertices
};

// Projection functions
void projectToTopView(const Object3D& object, Projection2D& projection);
void projectToFrontView(const Object3D& object, Projection2D& projection);
void projectToSideView(const Object3D& object, Projection2D& projection);

// Save projections as PNG
void saveCombinedProjectionAsImage(const std::string& filename,
                                   const Projection2D& topView,
                                   const Projection2D& frontView,
                                   const Projection2D& sideView);

// Save projections to text file
void saveProjectionsToTextFile(const std::string& filename,
                               const Projection2D& topView,
                               const Projection2D& frontView,
                               const Projection2D& sideView);

#endif

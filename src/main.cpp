#include "file_io.h"
#include "orthographic_projections.h"
#include "transformations.h"
#include <iostream>

int main() {
    Object3D object;

    // Read the 3D object from the input file
    read3DObjectFromFile("build/output/input3D.txt", object);

    // Perform transformations on the object if needed
    Transformations::rotateX(object, 30.0f);
    Transformations::rotateY(object, 45.0f);

    // Create Projection2D objects
    Projection2D topView, frontView, sideView;

    // Project the 3D object onto 2D planes
    projectToTopView(object, topView);   // Top view (XY plane)
    projectToFrontView(object, frontView); // Front view
    projectToSideView(object, sideView);  // Side view

    // Save the 2D projections as one combined PNG file
    saveCombinedProjectionAsImage("build/output/combined_views.png", topView, frontView, sideView);

    std::cout << "Projections saved to combined image." << std::endl;

    return 0;
}
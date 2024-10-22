#include "file_io.h"
#include "orthographic_projections.h"
#include "transformations.h"
#include <iostream>

int main() {
    Object3D object;

    // Read the 3D object from the input file
    read3DObjectFromFile("build/output/input3D.txt", object);



    // Create Projection2D objects
    Projection2D topView, frontView, sideView;

    // Project the 3D object onto 2D planes
    projectToTopView(object, topView);   // Top view (XY plane)
    projectToFrontView(object, frontView); // Front view
    projectToSideView(object, sideView);  // Side view

    // Save the 2D projections to text file
    saveProjectionsToTextFile("build/output/projections.txt", topView, frontView, sideView);

    // Save the projections as images
    saveProjectionAsImage("build/output/top_view.png", topView);
    saveProjectionAsImage("build/output/front_view.png", frontView);
    saveProjectionAsImage("build/output/side_view.png", sideView);

    std::cout << "Projections saved to files and images." << std::endl;

    return 0;
}

#include "file_io.h"
#include "orthographic_projections.h"
#include "transformations.h"
#include <iostream>
#include "slice.h"

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

    // Save the 2D projections as one combined PNG file
    saveCombinedProjectionAsImage("build/output/combined_views.png", topView, frontView, sideView);

    std::cout << "Projections saved to combined image." << std::endl;


     // Define the slicing plane (example: diagonal cut)
    Plane slicingPlane(1.0f, 0.0f, 0.0f, -2.0f);  // Ax + By + Cz = D

    // Create objects for the left and right sides after slicing
    Object3D leftSide, rightSide;

    // Slice the object
    sliceObject(object, slicingPlane, leftSide, rightSide);

    // Project the sliced objects
    Projection2D leftView, rightView;
    projectToTopView(leftSide, leftView);   // Top view (XY plane) for left side
    projectToTopView(rightSide, rightView); // Top view (XY plane) for right side

    // Combine the two projections into one PNG
    saveCombinedProjectionAsImage("build/output/sliced_views.png", leftView, rightView, rightView);

    std::cout << "Sliced object projections saved to file." << std::endl;

    return 0;
}
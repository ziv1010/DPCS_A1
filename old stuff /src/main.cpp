#include <iostream>
#include "vertex2d.h"
#include "edge2d.h"
#include "projection2d.h"
#include "wireframe.h"
#include <fstream>
#include <sstream>

void readProjectionFromFile(const std::string& filename, Projection2D& projection) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    int numVertices, numEdges;
    infile >> numVertices;
    for (int i = 0; i < numVertices; ++i) {
        float x, y;
        infile >> x >> y;
        projection.vertices.emplace_back(x, y);
    }

    infile >> numEdges;
    for (int i = 0; i < numEdges; ++i) {
        int v1, v2;
        infile >> v1 >> v2;
        projection.edges.emplace_back(v1, v2);
    }
    infile.close();
}

void saveWireframeToFile(const std::string& filename, const Wireframe& wireframe) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to write to file " << filename << std::endl;
        return;
    }

    // Output vertices
    int numVertices = wireframe.vertices3D.size();
    outfile << numVertices << std::endl;
    for (const auto& vertex : wireframe.vertices3D) {
        outfile << vertex.getX() << " " << vertex.getY() << " " << vertex.getZ() << std::endl;
    }

    // Output edges
    int numEdges = wireframe.edges3D.size();
    outfile << numEdges << std::endl;
    for (const auto& edge : wireframe.edges3D) {
        outfile << edge.v1_index << " " << edge.v2_index << std::endl;
    }

    outfile.close();
}

int main() {
    Projection2D frontView, topView, sideView;

    // Read projections from files
    readProjectionFromFile("build/output/front.txt", frontView);
    readProjectionFromFile("build/output/top.txt", topView);
    readProjectionFromFile("build/output/side.txt", sideView);

    // Create the wireframe model
    Wireframe wireframe(frontView, topView, sideView);

    // Save the wireframe to a file
    saveWireframeToFile("build/output/wireframe.txt", wireframe);

    std::cout << "3D wireframe reconstruction complete. Output saved to build/output/wireframe.txt" << std::endl;

    return 0;
}

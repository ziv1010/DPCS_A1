#include "file_io.h"
#include <iostream>
#include <fstream>

void read3DObjectFromFile(const std::string& filename, Object3D& object) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    int numVertices, numEdges;
    infile >> numVertices;
    for (int i = 0; i < numVertices; ++i) {
        float x, y, z;
        infile >> x >> y >> z;
        object.addVertex(Vertex(x, y, z));
    }
    
    infile >> numEdges;
    for (int i = 0; i < numEdges; ++i) {
        int v1, v2;
        infile >> v1 >> v2;
        object.addEdge(Edge(v1, v2));
    }
    infile.close();
}
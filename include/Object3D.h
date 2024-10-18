#ifndef OBJECT3D_H
#define OBJECT3D_H

#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include "Vertex.h"
#include "Edge.h"
#include "Loop.h"
#include "Face.h"
#include "Vector3D.h"
#include <GL/glew.h>  // Include GLEW for OpenGL functions

class Object3D {
public:
    GLuint VAO, VBO, EBO;  // Vertex Array Object, Vertex Buffer Object, Element Buffer Object
    std::unordered_map<int, Vertex> vertices;
    std::unordered_map<int, Edge> edges;
    std::unordered_map<int, Loop> loops;
    std::unordered_map<int, Face> faces;

    int vertexCounter;
    int edgeCounter;
    int loopCounter;
    int faceCounter;

public:
    Object3D();

    // Methods for adding and retrieving vertices
    int addVertex(float x, float y, float z);
    Vertex& getVertex(int id);

    // Vertex methods


    // Edge methods
    int addEdge(int startVertexID, int endVertexID);
    Edge& getEdge(int id);

    // Loop methods
    int addLoop();
    void addEdgeToLoop(int loopID, int edgeID);
    Loop& getLoop(int id);

    // Face methods
    int addFace();
    void addEdgeToFace(int faceID, int edgeID);
    Face& getFace(int id);

    // File loading
    void loadFromFile(const std::string& filename);

    // OpenGL initialization and drawing
    void initializeOpenGL();
    void draw();
};

#endif // OBJECT3D_H
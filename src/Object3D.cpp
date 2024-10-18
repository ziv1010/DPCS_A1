#include "Object3D.h"
#include <GL/glew.h>    // For OpenGL functions (after installing GLEW)
#include <glm/glm.hpp>  // For GLM vector and matrix functions
#include <glm/vec3.hpp> // Specifically include vec3 for 3D vectors

Object3D::Object3D() : vertexCounter(0), edgeCounter(0), loopCounter(0), faceCounter(0) {}

int Object3D::addVertex(float x, float y, float z) {
    int id = vertexCounter++;
    vertices.emplace(id, Vertex(id, x, y, z));
    return id;
}

Vertex& Object3D::getVertex(int id) {
    auto it = vertices.find(id);
    if (it == vertices.end())
        throw std::runtime_error("Vertex ID not found");
    return it->second;
}

const std::unordered_map<int, Vertex>& Object3D::getVertices() const {
    return vertices;
}

int Object3D::addEdge(int startVertexID, int endVertexID) {
    if (vertices.find(startVertexID) == vertices.end() || vertices.find(endVertexID) == vertices.end())
        throw std::runtime_error("One or both vertices not found");
    int id = edgeCounter++;
    edges.emplace(id, Edge(id, startVertexID, endVertexID));
    return id;
}

Edge& Object3D::getEdge(int id) {
    auto it = edges.find(id);
    if (it == edges.end())
        throw std::runtime_error("Edge ID not found");
    return it->second;
}

int Object3D::addLoop() {
    int id = loopCounter++;
    loops.emplace(id, Loop(id));
    return id;
}

void Object3D::addEdgeToLoop(int loopID, int edgeID) {
    if (loops.find(loopID) == loops.end())
        throw std::runtime_error("Loop ID not found");
    if (edges.find(edgeID) == edges.end())
        throw std::runtime_error("Edge ID not found");
    loops[loopID].addEdge(edgeID);
}

Loop& Object3D::getLoop(int id) {
    auto it = loops.find(id);
    if (it == loops.end())
        throw std::runtime_error("Loop ID not found");
    return it->second;
}

int Object3D::addFace() {
    int id = faceCounter++;
    faces.emplace(id, Face(id));
    return id;
}

void Object3D::addEdgeToFace(int faceID, int edgeID) {
    if (faces.find(faceID) == faces.end())
        throw std::runtime_error("Face ID not found");
    if (edges.find(edgeID) == edges.end())
        throw std::runtime_error("Edge ID not found");
    faces[faceID].addLoop(edgeID);  // This is where you should add loops to the face.
}

Face& Object3D::getFace(int id) {
    auto it = faces.find(id);
    if (it == faces.end())
        throw std::runtime_error("Face ID not found");
    return it->second;
}

void Object3D::loadFromFile(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    while (std::getline(infile, line)) {
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos)
            line = line.substr(0, commentPos);
        if (line.empty())
            continue;

        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            float x, y, z;
            iss >> x >> y >> z;
            addVertex(x, y, z);
        }
        else if (prefix == "e") {
            int start, end;
            iss >> start >> end;
            addEdge(start, end);
        }
        else if (prefix == "f") {
            std::vector<int> edgeIndices;
            int edgeID;
            while (iss >> edgeID) {
                edgeIndices.push_back(edgeID);
            }

            int faceID = addFace();
            for (int eid : edgeIndices) {
                addEdgeToFace(faceID, eid);
            }

            int loopID = addLoop();
            for (int eid : edgeIndices) {
                addEdgeToLoop(loopID, eid);
            }

            // Now add loop to the face
            faces[faceID].addLoop(loopID);
        }
    }

    infile.close();
}

void Object3D::initializeOpenGL() {
    std::vector<glm::vec3> vertexData;
    std::vector<unsigned int> indices;

    for (const auto& [id, vertex] : vertices) {
        vertexData.emplace_back(vertex.x, vertex.y, vertex.z);
    }

    for (const auto& [fid, face] : faces) {
        for (int eid : face.loopIDs) {
            Edge& edge = getEdge(eid);
            indices.push_back(edge.startVertex);
        }
    }

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(glm::vec3), vertexData.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void Object3D::draw() {
    glBindVertexArray(VAO);
    glDrawElements(GL_LINES, edges.size() * 2, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}
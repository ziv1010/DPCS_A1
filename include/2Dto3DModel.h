// 2Dto3DModel.h

#ifndef MODEL_3D_H
#define MODEL_3D_H

#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <glm/glm.hpp> // For GLM vector operations

#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "3Dto2DModel.h"

class Model_3D {
public:
    std::vector<Vertex> v;
    std::vector<Edge> e;
    std::vector<Surface> s;

    Model_2D topView;
    Model_2D frontView;
    Model_2D sideView;

    // Constructor
    Model_3D();

    // Function to generate the model from input file
    void generate_model(const std::string& inputFilePath);

    // Methods
    Model_2D basic_proj(int direction);
    void calc_2d(int direction);
    void calc_2d_1(int direction);
    void form_wireframe();
    void form_3dfaces();
    std::vector<int> get_neighbours(const Vertex& vn);
    Edge* is_edge(int a, int b);
    void model_3d_verification();
    void save2Dmodels(std::ofstream& out);
    void save3Dmodel(std::ofstream& out);
};

#endif // MODEL_3D_H
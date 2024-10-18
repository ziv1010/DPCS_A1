#ifndef ORTHOGRAPHICVIEW_H
#define ORTHOGRAPHICVIEW_H

#include "Object3D.h"
#include <GLFW/glfw3.h>  // Include GLFW and OpenGL

class OrthographicView {
public:
    void projectOntoXY(const Object3D& obj);
    void projectOntoXZ(const Object3D& obj);
    void projectOntoYZ(const Object3D& obj);
};

#endif
#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "Object3D.h"
#include "Vector3D.h"
#include <glm/glm.hpp>  // GLM for matrix transformations

class Transformation {
public:
    static void rotate(Object3D& obj, const Vector3D& axis, float angle);
    static void scale(Object3D& obj, float scaleFactor);
    static void translate(Object3D& obj, const Vector3D& translationVector);
};

#endif
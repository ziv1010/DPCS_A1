#ifndef GEOMETRICPROPERTIES_H
#define GEOMETRICPROPERTIES_H

#include "Object3D.h"

class GeometricProperties {
public:
    float calculateSurfaceArea(const Object3D& obj);
    float calculateVolume(const Object3D& obj);
};

#endif
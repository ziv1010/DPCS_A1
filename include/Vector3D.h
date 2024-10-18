#ifndef VECTOR3D_H
#define VECTOR3D_H

class Vector3D {
public:
    float x, y, z;

    Vector3D();
    Vector3D(float x, float y, float z);

    Vector3D operator+(const Vector3D& other) const;
    Vector3D operator*(float scalar) const;
};

#endif
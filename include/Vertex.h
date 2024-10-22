#ifndef VERTEX_H
#define VERTEX_H

class Vertex {
public:
    float x, y, z;

    Vertex(float x, float y, float z) : x(x), y(y), z(z) {}
    
    // Accessors
    float getX() const { return x; }
    float getY() const { return y; }
    float getZ() const { return z; }

    // Mutators
    void setX(float x) { this->x = x; }
    void setY(float y) { this->y = y; }
    void setZ(float z) { this->z = z; }
};

#endif // VERTEX_H
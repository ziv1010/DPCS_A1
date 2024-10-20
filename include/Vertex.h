// Vertex.h

#ifndef VERTEX_H
#define VERTEX_H

class Vertex {
public:
    // Attributes
    float x;
    float y;
    float z;
    bool isTrue;
    int vNo;

    // Constructors
    Vertex();
    Vertex(float x, float y, float z, int n);
    Vertex(const Vertex& v);

    // Methods
    Vertex getXY() const;
    Vertex getYZ() const;
    Vertex getXZ() const;
    bool same(const Vertex& a) const;
};

#endif // VERTEX_H
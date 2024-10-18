#ifndef VERTEX_H
#define VERTEX_H

struct Vertex {
    int id;         // Unique identifier
    float x, y, z;  // Coordinates

    Vertex(int id_, float x_, float y_, float z_)
        : id(id_), x(x_), y(y_), z(z_) {}
};

#endif // VERTEX_H
// Edge.h

#ifndef EDGE_H
#define EDGE_H

#include "Vertex.h"

class Edge {
public:
    // Attributes
    Vertex a;
    Vertex b;
    bool hidden;
    bool isTrue;
    int eno;

    // Constructors
    Edge();
    Edge(const Vertex& a, const Vertex& b, int c);

    // Methods
    Vertex getNeighbour(const Vertex& d) const;
    Vertex midpoint() const;
    Edge getEdge() const;
    bool overlap(const Edge& e2) const;
};

#endif // EDGE_H
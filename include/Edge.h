#ifndef EDGE_H
#define EDGE_H

struct Edge {
    int id;         // Unique identifier
    int startVertex; // ID of the start vertex
    int endVertex;   // ID of the end vertex

    Edge(int id_, int start, int end)
        : id(id_), startVertex(start), endVertex(end) {}
};

#endif // EDGE_H
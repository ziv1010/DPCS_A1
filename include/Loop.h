#ifndef LOOP_H
#define LOOP_H

#include <vector>

struct Loop {
    int id;                 // Unique identifier
    std::vector<int> edgeIDs; // IDs of edges forming the loop

    Loop(int id_) : id(id_) {}

    void addEdge(int edgeID) {
        edgeIDs.push_back(edgeID);
    }

    bool isClosed();
};

#endif // LOOP_H
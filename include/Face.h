#ifndef FACE_H
#define FACE_H

#include <vector>
#include "Loop.h"

struct Face {
    int id;                    // Unique identifier
    std::vector<int> loopIDs;  // IDs of loops forming the face

    Face(int id_) : id(id_) {}

    void addLoop(int loopID) {
        loopIDs.push_back(loopID);
    }
};

#endif // FACE_H
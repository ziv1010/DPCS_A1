#include "Transformation.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>  // For glm::rotate

void Transformation::rotate(Object3D& obj, const Vector3D& axis, float angle) {
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(angle), glm::vec3(axis.x, axis.y, axis.z));

    for (auto& vertex : obj.vertices) {
        glm::vec4 transformedVertex = rotationMatrix * glm::vec4(vertex.x, vertex.y, vertex.z, 1.0f);
        vertex = Vector3D(transformedVertex.x, transformedVertex.y, transformedVertex.z);
    }
}

void Transformation::scale(Object3D& obj, float scaleFactor) {
    for (auto& vertex : obj.vertices) {
        vertex.x *= scaleFactor;
        vertex.y *= scaleFactor;
        vertex.z *= scaleFactor;
    }
}

void Transformation::translate(Object3D& obj, const Vector3D& translationVector) {
    for (auto& vertex : obj.vertices) {
        vertex.x += translationVector.x;
        vertex.y += translationVector.y;
        vertex.z += translationVector.z;
    }
}
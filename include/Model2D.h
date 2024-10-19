// Model2D.h
#ifndef MODEL_2D_H
#define MODEL_2D_H

#include "model_2d_core.h"
#include "intersection.h"
#include "add_elements.h"
#include "display.h"
#include "geometry.h"

// Final Model2D class inheriting from all functional base classes
class Model2D : public Model2DIntersection,
               public Model2DAddElements,
               public Model2DDisplay,
               public Model2DGeometry {
public:
    // Default constructor
    Model2D()
        : Model2DIntersection(),
          Model2DAddElements(),
          Model2DDisplay(),
          Model2DGeometry() {}

    // Additional constructors or member functions can be added here if necessary
};

#endif // MODEL_2D_H
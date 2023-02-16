#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Shape.h"

using namespace std;

class Sphere: public Shape {   
    public:
        float radius;
        Eigen::Vector3f centre;

        Sphere(float x, float y, float z, float radius);
};
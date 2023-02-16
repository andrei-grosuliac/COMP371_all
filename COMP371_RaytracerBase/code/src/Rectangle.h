#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Shape.h"

using namespace std;

class Rectangle: public Shape {   
    public:
        Eigen::Vector3f p1, p2, p3, p4;

        Rectangle(float p1a, float p1b, float p1c, float p2a, float p2b, float p2c, float p3a, float p3b, float p3c, float p4a, float p4b, float p4c);
};
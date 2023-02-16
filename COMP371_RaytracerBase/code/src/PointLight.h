#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Light.h"

using namespace std;

class PointLight: public Light {
    public:
        Eigen::Vector3f centre;

        PointLight(float x, float y, float z);
};
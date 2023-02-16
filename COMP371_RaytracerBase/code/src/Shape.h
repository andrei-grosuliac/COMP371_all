#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

class Shape{
    public:
        float ka, kd, ks;
        Eigen::Vector3f ac, dc, sc;
        float pc;
        string type;

        void SetReflectionCoefficients(float ka, float kd, float ks);
        void SetReflectionColor(float ac1, float ac2, float ac3, float dc1, float dc2, float dc3, float sc1, float sc2, float sc3);
        void SetPhongCoefficients(float pc);
};
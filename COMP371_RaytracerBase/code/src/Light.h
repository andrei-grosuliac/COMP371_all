#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

class Light {
    public :
        Eigen::Vector3f id, is;

        void SetLightIntensity(float id1, float id2, float id3, float is1, float is2, float is3);
};
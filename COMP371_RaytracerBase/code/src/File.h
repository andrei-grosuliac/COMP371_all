#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

class File {
    public:       
        string filename;
        int size [2];
        int raysperpixel [2];
        bool globalIllumination;
        int maxBounces;
        float probTerminate;

        Vector3f lookat, up, center, ai, bkc;

        float fov, pc;
        File(){};
};
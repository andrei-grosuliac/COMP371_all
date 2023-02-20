#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using Vector3f = Eigen::Vector3f;

class Ray {
    public:
        Vector3f p;     // Point
        Vector3f d;     // Direction
        bool isNull = false;
        Ray() {isNull = true;}
        Ray(Vector3f point, Vector3f direction) {p = point; d = direction;}
        Ray(const Ray &ray) {p = ray.p; d = ray.d;}
        Vector3f at(float t) const {
            return p + t * d;
        }
};
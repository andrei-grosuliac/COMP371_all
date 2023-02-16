#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "Shape.h"
#include "Light.h"
#include "Ray.h"
#include "Sphere.h"
#include "Rectangle.h"

#include "../external/json.hpp"
#include "../external/simpleppm.h"

using json = nlohmann::json;
using Vector3f = Eigen::Vector3f;

class RayTracer {
    public:
        RayTracer(const json& scene);

        void run();
        void createGeometry(const json& scene);
        void createLights(const json& scene);
        void createOutputs(const json& scene);
        bool sphereIntersection(Ray ray, Sphere sphere, float& t);
        void rayCast(Ray ray);

        vector<Sphere> spheres;
        vector<Rectangle> rectangles;
        vector<Light> lights;
        string filename;
        int size [2];

        Vector3f lookat, up, center, ai, bkc;
        Vector3f color;

        float fov;        
};
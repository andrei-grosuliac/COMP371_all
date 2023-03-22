#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "Shape.h"
#include "Light.h"
#include "PointLight.h"
#include "AreaLight.h"
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
        bool rectangleIntersection(Ray ray, Rectangle rectangle, float& t);
        bool is_inside_triangle(const Vector3f& p, const Vector3f& v0, const Vector3f& v1, const Vector3f& v2);
        Vector3f rayCast(Ray ray);
        Vector3f compute_color(Shape shape, Vector3f intersection_pt, Vector3f normal, Ray ray);
        Vector3f clamp(Vector3f color, float min, float max);
        bool testShadow(Ray &ray, float max_distance);
        float random_float();
        Vector3f sample_hemisphere(const Vector3f& normal);
        bool intersect_scene(Ray ray, Shape& hit_shape, Vector3f& intersection_pt, Vector3f& normal);
        Vector3f path_trace(Ray& ray, int depth);
        Ray jitter_ray(const Ray& original_ray, float random_offset_x, float random_offset_y);
        void coordinate_system(Vector3f original_ray_d, Vector3f& u, Vector3f& v, Vector3f& w);

        vector<Sphere> spheres;
        vector<Rectangle> rectangles;
        vector<PointLight> pointLights;
        vector<AreaLight> areaLights;
        string filename;
        int size [2];
        int raysperpixel [2];
        bool globalIllumination;
        int maxBounces;
        float probTerminate;

        Vector3f lookat, up, center, ai, bkc;
        //Vector3f color;

        float fov, pc;        
};